## Flexible input readers: read_cnv / read_sv / read_purple_sv_vcf and the optional
## sample-name override that lets CN and SV files share a name despite differing
## file-naming conventions.

test_that("read_cnv reformats PURPLE-style CN, floors negatives, estimates ploidy", {
  cn <- data.frame(chromosome = c("chr1","chr1","chr9"),
                   start = c(1, 2000001, 1), end = c(2000000, 5000000, 3000000),
                   copyNumber = c(4, 8, -0.04),            # PURPLE emits occasional negatives
                   minorAlleleCopyNumber = c(2, 1, 0), majorAlleleCopyNumber = c(2, 7, 0))
  out <- read_cnv(cn, name = "DO1.purple.cnv.somatic.tsv")
  expect_true(all(c("sample","seqnames","start","end","copyNumber","ploidy",
                    "majorAlleleCopyNumber","minorAlleleCopyNumber") %in% names(out)))
  expect_equal(unique(out$sample), "DO1")                 # from the file name
  expect_true(all(out$copyNumber >= 0))                   # negatives clamped
  expect_true(is.finite(out$ploidy[1]) && out$ploidy[1] > 1.5)  # estimated, not defaulted
})

test_that("read_cnv sample and ploidy arguments override", {
  cn <- data.frame(seqnames = "chr7", start = 1, end = 1e6, copyNumber = 40,
                   minorAlleleCopyNumber = 1, majorAlleleCopyNumber = 39)
  out <- read_cnv(cn, sample = "TUMOUR_A", ploidy = 3.6, name = "whatever.tsv")
  expect_equal(unique(out$sample), "TUMOUR_A")
  expect_equal(unique(out$ploidy), 3.6)
})

test_that("read_sv keeps VF and JCN and normalises svclass", {
  sv <- data.frame(chrom1 = c("7","7"), start1 = c(54000000, 11352091),
                   chrom2 = c("7","13"), start2 = c(56000000, 50055509),
                   strand1 = c("-","-"), strand2 = c("+","+"),
                   svclass = c("DUP","TRA"), VF = c(1800, 33), JCN = c(14.6, 1.65),
                   name = c("a","b"), homlen = c(0, 0))
  out <- read_sv(sv, sample = "TUMOUR_A", name = "DO1.svs.bedpe")
  expect_equal(unique(out$sample), "TUMOUR_A")
  expect_true(all(c("VF","JCN") %in% names(out)))
  expect_equal(out$VF, c(1800, 33)); expect_equal(out$JCN, c(14.6, 1.65))
  expect_equal(out$svclass, c("DUP","TRA"))
  expect_true(all(grepl("^chr", out$chrom1)) && all(grepl("^chr", out$chrom2)))  # standardised 'chr' prefix
})

test_that("read_purple_sv_vcf extracts VF (INFO or FORMAT) and JCN, de-dups mates", {
  vcf <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR",
    "chr1\t8484797\tdup_o\tN\t]chr1:8484859]N\t9\tPASS\tPURPLE_JCN=1.47;VF=16;HOMLEN=1\tVF\t16",
    "chr1\t8484859\tdup_h\tN\tN[chr1:8484797[\t9\tPASS\tPURPLE_JCN=1.47;VF=16\tVF\t16",   # mate
    "chr2\t4884322\tdel_o\tN\tN[chr2:5108040[\t9\tPASS\tPURPLE_JCN=10.78\tVF\t273",        # VF only in FORMAT
    "chr1\t11352091\ttra_o\tN\t]chr13:50055509]N\t9\tPASS\tPURPLE_JCN=1.65;VF=33\tVF\t33")
  f <- tempfile(fileext = ".vcf"); writeLines(vcf, f)
  out <- read_purple_sv_vcf(f)
  expect_equal(nrow(out), 3)                              # mate de-duplicated
  expect_true(all(grepl("^chr", out$chrom1)))             # standardised 'chr' prefix
  dup <- out[out$chrom1 == "chr1" & out$start1 == 8484797, ]
  expect_equal(dup$svclass, "DUP"); expect_equal(dup$VF, 16); expect_equal(dup$JCN, 1.47)
  del <- out[out$svclass == "DEL", ]
  expect_equal(del$VF, 273)                               # pulled from the FORMAT column
  expect_true("TRA" %in% out$svclass)
})

test_that("prepare_amplicon_inputs aligns single-sample CN and SV with different names", {
  cn <- read_cnv(data.frame(seqnames = "chr7", start = c(1, 54000000, 56000001),
                            end = c(53999999, 56000000, 159000000), copyNumber = c(2, 42, 2),
                            minorAlleleCopyNumber = 1, majorAlleleCopyNumber = c(1, 41, 1)),
                 sample = "S")
  sv <- read_sv(data.frame(chrom1 = "7", start1 = 54000000, chrom2 = "7", start2 = 56000000,
                           strand1 = "-", strand2 = "+", svclass = "DUP", VF = 1800, JCN = 40),
                sample = "DIFFERENT_NAME")
  inp <- prepare_amplicon_inputs(cn, sv, sample = "S")    # SV name differs -> fallback to all rows
  expect_s4_class(inp$cnv_gr, "GRanges")
  expect_true(length(inp$breakpoints_gr) == 2)            # the DUP's two breakends
  expect_true("PURPLE_CN" %in% names(S4Vectors::mcols(inp$breakpoints_gr)))
})

test_that("gene_locus resolves a gene to a locus and fits to the CN event", {
  loc <- gene_locus("EGFR", genome = "hg38")
  expect_equal(nrow(loc), 1); expect_equal(loc$chr, "chr7")
  expect_true(loc$start < loc$end)
  # a homozygous deletion overlapping the gene grows the window to the deleted segment
  cn <- read_cnv(data.frame(seqnames = "chr9", start = c(1, 20000000, 23000001),
                            end = c(19999999, 23000000, 40000000), copyNumber = c(2, 0, 2),
                            minorAlleleCopyNumber = c(1, 0, 1), majorAlleleCopyNumber = c(1, 0, 1)),
                 sample = "S", ploidy = 2)
  hd <- gene_locus("CDKN2A", genome = "hg38", cnv = cn, sample = "S", target = "homdel")
  expect_true(hd$start <= 20000000 && hd$end >= 23000000)   # spans the deleted segment
  expect_warning(gene_locus(c("EGFR","NOTAGENE"), genome = "hg38"), "NOTAGENE")
})
