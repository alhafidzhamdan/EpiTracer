## Translocation-bridge amplification: the standalone call_translocation_bridge_amp()
## caller (Lee et al., Nature 2023). TBA = an amplicon whose boundary is an
## amplified inter-chromosomal translocation. tb_confident additionally requires
## the asymmetric bridge-arm LOH footprint (amplicon and/or partner arm).

## An amplicon on chr8 (55.0-55.5 Mb) with an amplified TRA at its distal edge
## connecting to chr17, optionally plus a short fold-back inversion at the
## proximal edge.
tba_inputs <- function(with_fbi = FALSE, partner_cn = 50) {
  ecdna_gr <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  rows <- data.table::data.table(
    seqnames = c("chr8", "chr17"),
    start = c(55500000, 38000000), end = c(55500000, 38000000),
    WGS_ID = "S1", event = "TRA1", svclass = "TRA",
    PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 800,
    PURPLE_CN = c(50, partner_cn), insLen = 0L, HOMLEN = 0L)
  if (with_fbi)  # short tail-to-tail fold-back at the proximal edge (2 kb span)
    rows <- rbind(rows, data.table::data.table(
      seqnames = "chr8", start = c(55000000, 55002000), end = c(55000000, 55002000),
      WGS_ID = "S1", event = "FBI1", svclass = "t2tINV",
      PURPLE_AF = 0.9, PURPLE_JCN = 30, VF = 600, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L))
  bp <- GenomicRanges::makeGRangesFromDataFrame(rows, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr8", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 140000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55100000), gene = "MYC")
  list(ecdna_gr = ecdna_gr, breakpoints_gr = bp, cnv_gr = cnv_gr, cgg = cgg)
}

test_that("call_translocation_bridge_amp flags an amplified boundary translocation", {
  d <- tba_inputs()
  r <- call_translocation_bridge_amp(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg,
                                     ext = 1e7, mc.cores = 1)
  expect_true(all(r$tba == "TRUE"))
  expect_true(all(r$n_boundary_tra >= 1))
  expect_true(all(grepl("chr17", r$tb_partner_chr)))
})

test_that("no boundary translocation -> tba FALSE (a clean episome is not TBA)", {
  d <- make_episome_inputs(flank_cn = 2)   # boundary DUP + DEL, no TRA
  r <- call_translocation_bridge_amp(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                                     d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(r$tba == "FALSE"))
})

test_that("TBA is independent of the episomal call (an amplicon can be both)", {
  # boundary DUP episome that ALSO carries an amplified boundary translocation
  d <- make_episome_inputs(flank_cn = 2)
  extra <- GenomicRanges::GRanges(c("chr7", "chr12"),
    IRanges::IRanges(c(55500000, 60000000), width = 1),
    WGS_ID = "S1", event = "TRAx", svclass = "TRA",
    PURPLE_AF = 0.9, PURPLE_JCN = 35, VF = 700, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  bp <- c(d$breakpoints_gr, extra)
  tba <- call_translocation_bridge_amp(d$ecdna_gr, bp, d$cnv_gr, d$cancer_genes_gr,
                                       ext = 1e7, mc.cores = 1)
  epi <- call_simple_excision(d$ecdna_gr, bp, d$cnv_gr, d$cancer_genes_gr,
                              ext = 1e7, mc.cores = 1)
  expect_true(all(tba$tba == "TRUE"))
  expect_true(all(epi$episomal == "TRUE"))
})

## Confirmatory footprint (Lee et al.): bridge-arm subtelomeric LOH + non-bridge
## arm spared -> tb_confident. Build a chr8 q-arm amplicon with a large
## subtelomeric LOH block on the q (bridge) arm; vary the p (non-bridge) arm.
tba_loh_inputs <- function(pattern = c("asymmetric", "symmetric")) {
  pattern <- match.arg(pattern)
  ecdna_gr <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = c("chr8", "chr17"), start = c(55500000, 38000000), end = c(55500000, 38000000),
    WGS_ID = "S1", event = "TRA1", svclass = "TRA", PURPLE_AF = 0.9, PURPLE_JCN = 40,
    VF = 800, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L), keep.extra.columns = TRUE)
  # p arm (non-bridge) 0-44 Mb: heterozygous by default; q arm has amplicon + a
  # 120-145 Mb subtelomeric LOH block reaching the q-telomere.
  p_minor <- if (pattern == "symmetric") 0 else 1   # symmetric -> p-arm also LOH
  cnv <- data.table::data.table(
    seqnames = "chr8",
    start = c(1,        45000000, 55000000, 55500001, 120000000),
    end   = c(44000000, 54999999, 55500000, 119999999, 145138636),
    sample = "S1",
    copyNumber = c(2, 2, 50, 2, 2), ploidy = 2,
    majorAlleleCopyNumber = c(2, 1, 49, 1, 2),
    minorAlleleCopyNumber = c(p_minor, 1, 1, 1, 0))
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55100000), gene = "MYC")
  list(ecdna_gr = ecdna_gr, breakpoints_gr = bp, cnv_gr = cnv_gr, cgg = cgg)
}

test_that("tb_confident TRUE for asymmetric bridge-arm LOH (non-bridge spared)", {
  d <- tba_loh_inputs("asymmetric")
  r <- call_translocation_bridge_amp(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg,
        ext = 1e7, centromeres = load_centromeres("hg38"),
        chrom_lengths = load_chrom_lengths("hg38"), mc.cores = 1)
  expect_true(all(r$tba == "TRUE"))
  expect_true(all(r$tb_bridge_arm_loh == "TRUE"))
  expect_true(all(r$tb_nonbridge_spared == "TRUE"))
  expect_true(all(r$tb_confident == "TRUE"))
})

test_that("tb_confident FALSE when both arms show LOH (symmetric = chromothripsis-like)", {
  d <- tba_loh_inputs("symmetric")
  r <- call_translocation_bridge_amp(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg,
        ext = 1e7, centromeres = load_centromeres("hg38"),
        chrom_lengths = load_chrom_lengths("hg38"), mc.cores = 1)
  expect_true(all(r$tba == "TRUE"))          # still a TBA (boundary translocation)
  expect_true(all(r$tb_nonbridge_spared == "FALSE"))
  expect_true(all(r$tb_confident == "FALSE"))
})

test_that("confidence columns are NA without centromeres/chrom_lengths", {
  d <- tba_loh_inputs("asymmetric")
  r <- call_translocation_bridge_amp(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg,
                                     ext = 1e7, mc.cores = 1)
  expect_true(all(r$tba == "TRUE"))
  expect_true(all(r$tb_confident == "FALSE"))
  expect_true(all(r$tb_bridge_arm_loh == "NA"))
})

test_that("tb_confident TRUE via the translocated PARTNER arm LOH (amplicon arm clean)", {
  # chr8 q-arm amplicon with a boundary TRA to chr17; chr8 carries NO LOH, but the
  # partner (chr17) bridge arm shows subtelomeric LOH -> confident via partner arm.
  ecdna_gr <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = c("chr8", "chr17"), start = c(55500000, 38000000), end = c(55500000, 38000000),
    WGS_ID = "S1", event = "TRA1", svclass = "TRA", PURPLE_AF = 0.9, PURPLE_JCN = 40,
    VF = 800, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L), keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = c("chr8","chr8","chr8","chr8","chr17","chr17"),
    start = c(1, 45000000, 55000000, 55500001, 1,        60000000),
    end   = c(44000000, 54999999, 55500000, 145138636, 59999999, 83257441),
    sample = "S1",
    copyNumber = c(2,2,50,2, 2,2), ploidy = 2,
    majorAlleleCopyNumber = c(1,1,49,1, 1,2),
    minorAlleleCopyNumber = c(1,1,1,1,  1,0))   # chr17 q-telomeric block is LOH
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr8", IRanges::IRanges(55000000, 55100000), gene = "MYC")
  r <- call_translocation_bridge_amp(ecdna_gr, bp, cnv_gr, cgg, ext = 1e7,
        centromeres = load_centromeres("hg38"), chrom_lengths = load_chrom_lengths("hg38"), mc.cores = 1)
  expect_true(all(r$tba == "TRUE"))
  expect_true(all(r$tb_bridge_arm_loh == "FALSE"))   # amplicon arm clean
  expect_true(all(r$tb_partner_arm_loh == "TRUE"))   # partner (chr17) arm LOH
  expect_true(all(r$tb_confident == "TRUE"))
})
