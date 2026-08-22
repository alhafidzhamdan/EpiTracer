## Breakage-replication/fusion: the standalone call_brf() caller (adjacent
## parallel breakpoints), independent of the episomal call.

test_that("call_brf flags adjacent parallel breakpoints", {
  # Two same-orientation breakends within 20 kb, from DISTINCT junctions, are
  # "adjacent parallel breakpoints" -> brf = TRUE (Mendez-Dorantes/Zhang/Pellman,
  # Nat Genet 2026). Orientation from the per-breakend strand when supplied.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7", start = c(55100000, 55110000, 55400000, 55600000),
    end   = c(55100000, 55110000, 55400000, 55600000),
    WGS_ID = "S1", event = c("J1", "J2", "J1", "J2"), svclass = "h2hINV",
    bp_strand = c("+", "+", "+", "+"),
    PURPLE_AF = 0.9, PURPLE_JCN = 30, VF = 500, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  r <- call_brf(ecdna_gr, breakpoints_gr, cnv_gr, cgg, ext = 1e7, mc.cores = 1)
  expect_true(all(r$brf == "TRUE"))              # J1@55.100 and J2@55.110 are parallel (10 kb, same +)
  expect_true(all(r$n_parallel_pairs >= 1))
})

test_that("call_brf is independent of episomal (an episome can also be BRF)", {
  # a clean boundary-DUP episome PLUS adjacent parallel breakpoints
  d <- make_episome_inputs(flank_cn = 2)
  extra <- GenomicRanges::GRanges("chr7", IRanges::IRanges(c(55150000, 55160000), width = 1),
    WGS_ID = "S1", event = c("P1", "P2"), svclass = "h2hINV",
    PURPLE_AF = 0.9, PURPLE_JCN = 20, VF = 100, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  d$breakpoints_gr <- c(d$breakpoints_gr, extra)
  brf <- call_brf(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  epi <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(brf$brf == "TRUE"))
  expect_true(all(epi$episomal == "TRUE"))
})

test_that("call_brf returns brf FALSE with no adjacent parallel breakpoints", {
  d <- make_episome_inputs(flank_cn = 2)   # only a boundary DUP + DEL, no inversions
  r <- call_brf(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(r$brf == "FALSE"))
})
