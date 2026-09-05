## Chimeric amplicon (cross-chromosome fusion): the standalone call_chimeric_amplicon()
## caller, and its independence from the episomal call (the combined `mechanism` label is
## assembled downstream by joining the caller outputs, so it is not tested here).

mn_inputs <- function(partner_cn = 50, boundary_dup = FALSE) {
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = c("chr7", "chr12", "chr7", "chr7"),
    start = c(55200000, 57200000, 55100000, 55400000),
    end   = c(55200000, 57200000, 55100000, 55400000),
    WGS_ID = "S1", event = c("TRA1", "TRA1", "D1", "D1"),
    svclass = c("TRA", "TRA", "DEL", "DEL"),
    PURPLE_AF = 0.9, PURPLE_JCN = c(50, 50, 10, 10), VF = c(3000, 3000, 100, 100),
    PURPLE_CN = c(50, partner_cn, 50, 50), insLen = 0L, HOMLEN = 0L)
  if (boundary_dup)  # add a boundary DUP spanning both edges -> also an episome
    bp <- rbind(bp, data.table::data.table(
      seqnames = "chr7", start = c(55000000, 55500000), end = c(55000000, 55500000),
      WGS_ID = "S1", event = "DUP1", svclass = "DUP", PURPLE_AF = 0.9, PURPLE_JCN = 40,
      VF = 2000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L))
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  list(ecdna_gr = ecdna_gr,
       breakpoints_gr = GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE),
       cnv_gr = GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE), cgg = cgg)
}

test_that("call_chimeric_amplicon flags a high-VF cross-chromosome amplified TRA", {
  d <- mn_inputs(partner_cn = 50)
  r <- call_chimeric_amplicon(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg, ext = 1e7, mc.cores = 1)
  expect_true(all(r$chimeric == "TRUE"))
})

test_that("call_chimeric_amplicon does NOT flag a TRA to a non-amplified partner", {
  d <- mn_inputs(partner_cn = 2)   # chr12 partner diploid
  r <- call_chimeric_amplicon(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg, ext = 1e7, mc.cores = 1)
  expect_true(all(r$chimeric == "FALSE"))
})

test_that("chimeric amplicon and episomal are independent (an amplicon can be both)", {
  d <- mn_inputs(partner_cn = 50, boundary_dup = TRUE)   # boundary DUP + cross-chr amplified TRA
  mn  <- call_chimeric_amplicon(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg, ext = 1e7, mc.cores = 1)
  epi <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg, ext = 1e7, mc.cores = 1)
  expect_true(all(mn$chimeric == "TRUE"))
  expect_true(all(epi$episomal == "TRUE"))
})
