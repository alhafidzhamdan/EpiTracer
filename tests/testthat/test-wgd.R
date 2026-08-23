## call_wgd(): PURPLE whole-genome-doubling rule (major allele > 1.5 over >= half
## the bases on >= 11 of 22 autosomes).

## One full-length segment per autosome; the first `n_dup` autosomes carry a
## duplicated major allele (2 > 1.5), the rest do not (1 < 1.5).
make_wgd_cnv <- function(n_dup) {
  chrs <- paste0("chr", 1:22)
  GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = chrs, start = 1L, end = 1e8L, sample = "S1",
    copyNumber = 4, ploidy = 2,
    majorAlleleCopyNumber = c(rep(2, n_dup), rep(1, 22 - n_dup)),
    minorAlleleCopyNumber = 1, stringsAsFactors = FALSE), keep.extra.columns = TRUE)
}

test_that("WGD is called at 11 doubled autosomes but not 10", {
  expect_true(call_wgd(make_wgd_cnv(11))$wgd)
  expect_equal(call_wgd(make_wgd_cnv(11))$n_wgd_autosomes, 11L)
  expect_false(call_wgd(make_wgd_cnv(10))$wgd)
})

test_that("the half-the-bases threshold is respected", {
  ## chr1 duplicated over only 40% of its bases -> does not count
  gr <- GenomicRanges::makeGRangesFromDataFrame(data.frame(
    seqnames = c("chr1", "chr1", paste0("chr", 2:22)),
    start = c(1L, 40000001L, rep(1L, 21)),
    end   = c(40000000L, 100000000L, rep(1e8L, 21)),
    sample = "S1", copyNumber = 4, ploidy = 2,
    majorAlleleCopyNumber = c(2, 1, rep(2, 21)),   # chr1: only 40 Mb / 100 Mb doubled
    minorAlleleCopyNumber = 1, stringsAsFactors = FALSE), keep.extra.columns = TRUE)
  res <- call_wgd(gr)
  expect_equal(res$n_wgd_autosomes, 21L)           # chr1 excluded, chr2-22 counted
  expect_true(res$wgd)
})

test_that("call_wgd works generically on a data.frame with no sample column", {
  gr <- make_wgd_cnv(12)
  df <- as.data.frame(gr); df$sample <- NULL          # drop sample -> single sample
  res <- call_wgd(df)
  expect_equal(nrow(res), 1L)
  expect_true(res$wgd)
})

test_that("call_wgd errors without a major-allele column", {
  gr <- make_wgd_cnv(11); gr$majorAlleleCopyNumber <- NULL
  expect_error(call_wgd(gr), "majorAlleleCopyNumber")
})

test_that("call_wgd returns one row per sample", {
  g1 <- make_wgd_cnv(11); g1$sample <- "A"
  g2 <- make_wgd_cnv(3);  g2$sample <- "B"
  res <- call_wgd(c(g1, g2))
  expect_setequal(res$sample, c("A", "B"))
  expect_equal(res[sample == "A"]$wgd, TRUE)
  expect_equal(res[sample == "B"]$wgd, FALSE)
})
