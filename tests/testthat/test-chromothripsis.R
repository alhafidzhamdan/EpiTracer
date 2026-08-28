## Chromothripsis within an (episomal) amplicon: ShatterSeek-style hallmarks
## scored over the amplified footprint -- prevalence of internal SVs, randomness
## of the four junction orientations, and oscillating copy number. A clean simple
## episome must be negative; a shattered footprint with balanced orientations and
## an oscillating CN profile must be positive; a fold-back-dominated (BFB-like)
## footprint must be negative on the random-joins test.

## Amplicon footprint chr7 [55.0, 55.5] Mb (the make_episome_inputs locus).
## - boundary DUP = the circularisation junction at lo/hi
## - `classes` = the internal-SV orientation classes to scatter across the footprint
##   (repeated `n_each` times); balanced 4-class -> random joins, single-class -> biased
## - `oscillate` = CN profile alternates between two high states (turning points)
##   vs a flat plateau
chromo_inputs <- function(classes = c("DEL", "DUP", "h2hINV", "t2tINV"),
                          n_each = 3L, oscillate = TRUE) {
  lo <- 55000000L; hi <- 55500000L
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(lo, hi),
                                     ID = "S1_amp1", WGS_ID = "S1")

  ## boundary circularisation DUP (two breakends at lo and hi)
  rows <- data.table::data.table(
    seqnames = "chr7", start = c(lo, hi), end = c(lo, hi), WGS_ID = "S1",
    event = "bDUP", svclass = "DUP", PURPLE_AF = 0.9, PURPLE_JCN = 40,
    VF = 1000, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L)

  ## internal SVs: each event gets two breakends inside the footprint
  span <- hi - lo; k <- 0L
  for (cls in classes) for (i in seq_len(n_each)) {
    k <- k + 1L
    a <- lo + as.integer(span * (k / (length(classes) * n_each + 1)))
    b <- a + 30000L
    rows <- rbind(rows, data.table::data.table(
      seqnames = "chr7", start = c(a, b), end = c(a, b), WGS_ID = "S1",
      event = paste0("INT", k), svclass = cls, PURPLE_AF = 0.9, PURPLE_JCN = 20,
      VF = 500, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L))
  }
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(rows, keep.extra.columns = TRUE)

  ## CN segments: diploid flanks + amplicon footprint (oscillating or flat)
  if (oscillate) {
    fs <- c(lo, lo + 100000L, lo + 200000L, lo + 300000L, lo + 400000L)
    fe <- c(fs[-1] - 1L, hi)
    fcn <- c(40, 18, 42, 16, 40)                 # 4 turning points
  } else {
    fs <- lo; fe <- hi; fcn <- 40
  }
  cnv <- data.table::data.table(
    seqnames = "chr7",
    start = c(1L, 40000000L, fs, hi + 1L),
    end   = c(39999999L, lo - 1L, fe, 159000000L),
    sample = "S1", copyNumber = c(2, 2, fcn, 2), ploidy = 2)
  cnv$majorAlleleCopyNumber <- pmax(1, cnv$copyNumber - 1)
  cnv$minorAlleleCopyNumber <- 1
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

  cancer_genes_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628),
                                            gene = "EGFR")
  list(ecdna_gr = ecdna_gr, breakpoints_gr = breakpoints_gr,
       cnv_gr = cnv_gr, cancer_genes_gr = cancer_genes_gr)
}

run_ct <- function(d, ...) call_chromothripsis(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                                               d$cancer_genes_gr, mc.cores = 1, ...)

test_that("a shattered footprint (balanced joins + oscillating CN) is high-confidence chromothripsis", {
  r <- run_ct(chromo_inputs())
  expect_true(all(r$chromothripsis == "TRUE"))
  expect_true(all(r$chromothripsis_conf == "high"))
  expect_true(all(r$n_internal_sv >= 6))
  expect_true(all(r$cn_oscillations >= 3))
})

test_that("a clean simple episome (only a boundary junction) is NOT chromothripsis", {
  r <- call_chromothripsis(make_episome_inputs()$ecdna_gr,
                           make_episome_inputs()$breakpoints_gr,
                           make_episome_inputs()$cnv_gr,
                           make_episome_inputs()$cancer_genes_gr, mc.cores = 1)
  expect_true(all(r$chromothripsis == "FALSE"))
  expect_true(all(r$n_internal_sv < 6))
})

test_that("a fold-back-dominated (BFB-like) footprint fails the random-joins test", {
  r <- run_ct(chromo_inputs(classes = "h2hINV", n_each = 12L))
  expect_true(all(r$chromothripsis == "FALSE"))
  ## many internal SVs, but a single dominant orientation -> joins test rejects
  expect_true(all(r$n_internal_sv >= 6))
})

test_that("balanced joins but a flat CN plateau is only low-confidence", {
  r <- run_ct(chromo_inputs(oscillate = FALSE))
  expect_true(all(r$chromothripsis == "TRUE"))
  expect_true(all(r$chromothripsis_conf == "low"))
  expect_true(all(r$cn_oscillations < 3))
})

test_that("too few internal SVs (below min_sv) is NOT chromothripsis", {
  r <- run_ct(chromo_inputs(n_each = 1L))   # 4 internal + boundary = 5 events < 6
  expect_true(all(r$chromothripsis == "FALSE"))
})
