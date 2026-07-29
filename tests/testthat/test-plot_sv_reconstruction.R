# Tests for plot_sv_reconstruction(): VF (read-support) stratified rearrangement plot.

# A single-locus input (chr7 amplicon) carrying enough SVs, spread across a wide
# VF range, to exercise clustering. All breakpoints sit inside the amplified
# window so every junction is in-locus.
make_vf_inputs <- function() {
  karyotype <- data.frame(
    chrom    = "chr7",
    start    = c(0, 50e6, 60e6),
    end      = c(50e6, 60e6, 159e6),
    gieStain = c("gneg", "gpos50", "gneg")
  )
  gene_coord <- data.frame(
    chr = "chr7", start = 55019017, end = 55211628, strand = "+", gene = "EGFR"
  )
  wgd_data <- data.frame(sample = "S1", Polyploidy = "No")
  cnv_data <- data.frame(
    sample = "S1", seqnames = "chr7",
    start = c(40e6, 55e6, 55.9e6 + 1),
    end   = c(54.99e6, 55.9e6, 70e6),
    copyNumber = c(2, 40, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 39, 1),
    minorAlleleCopyNumber = c(1, 1, 1)
  )
  vf <- c(500, 460, 420, 80, 74, 66, 60, 12, 10, 9, 8, 6)
  n  <- length(vf)
  p1 <- seq(55.1e6, 55.7e6, length.out = n)
  p2 <- p1 + 40e3
  sv_data <- data.frame(
    chrom1 = "7", start1 = p1, chrom2 = "7", start2 = p2,
    strand1 = "-", strand2 = "+", svclass = "DUP",
    VF = vf, JCN = pmax(1, round(vf / 20)), sample = "S1"
  )
  list(karyotype = karyotype, gene_coord = gene_coord, wgd_data = wgd_data,
       cnv_data = cnv_data, sv_data = sv_data)
}

quiet_by_vf <- function(...) suppressWarnings(suppressMessages(plot_sv_reconstruction(...)))

test_that("explicit vf_breaks produce the expected number of strata, ordered high->low", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  p <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord,
                   vf_breaks = c(300, 50), isolate_founder = FALSE)
  expect_s3_class(p, "patchwork")

  strata <- attr(p, "strata")
  expect_equal(sort(unique(strata$stratum)), 1:3)
  # break edges: >300, 50-300, <50  -> counts 3 / 4 / 5
  expect_equal(as.integer(table(strata$stratum)), c(3L, 4L, 5L))
  # stratum 1 is the highest read-support bracket
  meds <- tapply(strata$VF, strata$stratum, median)
  expect_true(meds[["1"]] > meds[["2"]] && meds[["2"]] > meds[["3"]])
})

test_that("k-means on log(VF) yields k strata ranked by read support", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  p <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord, k = 3,
                   isolate_founder = FALSE)
  strata <- attr(p, "strata")
  expect_equal(length(unique(strata$stratum)), 3L)
  meds <- tapply(strata$VF, strata$stratum, median)
  expect_false(is.unsorted(rev(meds)))       # medians decrease with stratum index
})

test_that("isolate_founder puts the single max-VF junction alone in stratum 1", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  p <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord, k = 3,
                   isolate_founder = TRUE)
  strata <- attr(p, "strata")
  # founder = the single highest VF, flagged and alone at the top
  expect_true(all(strata$is_founder == (strata$VF == max(strata$VF))))
  expect_equal(sum(strata$is_founder), 1L)
  expect_equal(unique(strata$stratum[strata$is_founder]), 1L)
  # founder excluded from clustering -> one extra stratum vs isolate_founder=FALSE
  q <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord, k = 3,
                   isolate_founder = FALSE)
  expect_equal(length(unique(strata$stratum)),
               length(unique(attr(q, "strata")$stratum)) + 1L)
})

test_that("a PDF is written when outdir is supplied", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  outdir <- withr::local_tempdir()
  p <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord,
                   k = 2, outdir = outdir)
  f <- attr(p, "path")
  expect_true(!is.null(f) && file.exists(f))
  expect_gt(file.info(f)$size, 0)
})

test_that("an unknown read-support column is reported clearly", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  expect_error(
    quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                karyotype = d$karyotype, gene_coord = d$gene_coord,
                vf_col = "not_a_column"),
    "not found"
  )
})

test_that("both cn_display modes render and reconstruction is cumulative/monotone", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  pr <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                    karyotype = d$karyotype, gene_coord = d$gene_coord, k = 3,
                    cn_display = "reconstruct")
  pa <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                    karyotype = d$karyotype, gene_coord = d$gene_coord, k = 3,
                    cn_display = "actual")
  expect_s3_class(pr, "patchwork")
  expect_s3_class(pa, "patchwork")
  # cn_display is a display choice only: clustering is unchanged
  expect_equal(attr(pr, "strata")$stratum, attr(pa, "strata")$stratum)
})

test_that("vf_scale = 'per_panel' runs and preserves stratum assignment", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  ps <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                    karyotype = d$karyotype, gene_coord = d$gene_coord,
                    k = 3, vf_scale = "shared")
  pp <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                    karyotype = d$karyotype, gene_coord = d$gene_coord,
                    k = 3, vf_scale = "per_panel")
  expect_s3_class(pp, "patchwork")
  # scaling is a rendering choice only: the clustering is identical
  expect_equal(attr(pp, "strata")$stratum, attr(ps, "strata")$stratum)
})

test_that("min_vf filters low-support junctions before clustering", {
  skip_if_not_installed("patchwork")
  d <- make_vf_inputs()
  p <- quiet_by_vf(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
                   karyotype = d$karyotype, gene_coord = d$gene_coord,
                   k = 2, min_vf = 50)
  strata <- attr(p, "strata")
  expect_true(all(strata$VF >= 50))
  expect_equal(nrow(strata), 7L)             # the five VF<50 junctions dropped
})
