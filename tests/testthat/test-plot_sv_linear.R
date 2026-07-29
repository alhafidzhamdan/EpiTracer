test_that("single focused locus (chromosome + range) returns a ggplot and writes files", {
  d <- make_plot_inputs()
  outdir <- withr::local_tempdir()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      outdir = outdir
    )
  ))
  expect_s3_class(p, "ggplot")
  out <- attr(p, "path")
  expect_length(out, 1)
  expect_match(out, "\\.pdf$")               # always PDF, never PNG
  expect_true(file.exists(out))
})

test_that("flank_pct widens the auto-detected window", {
  d <- make_multilocus_inputs()
  narrow <- suppressWarnings(suppressMessages(
    plot_sv_linear(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord,
      flank_pct = 0)))
  wide <- suppressWarnings(suppressMessages(
    plot_sv_linear(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord,
      flank_pct = 50)))
  # both build; the wider flank must not shrink the x-range
  rng_n <- diff(range(ggplot2::ggplot_build(narrow)$layout$panel_params[[1]]$x.range))
  rng_w <- diff(range(ggplot2::ggplot_build(wide)$layout$panel_params[[1]]$x.range))
  expect_gt(rng_w, rng_n)
})

test_that("plot builds without writing when outdir is NULL", {
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1)
    )
  ))
  expect_s3_class(p, "ggplot")
  expect_null(attr(p, "path"))
})

test_that("wgd_data is optional (title omits WGD status when not supplied)", {
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1)
    )
  ))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "S1")   # no "(Diploid)"/"(WGD)"
})

test_that("wgd_data annotates the title when supplied", {
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1)
    )
  ))
  expect_match(p$labels$title, "\\(Diploid\\)")
})

test_that("karyotype and gene_coord default to the bundled hg38 references", {
  d <- make_multilocus_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data)
  ))
  expect_s3_class(p, "ggplot")
})

test_that("wgd_data keyed by WGS_ID is accepted", {
  d <- make_plot_inputs()
  wgd_wgsid <- d$wgd_data
  names(wgd_wgsid)[names(wgd_wgsid) == "sample"] <- "WGS_ID"
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = wgd_wgsid,
      karyotype = d$karyotype, gene_coord = d$gene_coord, chromosome = "chr7"
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("auto-detects amplified loci across chromosomes and connects them", {
  d <- make_multilocus_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("explicit loci strings are accepted", {
  d <- make_multilocus_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      loci = c("chr7:54000000-57000000", "chr12:56000000-59000000")
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("events = 'homdel' errors when no homozygous deletion is present", {
  d <- make_multilocus_inputs()   # amplified loci only, no homdel
  expect_error(
    suppressWarnings(suppressMessages(
      plot_sv_linear(
        sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
        karyotype = d$karyotype, gene_coord = d$gene_coord, events = "homdel"
      )
    )),
    "homdel"
  )
})

test_that("displayExon draws exon models from cds_gr", {
  d <- make_plot_inputs()
  cds <- GenomicRanges::GRanges("chr7",
    IRanges::IRanges(c(55019017, 55142286), c(55019365, 55142437)), gene_name = "EGFR")
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      displayExon = TRUE, cds_gr = cds
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("displayExon = TRUE without cds_gr errors", {
  d <- make_plot_inputs()
  expect_error(
    suppressWarnings(suppressMessages(
      plot_sv_linear(
        sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
        karyotype = d$karyotype, gene_coord = d$gene_coord,
        chromosome = "chr7", displayExon = TRUE
      )
    )),
    "cds_gr"
  )
})

## Pull the SNV panel's plotted point data out of the stacked patchwork.
snv_points <- function(p) {
  snv_panel <- p[[2]]
  Reduce(function(acc, ly) {
    if (inherits(ly$geom, "GeomPoint")) rbind(acc, ly$data) else acc
  }, snv_panel$layers, init = NULL)
}

test_that("snv_data adds a stacked SNV panel (returns a patchwork)", {
  skip_if_not_installed("patchwork")
  d <- make_plot_inputs()
  outdir <- withr::local_tempdir()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      snv_data = make_snv_inputs(), outdir = outdir
    )
  ))
  expect_s3_class(p, "patchwork")            # stacked CN/SV + SNV panels
  out <- attr(p, "path")
  expect_match(out, "\\.pdf$")
  expect_true(file.exists(out))
})

test_that("without snv_data the result is a plain ggplot, not a patchwork", {
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1)
    )
  ))
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
})

test_that("intermutation-distance panel (default) is SNV-only, in-window, per-chr", {
  skip_if_not_installed("patchwork")
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      snv_data = make_snv_inputs()                       # snv_y defaults to "imd"
    )
  ))
  pts <- snv_points(p)
  ## 5 in-window S1 SNVs (the indel at 55.40 Mb and the OTHER sample are excluded);
  ## the first per chromosome (54.2 Mb) drops out as it has no previous SNV -> 4.
  expect_equal(nrow(pts), 4L)
  expect_true(all(is.finite(pts$yv)))
  ## y is log10(bp); the two 50 kb gaps (55.20->55.25, 55.25->55.30) are the min.
  expect_equal(min(pts$yv), log10(50000))
})

test_that("snv_y = 'vaf' plots VAF and drops out-of-range artefacts", {
  skip_if_not_installed("patchwork")
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      snv_data = make_snv_inputs(), snv_y = "vaf"
    )
  ))
  pts <- snv_points(p)
  ## 5 in-window S1 SNVs minus the artefactual VAF = 5 (55.25 Mb) -> 4, all in [0,1].
  expect_equal(nrow(pts), 4L)
  expect_true(all(pts$yv >= 0 & pts$yv <= 1))
})

test_that("snv_y = 'cn' plots the SNV copy number for every in-window SNV", {
  skip_if_not_installed("patchwork")
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      chromosome = "chr7", chromosome_range = matrix(c(54e6, 56e6), nrow = 1),
      snv_data = make_snv_inputs(), snv_y = "cn"
    )
  ))
  pts <- snv_points(p)
  ## All 5 in-window S1 SNVs are kept (no VAF/first-per-chr drop); the indel and
  ## the OTHER sample are excluded. y is variant_cn, whose in-window max is 8.
  expect_equal(nrow(pts), 5L)
  expect_equal(max(pts$yv), 8)
})

test_that("auto-detection errors when no amplification is present", {
  d <- make_multilocus_inputs()
  d$cnv_data$copyNumber <- 2
  d$cnv_data$majorAlleleCopyNumber <- 1
  expect_error(
    suppressWarnings(suppressMessages(
      plot_sv_linear(
        sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
        karyotype = d$karyotype, gene_coord = d$gene_coord
      )
    )),
    "amplified"
  )
})
