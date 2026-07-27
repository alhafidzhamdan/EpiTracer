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
  d <- make_recon_inputs()
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
  d <- make_recon_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data, wgd_data = d$wgd_data,
      karyotype = d$karyotype, gene_coord = d$gene_coord
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("explicit loci strings are accepted", {
  d <- make_recon_inputs()
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
  d <- make_recon_inputs()   # amplified loci only, no homdel
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

test_that("auto-detection errors when no amplification is present", {
  d <- make_recon_inputs()
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
