test_that("plot_sv_linear returns a ggplot and writes a non-empty PDF", {
  d <- make_plot_inputs()
  outdir <- withr::local_tempdir()

  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample     = "S1",
      chromosome = "chr7",
      karyotype  = d$karyotype,
      gene_coord = d$gene_coord,
      wgd_data   = d$wgd_data,
      cnv_data   = d$cnv_data,
      sv_data    = d$sv_data,
      outdir     = outdir
    )
  ))

  expect_s3_class(p, "ggplot")
  out <- attr(p, "path")
  expect_true(all(file.exists(out)))
  expect_true(any(grepl("\\.pdf$", out)))
  expect_true(any(grepl("\\.png$", out)))
  expect_true(all(file.info(out)$size > 1000))
})

test_that("plot_sv_linear builds a plot without writing when outdir is NULL", {
  d <- make_plot_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_sv_linear(
      sample = "S1", chromosome = "chr7",
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      wgd_data = d$wgd_data, cnv_data = d$cnv_data, sv_data = d$sv_data
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
      sample = "S1", chromosome = "chr7",
      karyotype = d$karyotype, gene_coord = d$gene_coord,
      wgd_data = wgd_wgsid, cnv_data = d$cnv_data, sv_data = d$sv_data
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("displayExon = TRUE requires cds_gr", {
  d <- make_plot_inputs()
  outdir <- withr::local_tempdir()

  expect_error(
    suppressWarnings(suppressMessages(
      plot_sv_linear(
        sample = "S1", chromosome = "chr7",
        karyotype = d$karyotype, gene_coord = d$gene_coord,
        wgd_data = d$wgd_data, cnv_data = d$cnv_data, sv_data = d$sv_data,
        outdir = outdir, displayExon = TRUE
      )
    )),
    "cds_gr"
  )
})
