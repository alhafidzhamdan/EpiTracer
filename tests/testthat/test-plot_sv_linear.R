test_that("plot_sv_linear writes a non-empty PDF for a simple locus", {
  d <- make_plot_inputs()
  outdir <- withr::local_tempdir()

  out <- suppressWarnings(suppressMessages(
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

  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 1000)
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
