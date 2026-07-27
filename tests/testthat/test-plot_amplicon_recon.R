test_that("plot_amplicon_recon auto-detects loci and returns a ggplot", {
  d <- make_recon_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_amplicon_recon(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord
    )
  ))
  expect_s3_class(p, "ggplot")
  expect_null(attr(p, "path"))
})

test_that("plot_amplicon_recon writes PDF and PNG when given outdir", {
  d <- make_recon_inputs()
  outdir <- withr::local_tempdir()
  p <- suppressWarnings(suppressMessages(
    plot_amplicon_recon(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord,
      outdir = outdir, dpi = 150
    )
  ))
  out <- attr(p, "path")
  expect_true(any(grepl("\\.pdf$", out)))
  expect_true(any(grepl("\\.png$", out)))
  expect_true(all(file.exists(out)))
})

test_that("plot_amplicon_recon accepts explicit loci strings", {
  d <- make_recon_inputs()
  p <- suppressWarnings(suppressMessages(
    plot_amplicon_recon(
      sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
      wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord,
      loci = c("chr7:54000000-57000000", "chr12:56000000-59000000")
    )
  ))
  expect_s3_class(p, "ggplot")
})

test_that("plot_amplicon_recon errors when no amplification is present", {
  d <- make_recon_inputs()
  d$cnv_data$copyNumber <- 2       # nothing above 3 x ploidy
  d$cnv_data$majorAlleleCopyNumber <- 1
  expect_error(
    suppressWarnings(suppressMessages(
      plot_amplicon_recon(
        sample = "S1", cnv_data = d$cnv_data, sv_data = d$sv_data,
        wgd_data = d$wgd_data, karyotype = d$karyotype, gene_coord = d$gene_coord
      )
    )),
    "amplified"
  )
})
