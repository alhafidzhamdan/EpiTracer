## plot_sv_circos(): whole-genome per-sample circos (requires circlize).

test_that("plot_sv_circos writes a PDF for a sample without error", {
  skip_if_not_installed("circlize")
  inp <- make_multilocus_inputs()
  td <- withr::local_tempdir()
  expect_no_error(
    plot_sv_circos("S1", inp$sv_data, inp$cnv_data, inp$wgd_data, outdir = td)
  )
  expect_true(file.exists(file.path(td, "S1_circos.pdf")))
})

test_that("plot_sv_circos errors clearly for an absent sample", {
  skip_if_not_installed("circlize")
  inp <- make_multilocus_inputs()
  expect_error(
    plot_sv_circos("NOPE", inp$sv_data, inp$cnv_data, inp$wgd_data, outdir = withr::local_tempdir()),
    "No structural variants"
  )
})

test_that("plot_sv_circos accepts highlight_events (by name)", {
  skip_if_not_installed("circlize")
  inp <- make_multilocus_inputs()
  sv <- inp$sv_data
  sv$name <- paste0("SV", seq_len(nrow(sv)))          # give the SVs ids to match
  td <- withr::local_tempdir()
  expect_no_error(
    plot_sv_circos("S1", sv, inp$cnv_data, inp$wgd_data, outdir = td,
                   highlight_events = "SV1", dim_unhighlighted = TRUE)
  )
  expect_true(file.exists(file.path(td, "S1_circos.pdf")))
})

test_that("highlight_events errors when no id column is present", {
  skip_if_not_installed("circlize")
  inp <- make_multilocus_inputs()                     # sv_data has no name/event col
  expect_error(
    plot_sv_circos("S1", inp$sv_data, inp$cnv_data, inp$wgd_data,
                   outdir = withr::local_tempdir(), highlight_events = "x"),
    "id column"
  )
})

test_that("plot_sv_circos respects overwrite = FALSE", {
  skip_if_not_installed("circlize")
  inp <- make_multilocus_inputs()
  td <- withr::local_tempdir()
  plot_sv_circos("S1", inp$sv_data, inp$cnv_data, inp$wgd_data, outdir = td)
  f <- file.path(td, "S1_circos.pdf")
  mtime <- file.info(f)$mtime
  plot_sv_circos("S1", inp$sv_data, inp$cnv_data, inp$wgd_data, outdir = td, overwrite = FALSE)
  expect_identical(file.info(f)$mtime, mtime)   # untouched
})
