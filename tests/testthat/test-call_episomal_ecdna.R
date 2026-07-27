test_that("input validation rejects non-GRanges and missing metadata", {
  gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(1, 2))
  expect_error(call_episomal_ecdna("nope", gr, gr, gr))
  # ecdna_gr without ID / WGS_ID
  expect_error(call_episomal_ecdna(gr, gr, gr, gr))
})

test_that("a textbook episome is classified as episomal with an excision scar", {
  d <- make_episome_inputs(flank_cn = 2)   # non-gained flanks
  res <- call_episomal_ecdna(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, ext = 1e7, mc.cores = 1)

  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 4L)

  # boundary DUP detected on both breakends and carries the highest VF
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary == "TRUE"))
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary_has_highest_VF == "TRUE"))

  # amplicon called episomal, with an excision scar
  expect_true(all(res$episomal == "TRUE"))
  expect_true(unique(res$has_excision_scar) == "TRUE")
})

test_that("gained flanks are NOT called episomal (validates line-199 fix)", {
  # Same boundary DUP + scar, but the flanking segments are amplified.
  # The original script's tautological flank test would still call this an
  # episome; the corrected `before & after` logic must not.
  d <- make_episome_inputs(flank_cn = 10)  # gained flanks
  res <- call_episomal_ecdna(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, ext = 1e7, mc.cores = 1)

  expect_true(all(res$episome_region == "FALSE"))
  expect_true(all(res$episomal == "FALSE"))
})
