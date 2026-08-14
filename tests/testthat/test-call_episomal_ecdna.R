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

test_that("differing boundary-DUP VFs do not error (length-1 highest-VF test)", {
  # Regression: give the two boundary-DUP breakends slightly different VFs, as
  # real read-support noise would. The highest still reaches the amplicon max,
  # so the amplicon is episomal — and the caller must not error on a length>1
  # condition.
  d <- make_episome_inputs(flank_cn = 2)
  d$breakpoints_gr$VF[d$breakpoints_gr$svclass == "DUP"] <- c(1000, 1100)
  res <- call_episomal_ecdna(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_s3_class(res, "data.table")
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary_has_highest_VF == "TRUE"))
  expect_true(all(res$episomal == "TRUE"))
})

test_that("detect_amplicon_seeds finds the amplified region from copy number", {
  d <- make_episome_inputs(flank_cn = 2)
  seeds <- detect_amplicon_seeds(d$cnv_gr)
  expect_s4_class(seeds, "GRanges")
  expect_gte(length(seeds), 1L)
  expect_false(is.null(seeds$ID))
  expect_false(is.null(seeds$WGS_ID))
  expect_true(all(seeds$WGS_ID == "S1"))

  # no amplification -> empty seeds
  d2 <- make_episome_inputs(flank_cn = 2)
  d2$cnv_gr$copyNumber <- 2
  expect_length(detect_amplicon_seeds(d2$cnv_gr), 0L)
})

test_that("ecdna_gr = NULL auto-detects seeds and still calls the episome", {
  d <- make_episome_inputs(flank_cn = 2)
  res <- call_episomal_ecdna(NULL, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_s3_class(res, "data.table")
  expect_true(any(res$episomal == "TRUE"))

  # errors clearly when there is nothing to seed from
  d0 <- make_episome_inputs(flank_cn = 2)
  d0$cnv_gr$copyNumber <- 2
  expect_error(
    call_episomal_ecdna(NULL, d0$breakpoints_gr, d0$cnv_gr, d0$cancer_genes_gr),
    "No focal amplicons"
  )
})
