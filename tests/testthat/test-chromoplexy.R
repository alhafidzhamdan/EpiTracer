## Chromoplexy caller: closed balanced-rearrangement cycles (Baca et al. 2013).

test_that("a closed 3-chromosome cycle is called as one chromoplexy", {
  inp <- make_chain_inputs(make_chromoplexy_junctions(close = TRUE), cn = 2)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)

  expect_equal(nrow(res), 1L)
  expect_equal(res$topology, "cycle")
  expect_equal(res$n_junctions, 3L)
  expect_equal(res$n_bridges, 3L)          # 3 deletion bridges (>= min_pairs)
  expect_equal(res$n_chromosomes, 3L)
  expect_equal(res$n_interchrom, 3L)
  expect_equal(res$frac_cp, 1)             # no intruding junctions
  expect_true(grepl("chr1", res$chromosomes) && grepl("chr12", res$chromosomes))
})

test_that("an OPEN chain (ends not rejoined) is NOT called", {
  inp <- make_chain_inputs(make_chromoplexy_junctions(close = FALSE), cn = 2)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)
  expect_equal(nrow(res), 0L)
})

test_that("an AMPLIFIED cycle is rejected by the low-copy filter", {
  ## same closed topology, but every junction sits at 50 copies (ploidy 2) -> not
  ## a balanced chromoplexy; it is an amplicon and must not be called.
  inp <- make_chain_inputs(make_chromoplexy_junctions(close = TRUE), cn = 50)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)
  expect_equal(nrow(res), 0L)
})

test_that("high copy is tolerated on a polysomic (high-ploidy) baseline", {
  ## cn = 5 exceeds the diploid cutoff (3) but not the ploidy-4 cutoff (max_cn *
  ## 4/2 = 6): a balanced chain on a polysomic chromosome must still be called.
  inp <- make_chain_inputs(make_chromoplexy_junctions(close = TRUE), cn = 5, ploidy = 4)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)
  expect_equal(nrow(res), 1L)
})

test_that("fewer than min_pairs junctions is not a chromoplexy", {
  jun <- make_chromoplexy_junctions(close = TRUE)[1:2, ]   # only 2 junctions
  inp <- make_chain_inputs(jun, cn = 2)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)
  expect_equal(nrow(res), 0L)
})

test_that("an intruding junction lowers the cleanliness score (frac_cp < 1)", {
  jun <- make_chromoplexy_junctions(close = TRUE)
  ## add a short balanced intrachromosomal DEL whose near end lands inside the
  ## chr1 bridge locus but which is too short (1 Mb < min_span) to be a chain edge
  jun <- rbind(jun, data.frame(chr1 = "chr1", pos1 = 10003000,
                               chr2 = "chr1", pos2 = 11003000, stringsAsFactors = FALSE))
  inp <- make_chain_inputs(jun, cn = 2)
  res <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)

  expect_equal(nrow(res), 1L)             # still a clean 3-junction cycle
  expect_equal(res$n_junctions, 3L)       # the intruder is not a cycle edge
  expect_equal(res$num_other, 1L)         # its near breakend intrudes on a locus
  expect_equal(res$frac_cp, 0.75)         # 3 / (3 + 1)
})

test_that("a 2-chromosome closed cycle is rejected by default (needs >= 3 chromosomes)", {
  ## 3-junction closed triangle spanning only chr1 and chr5 (two loci on chr1,
  ## joined by a long intrachromosomal edge): a valid cycle, but not chromoplexy.
  jun <- data.frame(
    chr1 = c("chr1",    "chr5",     "chr1"),    pos1 = c(10000000, 20003000, 10003000),
    chr2 = c("chr5",    "chr1",     "chr1"),    pos2 = c(20000000, 50000000, 50003000),
    stringsAsFactors = FALSE)
  inp <- make_chain_inputs(jun, cn = 2)
  expect_equal(nrow(call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr)), 0L)   # default 3
  res2 <- call_chromoplexy(inp$breakpoints_gr, inp$cnv_gr, min_chromosomes = 2L)
  expect_equal(nrow(res2), 1L)
  expect_equal(res2$n_chromosomes, 2L)
})

test_that("empty breakpoints return the typed 0-row schema", {
  empty_bp <- GenomicRanges::GRanges()
  S4Vectors::mcols(empty_bp) <- S4Vectors::DataFrame(
    WGS_ID = character(), event = character(), svclass = character(),
    PURPLE_CN = numeric())
  cnv_gr <- make_chain_inputs(make_chromoplexy_junctions())$cnv_gr
  res <- call_chromoplexy(empty_bp, cnv_gr)
  expect_equal(nrow(res), 0L)
  expect_true(all(c("WGS_ID", "chromoplexy_id", "topology", "n_junctions",
                    "frac_cp", "events") %in% names(res)))
})
