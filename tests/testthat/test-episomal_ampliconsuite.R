# Tests for the AmpliconSuite episomal caller (R/episomal_ampliconsuite.R).
# All fixtures are synthesised into tempfiles so the tests are fully portable
# (no external AmpliconArchitect / AmpliconClassifier data required).

## ---- fixture writers ------------------------------------------------------
# Write an AmpliconArchitect *_graph.txt from segment/junction data.frames.
#   segments : chr, start, end, cn
#   junctions: cA, pA, sA, cB, pB, sB, jcn, pe, homlen (NA -> "None")
aa_graph_file <- function(segments = NULL, junctions = NULL,
                          path = tempfile(fileext = "_graph.txt")) {
  con <- file(path, "w"); on.exit(close(con))
  writeLines(paste("SequenceEdge: StartPosition, EndPosition, PredictedCopyCount,",
                   "AverageCoverage, Size, NumberReadsMapped"), con)
  if (!is.null(segments)) for (i in seq_len(nrow(segments))) {
    s <- segments[i, ]
    writeLines(sprintf("sequence\t%s:%.0f-\t%s:%.0f+\t%s\t100\t%.0f\t100",
                       s$chr, s$start, s$chr, s$end,
                       format(s$cn, scientific = FALSE), s$end - s$start), con)
  }
  writeLines(paste("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount,",
                   "NumberOfReadPairs, HomologySizeIfAvailable, Homology/InsertionSequence"), con)
  if (!is.null(junctions)) for (i in seq_len(nrow(junctions))) {
    j <- junctions[i, ]
    hom <- if (is.na(j$homlen)) "None" else as.character(j$homlen)
    writeLines(sprintf("discordant\t%s:%.0f%s->%s:%.0f%s\t%s\t%.0f\t%s\tNone",
                       j$cA, j$pA, j$sA, j$cB, j$pB, j$sB,
                       format(j$jcn, scientific = FALSE), j$pe, hom), con)
  }
  path
}

# Write an AmpliconClassifier ecDNA interval BED. `intervals` gives INCLUSIVE
# coordinates; AC's exclusive-end convention (end = last base + 1) is applied
# here, so the caller reads back the same inclusive extent.
ac_bed_file <- function(intervals, path = tempfile(fileext = "_ecDNA_1_intervals.bed")) {
  utils::write.table(
    data.frame(intervals$chr, as.integer(intervals$start), as.integer(intervals$end) + 1L),
    path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  path
}

jx <- function(cA, pA, sA, cB, pB, sB, jcn, pe, homlen = NA_real_)
  data.frame(cA = cA, pA = pA, sA = sA, cB = cB, pB = pB, sB = sB,
             jcn = jcn, pe = pe, homlen = homlen, stringsAsFactors = FALSE)

## ---- read_ampliconsuite ---------------------------------------------------
test_that("read_ampliconsuite parses segments and junctions", {
  g <- aa_graph_file(
    segments  = data.frame(chr = "chr7", start = 55e6, end = 55.5e6, cn = 60),
    junctions = jx("chr7", 55.5e6, "+", "chr7", 55e6, "-", jcn = 50, pe = 200, homlen = 3))
  gr <- read_ampliconsuite(g)
  expect_s3_class(gr$segments, "data.table")
  expect_equal(nrow(gr$segments), 1L)
  expect_equal(gr$segments$cn, 60)
  expect_equal(nrow(gr$junctions), 1L)
  expect_equal(gr$junctions$pe_support, 200)   # col 4 = NumberOfReadPairs
  expect_equal(gr$junctions$jcn, 50)           # col 3 = PredictedCopyCount
  expect_equal(gr$junctions$homlen, 3)         # col 5 numeric
})

test_that("svclass is coordinate-invariant regardless of edge direction", {
  # each junction written 'hi -> lo' (arbitrary AA direction)
  j <- rbind(
    jx("chr7", 200, "+", "chr7", 100, "-", 1, 10),  # lo-,hi+ => DUP
    jx("chr7", 100, "+", "chr7", 200, "-", 1, 10),  # lo+,hi- => DEL
    jx("chr7", 100, "+", "chr7", 200, "+", 1, 10),  # same    => INV
    jx("chr7", 100, "+", "chr8", 200, "-", 1, 10))  # inter   => TRA
  gr <- read_ampliconsuite(aa_graph_file(junctions = j))
  expect_equal(gr$junctions$svclass, c("DUP", "DEL", "INV", "TRA"))
  expect_equal(gr$junctions$lo[1], 100)          # breakends sorted lo/hi
  expect_equal(gr$junctions$hi[1], 200)
})

test_that("read_ampliconsuite handles 'None' homology and missing files", {
  g <- aa_graph_file(junctions = jx("chr7", 200, "+", "chr7", 100, "-", 1, 10, homlen = NA))
  expect_true(is.na(read_ampliconsuite(g)$junctions$homlen))
  empty <- read_ampliconsuite(tempfile())        # non-existent path
  expect_equal(nrow(empty$segments), 0L)
  expect_equal(nrow(empty$junctions), 0L)
})

test_that(".read_ac_intervals converts exclusive end to inclusive", {
  bed <- ac_bed_file(data.frame(chr = "chr7", start = 55e6, end = 55.5e6))
  iv <- .read_ac_intervals(bed)
  expect_equal(iv$start, 55e6)
  expect_equal(iv$end, 55.5e6)                    # end-1 of the written 55.5e6+1
})

## ---- call_episomal_from_ampliconsuite: core calls -------------------------
seg1 <- function(cn = 60) data.frame(chr = "chr7", start = 55e6, end = 55.5e6, cn = cn)
bed1 <- function() ac_bed_file(data.frame(chr = "chr7", start = 55e6, end = 55.5e6))

test_that("a clean single-chr episome with a bracketing boundary DUP is episomal", {
  g <- aa_graph_file(seg1(), jx("chr7", 55e6, "-", "chr7", 55.5e6, "+", jcn = 50, pe = 200))
  r <- call_episomal_from_ampliconsuite(g, bed1())
  expect_true(r$episomal)
  expect_equal(r$status, "episomal")
  expect_false(r$chimeric)
  expect_equal(r$n_chr, 1L)
  expect_equal(r$bd_pe_support, 200)
  expect_equal(r$bd_occupancy, 1, tolerance = 1e-6)
})

test_that("no bracketing DUP -> complex", {
  # a DUP that sits INSIDE the extent (does not reach either edge)
  g <- aa_graph_file(seg1(), jx("chr7", 55.1e6, "-", "chr7", 55.4e6, "+", jcn = 50, pe = 200))
  r <- call_episomal_from_ampliconsuite(g, bed1())
  expect_false(r$episomal)
  expect_equal(r$status, "complex")
  expect_true(is.na(r$bd_chr))
})

test_that("boundary DUP is the highest-pe SPANNING DUP, not beaten by a higher-pe DEL", {
  g <- aa_graph_file(seg1(), rbind(
    jx("chr7", 55e6,  "-", "chr7", 55.5e6, "+", jcn = 50, pe = 200),   # boundary DUP
    jx("chr7", 55.1e6, "+", "chr7", 55.2e6, "-", jcn = 90, pe = 300))) # internal DEL, higher pe
  r <- call_episomal_from_ampliconsuite(g, bed1())
  expect_true(r$episomal)                # DEL does not gate the call
  expect_equal(r$bd_pe_support, 200)     # selector picks the DUP
})

## ---- chimeric -------------------------------------------------------------
test_that("multi-chr with a boundary DUP on each chr + a TRA is chimeric episomal", {
  segs <- rbind(data.frame(chr = "chr7", start = 55e6,  end = 55.5e6,  cn = 60),
                data.frame(chr = "chr8", start = 128e6, end = 128.5e6, cn = 60))
  j <- rbind(
    jx("chr7", 55e6,  "-", "chr7", 55.5e6,  "+", jcn = 50, pe = 200),   # chr7 boundary DUP
    jx("chr8", 128e6, "-", "chr8", 128.5e6, "+", jcn = 50, pe = 180),   # chr8 boundary DUP
    jx("chr7", 55.25e6, "+", "chr8", 128.25e6, "-", jcn = 40, pe = 150))# linking TRA
  bed <- ac_bed_file(rbind(data.frame(chr = "chr7", start = 55e6,  end = 55.5e6),
                           data.frame(chr = "chr8", start = 128e6, end = 128.5e6)))
  r <- call_episomal_from_ampliconsuite(aa_graph_file(segs, j), bed)
  expect_true(r$episomal)
  expect_true(r$chimeric)
  expect_equal(r$n_chr, 2L)
})

test_that("multi-chr WITHOUT a TRA is not episomal", {
  segs <- rbind(data.frame(chr = "chr7", start = 55e6,  end = 55.5e6,  cn = 60),
                data.frame(chr = "chr8", start = 128e6, end = 128.5e6, cn = 60))
  j <- rbind(jx("chr7", 55e6,  "-", "chr7", 55.5e6,  "+", jcn = 50, pe = 200),
             jx("chr8", 128e6, "-", "chr8", 128.5e6, "+", jcn = 50, pe = 180))
  bed <- ac_bed_file(rbind(data.frame(chr = "chr7", start = 55e6,  end = 55.5e6),
                           data.frame(chr = "chr8", start = 128e6, end = 128.5e6)))
  r <- call_episomal_from_ampliconsuite(aa_graph_file(segs, j), bed)
  expect_false(r$episomal)
  expect_equal(r$status, "complex")
})

## ---- bridged / ambiguous downgrade ---------------------------------------
test_that("same-chr highest-pe junction reaching an amplified far region downgrades to bridged", {
  segs <- rbind(data.frame(chr = "chr7", start = 55e6,  end = 55.5e6,  cn = 60),
                data.frame(chr = "chr7", start = 57.9e6, end = 58.1e6, cn = 60)) # amplified far block
  j <- rbind(
    jx("chr7", 55e6,   "-", "chr7", 55.5e6, "+", jcn = 50, pe = 100),  # boundary DUP
    jx("chr7", 55.25e6,"-", "chr7", 58e6,   "+", jcn = 40, pe = 500))  # highest-pe, far end out+amplified
  r <- call_episomal_from_ampliconsuite(aa_graph_file(segs, j), bed1())
  expect_false(r$episomal)
  expect_equal(r$status, "bridged")
  expect_equal(r$bridge_pe_support, 500)
  expect_equal(r$bridge_partner, "chr7:58000000")
})

test_that("a cross-chr TRA far end is the chimeric signature and is NOT downgraded", {
  g <- aa_graph_file(seg1(), rbind(
    jx("chr7", 55e6,   "-", "chr7", 55.5e6, "+", jcn = 50, pe = 100),  # boundary DUP
    jx("chr7", 55.25e6,"+", "chr9", 10e6,   "-", jcn = 40, pe = 500))) # cross-chr, highest pe
  r <- call_episomal_from_ampliconsuite(g, bed1())
  expect_true(r$episomal)
  expect_equal(r$status, "episomal")
})

## ---- occupancy guard ------------------------------------------------------
test_that("a DUP spanning a mostly non-amplified extent is rejected by the occupancy guard", {
  # two tiny amplified intervals ~120 Mb apart; a lone DUP brackets the outer span
  segs <- rbind(data.frame(chr = "chr12", start = 4e6,     end = 4.05e6,   cn = 60),
                data.frame(chr = "chr12", start = 124e6,   end = 124.05e6, cn = 60))
  j <- jx("chr12", 4e6, "-", "chr12", 124.05e6, "+", jcn = 50, pe = 50)
  bed <- ac_bed_file(rbind(data.frame(chr = "chr12", start = 4e6,   end = 4.05e6),
                           data.frame(chr = "chr12", start = 124e6, end = 124.05e6)))
  g <- aa_graph_file(segs, j)

  r <- call_episomal_from_ampliconsuite(g, bed)                    # default 5% guard
  expect_false(r$episomal)
  expect_equal(r$status, "complex")
  expect_lt(r$bd_occupancy, 0.05)

  r0 <- call_episomal_from_ampliconsuite(g, bed, min_occupancy = 0) # guard off
  expect_true(r0$episomal)                                          # the guard is what flips it
})

## ---- aa_cycle_intervals ---------------------------------------------------
test_that("aa_cycle_intervals reads the dominant cycle and drops the source vertex", {
  cyc <- tempfile(fileext = "_cycles.txt")
  writeLines(c(
    "Segment\t1\tchr7\t55000000\t55200000",
    "Segment\t2\tchr7\t55200001\t55500000",
    "Cycle=1;Copy_count=100;Segments=0+,1+,2+",
    "Cycle=2;Copy_count=10;Segments=0+,1+"), cyc)
  iv <- aa_cycle_intervals(cyc)                       # dominant = cycle 1
  expect_equal(nrow(iv), 2L)                          # segment 0 dropped
  expect_equal(sort(iv$chr), c("chr7", "chr7"))
  expect_equal(min(iv$start), 55000000)
  expect_equal(max(iv$end), 55500001)                # end + 1 (exclusive convention)
})

test_that("a *_cycles.txt path can be passed directly as the intervals argument", {
  cyc <- tempfile(fileext = "_cycles.txt")
  writeLines(c(
    "Segment\t1\tchr7\t55000000\t55500000",
    "Cycle=1;Copy_count=100;Segments=0+,1+"), cyc)
  g <- aa_graph_file(seg1(), jx("chr7", 55e6, "-", "chr7", 55.5e6, "+", jcn = 50, pe = 200))
  r <- call_episomal_from_ampliconsuite(g, cyc)       # caller detects _cycles.txt
  expect_true(r$episomal)
  expect_equal(r$n_chr, 1L)
})

## ---- bundled gold standards ----------------------------------------------
test_that("bundled GBM39 (EGFR) is episomal and COLO320DM (MYC) is complex", {
  aa <- function(f) system.file("extdata", "ampliconsuite", f, package = "EpiTracer")
  gbm39 <- call_episomal_from_ampliconsuite(
    aa("GBM39_EGFR_graph.txt"), aa("GBM39_EGFR_ecDNA_1_intervals.bed"))
  colo <- call_episomal_from_ampliconsuite(
    aa("COLO320DM_MYC_graph.txt"), aa("COLO320DM_MYC_ecDNA_1_intervals.bed"))

  expect_true(gbm39$episomal)                 # EGFR simple-excision circle
  expect_equal(gbm39$status, "episomal")
  expect_equal(gbm39$bd_chr, "chr7")
  expect_gt(gbm39$bd_occupancy, 0.9)

  expect_false(colo$episomal)                 # MYC ecDNA is complex
  expect_equal(colo$status, "complex")
  expect_true(is.na(colo$bd_chr))
})

## ---- degenerate inputs ----------------------------------------------------
test_that("missing / empty intervals return an NA (uncallable) row", {
  g <- aa_graph_file(seg1(), jx("chr7", 55e6, "-", "chr7", 55.5e6, "+", 50, 200))
  r <- call_episomal_from_ampliconsuite(g, NA_character_)
  expect_true(is.na(r$episomal))
  expect_equal(nrow(r), 1L)
})

## ---- interval trimming to the amplified extent -----------------------------
test_that("trim_to_amplified rescues an episome whose AC interval is padded with diploid flanks", {
  # amplified EGFR-like body (CN 60) flanked by ~400 kb of diploid (CN 2) either side;
  # the boundary DUP brackets the body but NOT the padded AC interval.
  segs <- rbind(
    data.frame(chr = "chr7", start = 54.6e6,       end = 55.0e6 - 1, cn = 2),   # left diploid flank
    data.frame(chr = "chr7", start = 55.0e6,       end = 55.5e6,     cn = 60),  # amplified body
    data.frame(chr = "chr7", start = 55.5e6 + 1,   end = 56.15e6,    cn = 2))   # right diploid flank
  g   <- aa_graph_file(segs, jx("chr7", 55.0e6, "-", "chr7", 55.5e6, "+", jcn = 55, pe = 200))
  bed <- ac_bed_file(data.frame(chr = "chr7", start = 54.6e6, end = 56.15e6))   # padded ~400 kb each side

  r_trim <- call_episomal_from_ampliconsuite(g, bed)                       # default trim on
  expect_true(r_trim$episomal)
  expect_equal(r_trim$status, "episomal")
  expect_gt(r_trim$bd_occupancy, 0.9)

  r_notrim <- call_episomal_from_ampliconsuite(g, bed, trim_to_amplified = FALSE)
  expect_false(r_notrim$episomal)                                          # DUP fails to bracket the padded interval
  expect_equal(r_notrim$status, "complex")
})

test_that("trimming keeps internal segments and does not fire on a tiny high-CN spike", {
  # a genuine CN~7 body (>=2 kb) with a 300 bp CN~250 spike must not set the peak
  # and trim the body away (regression guard for the sub-kb-spike case).
  segs <- rbind(
    data.frame(chr = "chr8", start = 68.0e6, end = 69.4e6,          cn = 7),
    data.frame(chr = "chr8", start = 82.7540e6, end = 82.7543e6,    cn = 250),  # 300 bp spike
    data.frame(chr = "chr8", start = 82.7543e6 + 1, end = 85.0e6,   cn = 7))
  g   <- aa_graph_file(segs, jx("chr8", 68.0e6, "-", "chr8", 85.0e6, "+", jcn = 6, pe = 40))
  bed <- ac_bed_file(data.frame(chr = "chr8", start = 68.0e6, end = 85.0e6))
  # the whole CN-7 region is retained (spike ignored), so the extent is unchanged
  r <- call_episomal_from_ampliconsuite(g, bed)
  expect_equal(round(r$bd_lo), 68.0e6)   # boundary DUP still spans the full region
  expect_equal(round(r$bd_hi), 85.0e6)
})
