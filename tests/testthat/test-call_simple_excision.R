test_that("input validation rejects non-GRanges and missing metadata", {
  gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(1, 2))
  expect_error(call_simple_excision("nope", gr, gr, gr))
  # ecdna_gr without ID / WGS_ID
  expect_error(call_simple_excision(gr, gr, gr, gr))
})

test_that("a textbook episome is classified as episomal with an excision scar", {
  d <- make_episome_inputs(flank_cn = 2)   # non-gained flanks
  res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
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
  res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
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
  res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
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

test_that("SV-aware seeds link two boxes across a gap only when both breakends are above diploid", {
  # Two amplified boxes 5 Mb apart (> merge_gap), so CN-only detection keeps them
  # separate. A junction with BOTH breakends above diploid links them into one
  # amplicon; a junction with a diploid breakend does not.
  cnv <- GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7", start = c(1, 50000000, 52000001, 57000001),
    end = c(49999999, 52000000, 57000000, 70000000),
    sample = "S1", copyNumber = c(2, 40, 2, 40), ploidy = 2,
    majorAlleleCopyNumber = c(1, 39, 1, 39), minorAlleleCopyNumber = 1),
    keep.extra.columns = TRUE)
  mkbp <- function(cnB) GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7", start = c(51000000, 58000000), end = c(51000000, 58000000),
    WGS_ID = "S1", event = "J1", svclass = "DUP", PURPLE_AF = NA_real_, PURPLE_JCN = 5,
    VF = 100, PURPLE_CN = c(40, cnB), insLen = 0L, HOMLEN = 0L), keep.extra.columns = TRUE)

  s_link <- detect_amplicon_seeds(cnv, breakpoints = mkbp(40))   # both above diploid -> 1 seed
  expect_length(s_link, 1L)
  expect_false(s_link$cn_only[1])                                # has a junction inside
  expect_length(detect_amplicon_seeds(cnv, breakpoints = mkbp(2)), 2L)  # diploid breakend -> 2 seeds

  # cn_only: an amplified box with no breakend anywhere inside it is flagged
  bp_far <- GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7", start = c(51000000, 51500000), end = c(51000000, 51500000),
    WGS_ID = "S1", event = "J1", svclass = "DUP", PURPLE_AF = NA_real_, PURPLE_JCN = 5,
    VF = 100, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L), keep.extra.columns = TRUE)
  s2 <- detect_amplicon_seeds(cnv, breakpoints = bp_far)         # junction only in box1
  expect_true(s2$cn_only[s2$ID == "S1_amp2"])                    # box2 (57-70 Mb) has no junction
})

test_that("detect_amplicon_seeds merges similar-CN segments across a gap but not different-CN ones", {
  # Two amplified segments separated by a ~2 Mb diploid gap (> the 1 Mb `gap`,
  # < the 3 Mb `merge_gap`). Similar CN -> one amplicon broken by an internal
  # deletion (merge); very different CN -> two separate amplicons (keep apart).
  mk <- function(cn2) GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7",
    start = c(1,        50000000, 52500001, 55000001),
    end   = c(49999999, 52500000, 54999999, 60000000),  # amp, gap, amp
    sample = "S1", copyNumber = c(2, 20, 2, cn2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 19, 1, cn2 - 1), minorAlleleCopyNumber = 1),
    keep.extra.columns = TRUE)
  expect_length(detect_amplicon_seeds(mk(22)), 1L)   # CN 20 vs 22 -> merged
  expect_length(detect_amplicon_seeds(mk(60)), 2L)   # CN 20 vs 60 -> separate
})

test_that("SV-aware linking needs an AMPLIFIED junction and tolerates an edge-anchored breakend", {
  # Two amplified boxes (CN 40) 5 Mb apart (> merge_gap), so CN-only keeps them
  # separate. The linking DUP must sit in AMPLIFIED copy number on BOTH breakends
  # to fuse them (> min_cn_ratio * ploidy = 6 here): a merely mildly-gained
  # breakend (above diploid but below the amplified bar) does NOT link. And a
  # boundary-anchoring breakend landing a base OUTSIDE the reduced seed edge still
  # maps to its seed (link_tol), so an edge-defining junction is not lost to an
  # off-by-one (the DO12952T1 EGFR pattern).
  cnv <- GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7", start = c(1, 50000000, 52000001, 57000000),
    end = c(49999999, 52000000, 56999999, 70000000),
    sample = "S1", copyNumber = c(2, 40, 2, 40), ploidy = 2,
    majorAlleleCopyNumber = c(1, 39, 1, 39), minorAlleleCopyNumber = 1),
    keep.extra.columns = TRUE)
  mkbp <- function(cnB, p1 = 51000000) GenomicRanges::makeGRangesFromDataFrame(data.table::data.table(
    seqnames = "chr7", start = c(p1, 58000000), end = c(p1, 58000000),
    WGS_ID = "S1", event = "J1", svclass = "DUP", PURPLE_AF = NA_real_, PURPLE_JCN = 5,
    VF = 100, PURPLE_CN = c(40, cnB), insLen = 0L, HOMLEN = 0L), keep.extra.columns = TRUE)

  expect_length(detect_amplicon_seeds(cnv, breakpoints = mkbp(40)), 1L)  # both amplified -> linked
  expect_length(detect_amplicon_seeds(cnv, breakpoints = mkbp(5)),  2L)  # CN 5: gained, NOT amplified -> not linked
  expect_length(detect_amplicon_seeds(cnv, breakpoints = mkbp(2)),  2L)  # diploid -> not linked
  # edge-anchored: breakend 1 bp before box1's reduced start (50000000) still links
  expect_length(detect_amplicon_seeds(cnv, breakpoints = mkbp(40, p1 = 49999999)), 1L)
})

test_that("ecdna_gr = NULL auto-detects seeds and still calls the episome", {
  d <- make_episome_inputs(flank_cn = 2)
  res <- call_simple_excision(NULL, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_s3_class(res, "data.table")
  expect_true(any(res$episomal == "TRUE"))

  # errors clearly when there is nothing to seed from
  d0 <- make_episome_inputs(flank_cn = 2)
  d0$cnv_gr$copyNumber <- 2
  expect_error(
    call_simple_excision(NULL, d0$breakpoints_gr, d0$cnv_gr, d0$cancer_genes_gr),
    "No focal amplicons"
  )
})

test_that("focal episome on a polysomic chromosome: flank_baseline recovers it", {
  # chr7 whole-chromosome baseline of 4 (polysomy, as in glioblastoma), with a
  # focal EGFR amplicon (CN 50) whose flanks sit at the chr7 baseline (4), and a
  # boundary DUP + excision scar. This IS an episome excised from the polysomic
  # chromosome, but the flanks (4) exceed 1.4 * tumour ploidy (2).
  d <- make_episome_inputs(flank_cn = 4, chr_baseline = 4)

  # default ("chromosome"): flanks compared to the local chr7 baseline (4) -> episomal
  loc <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, mc.cores = 1)
  expect_true(all(loc$episomal == "TRUE"))
  expect_true(all(loc[loc$event == "DUP1", ]$duplication_at_boundary == "TRUE"))

  # "ploidy": flanks compared to global ploidy (2) -> the polysomic flanks read
  # as gained and the episome is missed (the pre-fix behaviour)
  glo <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, mc.cores = 1,
                             flank_baseline = "ploidy")
  expect_true(all(glo$episomal == "FALSE"))
})

test_that("a diploid-flank episome is episomal under BOTH flank_baseline modes", {
  d <- make_episome_inputs(flank_cn = 2, chr_baseline = 2)
  for (fb in c("chromosome", "ploidy")) {
    res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                               d$cancer_genes_gr, mc.cores = 1, flank_baseline = fb)
    expect_true(all(res$episomal == "TRUE"), info = fb)
  }
})

test_that("a flank abutting the boundary DUP breakend is not dropped (<=/>= fix)", {
  # DKFZ-GKS1T1 EGFR pattern: the boundary DUP breakends land exactly on the
  # flank segment edges (prox == left-flank end, dist == right-flank start).
  # A strict < / > in the flank test dropped the diploid flanks -> false negative.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7", start = c(54999999, 55500001), end = c(54999999, 55500001),
    WGS_ID = "S1", event = "DUP1", svclass = "DUP",
    PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")

  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$duplication_at_boundary == "TRUE"))
  expect_true(all(res$episomal == "TRUE"))
})

test_that("a sub-kb segmentation shoulder at the boundary does not mask an episome", {
  # NYGC29T1 / HMF001333T2 pattern: a narrow (sub-kb) intermediate CN segment
  # abuts the boundary DUP breakend at a gained level, while the substantive
  # flank just beyond it is diploid. Reading the sliver as THE flank drops the
  # episome; min_flank_width must skip it.
  d <- make_episome_inputs(flank_cn = 2, chr_baseline = 2,
                           prox_shoulder_cn = 8, prox_shoulder_width = 999L)

  # default min_flank_width (2000) skips the 999 bp shoulder -> episomal
  res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, mc.cores = 1)
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary_has_highest_VF == "TRUE"))
  expect_true(any(res$episomal == "TRUE"))

  # disabling the width filter reads the gained sliver as the flank -> missed
  res0 <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                              d$cancer_genes_gr, mc.cores = 1, min_flank_width = 0)
  expect_true(all(res0$episomal == "FALSE"))
})

test_that("a boundary DUP bridging a large low-copy gap is flagged but STILL episomal", {
  # HMF000387T2 / DO12484T1 pattern: two amplicon-level segments fused by a
  # junction spanning a > seed_gap non-amplified gap. The bridging flag is set for
  # review, but a boundary DUP with non-gained flanks is still called episomal
  # (the flag no longer disqualifies).
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(54000000, 56000000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7", start = c(54000000, 56000000), end = c(54000000, 56000000),
    WGS_ID = "S1", event = "DUP1", svclass = "DUP",
    PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7",
    start = c(1,        54000000, 54400001, 55600000, 56000001),
    end   = c(53999999, 54400000, 55599999, 55999999, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2, 50, 2), ploidy = 2,  # 1.2 Mb diploid gap (> seed_gap)
    majorAlleleCopyNumber = c(1, 49, 1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$duplication_at_boundary == "TRUE"))
  expect_true(all(res$flag_bridging_amplicon == "TRUE"))        # bridging surfaced for review
  expect_true(all(res$episomal == "TRUE"))                      # but not disqualified
})

test_that("a boundary defined by a dominant inversion is excluded (INV out-JCNs the boundary DUP)", {
  # DKFZ-4PGFT2 distal-edge pattern: a fold-back inversion whose breakend sits at
  # the amplicon edge and whose JCN exceeds the boundary DUP's -> not a
  # circularisation -> excluded. A subordinate (lower-JCN) boundary inversion on
  # an otherwise clean boundary-DUP episome does NOT flag.
  d <- make_episome_inputs(flank_cn = 2)   # boundary DUP JCN 40; distal edge 55.5 Mb
  add_inv <- function(jcn) {
    gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(c(55200000, 55500000), width = 1L),
      WGS_ID = "S1", event = "INVB", svclass = "t2tINV", PURPLE_AF = 0.9,
      PURPLE_JCN = jcn, VF = 500, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
    x <- d; x$breakpoints_gr <- c(d$breakpoints_gr, gr); x
  }
  hi <- add_inv(60)                        # INV JCN 60 > boundary DUP JCN 40
  r_hi <- call_simple_excision(hi$ecdna_gr, hi$breakpoints_gr, hi$cnv_gr, hi$cancer_genes_gr, mc.cores = 1)
  expect_true(all(r_hi$flag_inv_at_boundary == "TRUE"))
  expect_true(all(r_hi$episomal == "FALSE"))

  lo <- add_inv(20)                        # INV JCN 20 < boundary DUP JCN 40 -> subordinate
  r_lo <- call_simple_excision(lo$ecdna_gr, lo$breakpoints_gr, lo$cnv_gr, lo$cancer_genes_gr, mc.cores = 1)
  expect_true(all(r_lo$flag_inv_at_boundary == "FALSE"))
  expect_true(all(r_lo$episomal == "TRUE"))
})

test_that("an amplicon terminating within a centromere is a chromosomal bridge (excluded)", {
  # NYGC20T1 pattern: a clean boundary-DUP episome whose distal amplified edge
  # falls INSIDE the centromere -> dicentric / chromosomal-bridge event, excluded.
  d <- make_episome_inputs(flank_cn = 2)  # textbook episome, amplicon 55.0-55.5 Mb
  cen_hit  <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55400000, 56000000))  # covers the distal edge
  cen_miss <- GenomicRanges::GRanges("chr7", IRanges::IRanges(58100000, 62100000))  # true chr7 centromere, far away

  r_hit <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cancer_genes_gr,
                               mc.cores = 1, centromeres = cen_hit)
  expect_true(all(r_hit$flag_chromosomal_bridge == "TRUE"))
  expect_true(all(r_hit$episomal == "FALSE"))

  # centromere elsewhere -> not a bridge -> still episomal
  r_miss <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cancer_genes_gr,
                                mc.cores = 1, centromeres = cen_miss)
  expect_true(all(r_miss$flag_chromosomal_bridge == "FALSE"))
  expect_true(all(r_miss$episomal == "TRUE"))
})

test_that("two different boundary DUPs (not the same SV) are episomal type 2, not episomal", {
  # C3N-01334T1 pattern: the proximal edge is reached by one DUP and the distal
  # edge by a DIFFERENT, interleaving DUP. The boundary is not a single connected
  # junction, so it is not straightforward episomal -- flagged "episomal type 2".
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7", start = c(55000000, 55300000, 55200000, 55500000),
    end = c(55000000, 55300000, 55200000, 55500000),
    WGS_ID = "S1", event = c("DUPA", "DUPA", "DUPB", "DUPB"), svclass = "DUP",
    PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$episomal_type2 == "TRUE"))
  expect_true(all(res$episomal == "FALSE"))

  # control: a SINGLE DUP spanning both edges is connected -> straightforward episomal
  res1 <- call_simple_excision(make_episome_inputs(flank_cn = 2)$ecdna_gr,
                              make_episome_inputs(flank_cn = 2)$breakpoints_gr,
                              make_episome_inputs(flank_cn = 2)$cnv_gr,
                              make_episome_inputs(flank_cn = 2)$cancer_genes_gr, mc.cores = 1)
  expect_true(all(res1$episomal == "TRUE"))
  expect_true(all(res1$episomal_type2 == "FALSE"))
})

test_that("a passenger DUP reaching into diploid sequence is not a boundary candidate (still episomal)", {
  # DO12952T1 EGFR pattern: one clean boundary DUP spans both amplicon edges, but
  # a low-JCN passenger tandem duplication pokes a breakend just PAST the distal
  # edge (amplified there, within 10 kb) while its FAR breakend sits out in diploid
  # sequence. Because that passenger is not amplified on BOTH ends it is not an
  # intra-amplicon junction, so it must not be chosen as the distal border (the
  # outermost-breakend rule would otherwise pick it and mis-call the connected
  # episome as type 2). The episome stays type-1 episomal.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7",
    start = c(55000000, 55500000,   # boundary DUP: both edges, both amplified
              55500100, 56000000),  # passenger DUP: near edge (amplified) + far diploid
    end   = c(55000000, 55500000, 55500100, 56000000),
    WGS_ID = "S1",
    event   = c("DUP1", "DUP1", "PDUP", "PDUP"),
    svclass = "DUP",
    PURPLE_AF = 0.9,
    PURPLE_JCN = c(40, 40, 1, 1),
    VF = c(1000, 1000, 20, 20),
    PURPLE_CN = c(50, 50, 50, 2),   # passenger far breakend is diploid (CN 2)
    insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000000, 55500001), end = c(54999999, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$episomal == "TRUE"))
  expect_true(all(res$episomal_type2 == "FALSE"))
  # the boundary DUP is flagged as the boundary, the passenger is not
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary == "TRUE"))
  expect_true(all(res[res$event == "PDUP", ]$duplication_at_boundary == "FALSE"))
})

test_that("a non-amplified flank does not disqualify an episome (lost chromosome / modest gain)", {
  # In the per-chromosome flank test a flank counts as "gained" only if it is
  # itself amplified (>= min_cn_ratio * ploidy) or above the chromosome-gain
  # level. Two failure modes this guards against:
  #   (i) LOST chromosome -- the background is below diploid, so a near-diploid
  #       flank would otherwise read as gained vs a depressed baseline
  #       (DO12742T1 PDGFRA, chr4 ~CN 1.3);
  #   (ii) a MODEST flank a little above diploid abutting a very-high-CN circle
  #       (C3L-02504T1 PDGFRA, distal flank ~CN 3.4 beside a CN-30 amplicon).
  # Control: a flank that IS amplified still disqualifies.
  lost <- make_episome_inputs(flank_cn = 2, chr_baseline = 1)   # diploid flanks, lost-chr background
  r_lost <- call_simple_excision(lost$ecdna_gr, lost$breakpoints_gr, lost$cnv_gr,
                                lost$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(r_lost$episomal == "TRUE"))

  modest <- make_episome_inputs(flank_cn = 4, chr_baseline = 2) # flank CN 4 (< 3*ploidy = 6)
  r_modest <- call_simple_excision(modest$ecdna_gr, modest$breakpoints_gr, modest$cnv_gr,
                                  modest$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(r_modest$episomal == "TRUE"))

  amp <- make_episome_inputs(flank_cn = 8, chr_baseline = 2)    # flank CN 8 (>= 6) -> still gained
  r_amp <- call_simple_excision(amp$ecdna_gr, amp$breakpoints_gr, amp$cnv_gr,
                               amp$cancer_genes_gr, ext = 1e7, mc.cores = 1)
  expect_true(all(r_amp$episomal == "FALSE"))
})

test_that("a boundary DUP is anchored at the edge, not to a neighbouring amplified block", {
  # DO12370T1 chr6 pattern: a focal episome (55.0-55.5 Mb) with a single boundary
  # DUP, and a SEPARATE amplified block 4.5 Mb away (60.0-60.2 Mb) that carries
  # its own amplified DUP. Because breakpoints are gathered from a +/- ext window,
  # the neighbouring block's DUP breakend is in scope; an unbounded "past the
  # distal edge" border test would pick it and mis-call the clean single-boundary
  # episome as type 2. Anchoring the border to within 10 kb of the amplicon edge
  # keeps only the true boundary DUP -> episomal.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000000, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7",
    start = c(55000000, 55500000,   # focal boundary DUP (both edges)
              60000000, 60200000),  # neighbouring block's amplified DUP, 4.5 Mb away
    end   = c(55000000, 55500000, 60000000, 60200000),
    WGS_ID = "S1",
    event   = c("DUP1", "DUP1", "NDUP", "NDUP"),
    svclass = "DUP",
    PURPLE_AF = 0.9, PURPLE_JCN = c(40, 40, 20, 20), VF = c(1000, 1000, 500, 500),
    PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)                    # all four breakends amplified
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cnv <- data.table::data.table(
    seqnames = "chr7",
    start = c(1,        55000000, 55500001, 60000000, 60200001),
    end   = c(54999999, 55500000, 59999999, 60200000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, ext = 1e7, mc.cores = 1)
  expect_true(all(res$episomal == "TRUE"))
  expect_true(all(res$episomal_type2 == "FALSE"))
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary == "TRUE"))
})

test_that("interleaved high-JCN founder DUP+INV across a gap is flagged a micronucleus (annotation only, still episomal)", {
  # BCCA29T1 pattern: two amplified boxes fused by a founder DUP that spans the
  # diploid gap between them, with a high-JCN inversion starting INSIDE that DUP.
  # Interleaved founder DUP+INV across a gap = intra-micronucleus chromothripsis;
  # this is RECORDED (flag_micronucleus) but is annotation only and does NOT
  # disqualify -- an episome can be internally rearranged. The call stays episomal.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(54000000, 56000000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  cnv <- data.table::data.table(
    seqnames = "chr7",
    start = c(1,        54000000, 54400001, 55600000, 56000001),
    end   = c(53999999, 54400000, 55599999, 55999999, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2, 50, 2), ploidy = 2,  # 1.2 Mb gap between the boxes
    majorAlleleCopyNumber = c(1, 49, 1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  mkbp <- function(with_inv) {
    bp <- data.table::data.table(
      seqnames = "chr7", start = c(54000000, 56000000), end = c(54000000, 56000000),
      WGS_ID = "S1", event = "DUP1", svclass = "DUP",
      PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
    if (with_inv)
      bp <- rbind(bp, data.table::data.table(
        seqnames = "chr7", start = c(54800000, 55700000), end = c(54800000, 55700000),
        WGS_ID = "S1", event = "INV1", svclass = "h2hINV",
        PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 900, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L))
    GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  }
  res <- call_simple_excision(ecdna_gr, mkbp(TRUE), cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$flag_micronucleus == "TRUE"))
  expect_true(all(res$episomal == "TRUE"))   # annotation only; no longer disqualifying

  # same gap-spanning boundary DUP WITHOUT the interleaving inversion -> episomal
  res0 <- call_simple_excision(ecdna_gr, mkbp(FALSE), cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res0$flag_micronucleus == "FALSE"))
  expect_true(all(res0$episomal == "TRUE"))
})

test_that("an internal inversion out-copying the boundary DUP is flagged (annotation only, still episomal)", {
  # An internal inversion carrying a higher VF than the boundary DUP is the
  # inverted-duplication signature (GMPS/HMF001167T2, HMF001611T2, DO11794T2
  # pattern). This is RECORDED (flag_internal_inversion) but is annotation only
  # and does NOT disqualify -- an episome can be internally rearranged -- so the
  # amplicon stays episomal.
  d_inv <- add_internal_sv(make_episome_inputs(flank_cn = 2), cls = "h2hINV", vf = 1200)
  r_inv <- call_simple_excision(d_inv$ecdna_gr, d_inv$breakpoints_gr, d_inv$cnv_gr, d_inv$cancer_genes_gr, mc.cores = 1)
  expect_true(all(r_inv$flag_internal_inversion == "TRUE"))
  expect_true(all(r_inv$episomal == "TRUE"))   # annotation only; no longer disqualifying
})

test_that("an internal deletion out-VFing the boundary DUP is ancestral (still episomal)", {
  # A deletion REMOVES sequence: a high-VF internal deletion is an ancestral /
  # smaller-circle event, not a rival formation mechanism, so the amplicon stays
  # episomal; the flag is informational.
  d_del <- add_internal_sv(make_episome_inputs(flank_cn = 2), cls = "DEL", vf = 1200)
  r_del <- call_simple_excision(d_del$ecdna_gr, d_del$breakpoints_gr, d_del$cnv_gr, d_del$cancer_genes_gr, mc.cores = 1)
  expect_true(all(r_del$flag_internal_sv_high_vf == "TRUE"))
  expect_true(all(r_del$flag_internal_inversion == "FALSE"))
  expect_true(all(r_del$episomal == "TRUE"))
})

test_that("a flank ending 1-2 bp past the boundary breakend is not dropped (amplicon-edge anchoring)", {
  # C3L-03405T1 / C3N-02785T1 pattern: the boundary DUP breakend lands a base or
  # two INSIDE the flanking CN segment (the segment ends at breakend+1), so keying
  # the flank test off the breakend coordinate missed the diploid flank on an
  # off-by-one. Anchoring the flank test to the CN-amplicon edge fixes it.
  ecdna_gr <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55000002, 55500000),
                                     ID = "S1_amp1", WGS_ID = "S1")
  bp <- data.table::data.table(
    seqnames = "chr7", start = c(55000000, 55500001), end = c(55000000, 55500001),
    WGS_ID = "S1", event = "DUP1", svclass = "DUP",
    PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  # proximal flank ends at breakend+1 (55000001); distal flank starts at breakend+1 (55500002)
  cnv <- data.table::data.table(
    seqnames = "chr7", start = c(1, 55000002, 55500002), end = c(55000001, 55500000, 70000000),
    sample = "S1", copyNumber = c(2, 50, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  res <- call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cgg, mc.cores = 1)
  expect_true(all(res$duplication_at_boundary == "TRUE"))
  expect_true(all(res$episomal == "TRUE"))
})

test_that("boundary-DUP breakends in reversed row order are still called episomal", {
  # C3L-03405T1 / C3N-02785T1 pattern: breakpoints_gr stacks all chrom1 then all
  # chrom2 breakends, so the first-listed boundary breakend is not necessarily
  # the lower coordinate. Assigning prox/dist by row order put the "before" flank
  # inside the amplicon (gained) and dropped the episome; prox/dist must be taken
  # by genomic position.
  d <- make_episome_inputs(reverse_boundary = TRUE)
  res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                             d$cancer_genes_gr, mc.cores = 1)
  expect_true(all(res[res$event == "DUP1", ]$duplication_at_boundary == "TRUE"))
  expect_true(all(res$episomal == "TRUE"))
})

test_that("circularisation-junction microhomology classifies NHEJ / MMEJ / HR", {
  # The boundary DUP is the self-ligation junction; its breakpoint homology
  # length classifies the inferred DSB-repair pathway (Eugen-Olsen et al.,
  # Nucleic Acids Res 2025, doi:10.1093/nar/gkaf122): < 2 bp = NHEJ (blunt),
  # 2..<14 bp = MMEJ (microhomology / alt-EJ), >= 14 bp = HR (long homology).
  for (hl_class in list(c(0L, "NHEJ"), c(5L, "MMEJ"), c(20L, "HR"))) {
    d <- make_episome_inputs(boundary_homlen = as.integer(hl_class[[1]]))
    res <- call_simple_excision(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr,
                               d$cancer_genes_gr, mc.cores = 1)
    expect_true(all(res$episomal == "TRUE"))
    bd <- res[res$duplication_at_boundary == "TRUE", ]
    expect_true(all(bd$boundary_homology == as.numeric(hl_class[[1]])))
    expect_true(all(bd$junction_homology_class == hl_class[[2]]))
  }
})

test_that("junction_homology_class is NA when there is no boundary DUP", {
  # No boundary DUP -> not episomal, and the homology annotation stays NA.
  d <- make_episome_inputs(flank_cn = 2)
  no_dup <- d$breakpoints_gr[d$breakpoints_gr$svclass != "DUP"]
  res <- call_simple_excision(d$ecdna_gr, no_dup, d$cnv_gr,
                             d$cancer_genes_gr, mc.cores = 1)
  expect_true(all(res$episomal == "FALSE"))
  expect_true(all(is.na(res$junction_homology_class)))
})
