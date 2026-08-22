## Breakage-fusion-bridge (classical, staircase-gated): an intrachromosomal
## fold-back amplicon, TERMINAL on one arm with a DISTAL DELETION, and showing a
## copy-number STAIRCASE (fold-backs spread across several stepped CN levels) --
## the feature that separates iterative BFB from a single BRF event / focal ecDNA.

CL  <- c(chr7 = 159345973)
CEN <- GenomicRanges::GRanges("chr7", IRanges::IRanges(58100000, 62100000))

## p-arm amplicon [lo,up] below the centromere.
##   distal_cn : CN of the terminal region [0,lo)  (1 = deleted/BFB; 2 = retained)
##   flat      : single CN level across the amplicon (no staircase)
##   clustered : the two fold-backs sit close together (BRF-like, low spread)
##   foldback  : inversions (TRUE) vs a DUP (FALSE)
##   tra_partner_cn : if > 0, add a chr12 TRA partner at this CN (amplified => not intrachromosomal)
bfb_inputs <- function(lo = 10e6, up = 20e6, distal_cn = 1, flat = FALSE,
                       clustered = FALSE, foldback = TRUE, tra_partner_cn = 0) {
  w <- up - lo; b1 <- lo + round(w / 3); b2 <- lo + round(2 * w / 3)
  amp <- if (flat) data.frame(start = lo, end = up, cn = 40)
         else      data.frame(start = c(lo, b1 + 1, b2 + 1), end = c(b1, b2, up), cn = c(40, 25, 12))
  cnv <- data.frame(seqnames = "chr7",
    start = c(1, amp$start, up + 1), end = c(lo - 1, amp$end, 159000000),
    copyNumber = c(distal_cn, amp$cn, 2))
  cnv$sample <- "S1"; cnv$ploidy <- 2
  cnv$majorAlleleCopyNumber <- pmax(0, cnv$copyNumber - 1)
  cnv$minorAlleleCopyNumber <- ifelse(cnv$copyNumber < 1.5, 0, 1)
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

  cls <- if (foldback) "h2hINV" else "DUP"
  fbp <- if (clustered) c(lo + 0.10 * w, lo + 0.14 * w) else c(lo + 0.20 * w, lo + 0.80 * w)
  bp <- data.table::data.table(seqnames = "chr7", start = fbp, end = fbp, WGS_ID = "S1",
    event = c("FB1", "FB2"), svclass = cls, PURPLE_AF = 0.9, PURPLE_JCN = 30,
    VF = 500, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L)
  if (tra_partner_cn > 0)
    bp <- rbind(bp, data.table::data.table(seqnames = c("chr7", "chr12"),
      start = c((lo + up) / 2, 57e6), end = c((lo + up) / 2, 57e6), WGS_ID = "S1",
      event = "TRA1", svclass = "TRA", PURPLE_AF = 0.9, PURPLE_JCN = 30, VF = 500,
      PURPLE_CN = c(40, tra_partner_cn), insLen = 0L, HOMLEN = 0L))
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)
  cgg <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR")
  list(ecdna_gr = GenomicRanges::GRanges("chr7", IRanges::IRanges(lo, up), ID = "S1_amp1", WGS_ID = "S1"),
       breakpoints_gr = breakpoints_gr, cnv_gr = cnv_gr, cgg = cgg)
}
run <- function(d) call_bfb(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg,
                            centromeres = CEN, chrom_lengths = CL, mc.cores = 1)

test_that("a terminal, distally-deleted, staircase fold-back amplicon is BFB", {
  r <- run(bfb_inputs())
  expect_true(all(r$bfb == "TRUE"))
  expect_true(all(r$bfb_anchor == "p-telomere"))
  expect_true(all(r$n_foldbacks >= 2))
})

test_that("clustered fold-backs (BRF-like, no staircase spread) are NOT BFB", {
  r <- run(bfb_inputs(clustered = TRUE))
  expect_true(all(r$bfb == "FALSE"))
})

test_that("a flat plateau (single CN level, no staircase) is NOT BFB", {
  r <- run(bfb_inputs(flat = TRUE))
  expect_true(all(r$bfb == "FALSE"))
})

test_that("no distal deletion (diploid terminus) is NOT BFB", {
  r <- run(bfb_inputs(distal_cn = 2))
  expect_true(all(r$bfb == "FALSE"))
})

test_that("an amplified TRA partner (not intrachromosomal) is NOT BFB", {
  r <- run(bfb_inputs(tra_partner_cn = 50))
  expect_true(all(r$bfb == "FALSE"))
})

test_that("no fold-back (DUP only) is NOT BFB", {
  r <- run(bfb_inputs(foldback = FALSE))
  expect_true(all(r$bfb == "FALSE"))
  expect_true(all(r$n_foldbacks == 0))
})

test_that("without chrom_lengths/centromeres the BFB annotation is disabled", {
  d <- bfb_inputs()
  r <- call_bfb(d$ecdna_gr, d$breakpoints_gr, d$cnv_gr, d$cgg, mc.cores = 1)
  expect_true(all(r$bfb == "FALSE"))
})
