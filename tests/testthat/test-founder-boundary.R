## Founder-boundary nomination (call_founder_boundary): VF-stratified test for a
## circularisation DUP. A clean episome's founder (max-VF) junction is the
## boundary DUP; if an internal fold-back inversion out-copies it, the boundary
## DUP is demoted to cluster 1 (stratum 2) and the founder becomes the inversion
## -- the inverted-duplication signature that disqualifies simple excision.

agg <- function(inputs) {
  r <- call_founder_boundary(inputs$ecdna_gr, inputs$breakpoints_gr, inputs$cnv_gr,
                             inputs$cancer_genes_gr, mc.cores = 1)
  d <- data.table::as.data.table(as.data.frame(r))
  d[!is.na(founder_class),
    .(founder_class = founder_class[1], early_boundary_dup = early_boundary_dup[1],
      early_dup_stratum = early_dup_stratum[1], interlocus_tra = interlocus_tra[1]), by = ID]
}

test_that("a clean episome has a DUP founder that is the early boundary (stratum 1)", {
  a <- agg(make_episome_inputs())
  expect_equal(a$founder_class, "DUP")
  expect_true(all(a$early_boundary_dup == "TRUE"))
  expect_equal(a$early_dup_stratum, 1L)
})

test_that("an internal fold-back that out-copies the boundary makes the founder an inversion", {
  ## boundary DUP VF is 1000 (make_episome_inputs); add an internal h2h inversion at VF 2000
  a <- agg(add_internal_sv(make_episome_inputs(), "h2hINV", vf = 2000))
  expect_true(a$founder_class %in% c("h2hINV", "t2tINV"))   # founder is now the inversion
  ## the boundary DUP still exists but is demoted out of the founder slot
  expect_true(all(a$early_boundary_dup == "TRUE"))
  expect_gt(a$early_dup_stratum, 1L)                        # cluster 1, not founder
})

test_that("a low-VF internal DUP does not displace the boundary DUP founder", {
  a <- agg(add_internal_sv(make_episome_inputs(), "DUP", vf = 100))
  expect_equal(a$founder_class, "DUP")
  expect_equal(a$early_dup_stratum, 1L)
})

test_that("a high-VF deletion is ancestral, not the founder", {
  ## a deletion removes sequence and cannot found an amplification; even at a VF
  ## exceeding the boundary DUP it must not be called the founder
  a <- agg(add_internal_sv(make_episome_inputs(), "DEL", vf = 5000))
  expect_equal(a$founder_class, "DUP")             # boundary DUP founds it; the DEL is ancestral
  expect_true(all(a$early_boundary_dup == "TRUE"))
})

test_that("an inter-locus translocation is excluded from the founder call and flagged", {
  ## a high-VF TRA to another locus must NOT become the founder; the amplicon is
  ## still assessed on its own boundary DUP, and the TRA is reported separately
  inp <- make_episome_inputs()
  tra <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55250000, width = 1L),
    WGS_ID = "S1", event = "TRA1", svclass = "TRA", PURPLE_AF = 0.9, PURPLE_JCN = 30,
    VF = 5000, PURPLE_CN = 50, insLen = 0L, HOMLEN = 0L)
  inp$breakpoints_gr <- c(inp$breakpoints_gr, tra)
  a <- agg(inp)
  expect_equal(a$founder_class, "DUP")                 # TRA excluded; boundary DUP is founder
  expect_true(all(a$early_boundary_dup == "TRUE"))
  expect_true(all(a$interlocus_tra == "TRUE"))         # TRA flagged for the recombination hypothesis
})
