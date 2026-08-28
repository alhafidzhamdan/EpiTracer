## Tests for the ecDNA simulation suite. The contract being asserted is that the
## simulator emits, by construction, the signatures the callers look for -- so a
## caller regression shows up here as a simulated amplicon that stops being
## recovered, and a simulator regression shows up as a lost signature.

q <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("sim_episome builds a circle with the boundary-duplication signature", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T1", copies = 50, seed = 1)

  expect_s3_class(ec, "ecdna_sim")
  expect_true(ec$circular)
  expect_length(ec$fragments, 1L)

  ## the circularisation junction is duplication-orientation (-/+) and spans the
  ## amplicon; the excision scar is deletion-orientation (+/-) just outside it
  amp <- ec$junctions[origin == "amplicon"]
  expect_equal(nrow(amp), 1L)
  expect_equal(amp$svclass, "DUP")
  expect_equal(amp$strand1, "-"); expect_equal(amp$strand2, "+")
  expect_equal(amp$start1, GenomicRanges::start(ec$fragments)[1])
  expect_equal(amp$start2, GenomicRanges::end(ec$fragments)[1])

  scar <- ec$junctions[origin == "host"]
  expect_equal(scar$svclass, "DEL")
  expect_lt(scar$start1, amp$start1)
  expect_gt(scar$start2, amp$start2)
})

test_that("scar and source retention can be switched off", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T1", scar = FALSE, seed = 1)
  expect_equal(nrow(ec$junctions[origin == "host"]), 0L)

  ec2 <- sim_episome(seed_locus("EGFR"), sample = "T1", retain_source = FALSE, seed = 1)
  expect_null(ec2$excised)
})

test_that("noiseless emission reproduces the ground-truth copy number exactly", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T2", copies = 50,
                    host_ploidy = 2, seed = 1)
  inp <- sim_to_epitracer(ec, noise = sim_noise(0), seed = 1)
  cn <- gr2dt(inp$cnv_gr)

  ## flanks sit at the host baseline, the amplicon at copies (one allele excised,
  ## the other retained: 50 + 1)
  expect_equal(sort(round(unique(cn$copyNumber))), c(2, 51))
  amp <- cn[copyNumber > 10]
  expect_equal(amp$minorAlleleCopyNumber, 1)
  expect_equal(amp$majorAlleleCopyNumber, 50)

  ## no false positives, no dropout at zero noise: exactly the two true junctions
  bp <- gr2dt(inp$breakpoints_gr)
  expect_equal(length(unique(bp$event)), 2L)
  expect_setequal(unique(bp$svclass), c("DUP", "DEL"))
})

test_that("a simulated episome is called episomal by call_simple_excision", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T3", copies = 50, seed = 1)
  inp <- sim_to_epitracer(ec, seed = 1)
  res <- q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                                inp$cnv_gr, inp$cancer_genes_gr))

  expect_gt(nrow(res), 0L)
  expect_true(any(res$episomal == "TRUE"))
  expect_true(any(res$duplication_at_boundary == "TRUE"))
  expect_true(any(res$has_excision_scar == "TRUE"))
})

test_that("the repair pathway sets the reported breakpoint homology", {
  for (p in c("nhej", "mmej", "nahr")) {
    ec <- sim_episome(seed_locus("EGFR"), sample = "T4", homology = p, seed = 2)
    bp <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(0), seed = 2)$breakpoints_gr)
    hom <- bp[svclass == "DUP"]$HOMLEN
    switch(p,
      nhej = expect_true(all(hom <= 1)),
      mmej = expect_true(all(hom >= 2 & hom <= 10)),
      nahr = expect_true(all(hom >= 14)))
  }
})

test_that("micronucleation shatters the circle and randomises junction orientation", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T5", copies = 40, seed = 3)
  mn <- sim_micronucleation(ec, n_breaks = 30, seed = 3)

  expect_gt(length(mn$fragments), length(ec$fragments))
  expect_gt(nrow(mn$junctions), nrow(ec$junctions))
  expect_equal(mn$round, 1L)
  expect_true(mn$circular)

  ## all four intrachromosomal orientations appear -- the random-joins hallmark
  expect_true(all(c("DEL", "DUP", "h2hINV", "t2tINV") %in% mn$junctions$svclass))
})

test_that("the oncogene survives every round of shattering", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T6", seed = 4)
  gene <- ec$preserved
  for (r in 1:5) {
    ec <- sim_micronucleation(ec, n_breaks = 25, del_p = 0.6, seed = 100 + r)
    expect_gt(sum(GenomicRanges::countOverlaps(ec$fragments, gene)), 0L)
  }
})

test_that("copies per cell are capped so peak copy number respects max_cn", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T7", copies = 60, seed = 5)
  ev <- sim_evolve(ec, rounds = 5, n_breaks = 20, amplify = 2, max_cn = 150, seed = 5)
  expect_lte(summary(ev)$max_cn, 150 + 1e-6)
})

test_that("sim_evolve returns the whole trajectory with keep_all", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T8", seed = 6)
  traj <- sim_evolve(ec, rounds = 3, n_breaks = 10, keep_all = TRUE, seed = 6)

  expect_length(traj, 4L)
  expect_equal(vapply(traj, function(z) z$round, integer(1)), 0:3)
  ## junction burden accumulates with rounds
  nj <- vapply(traj, function(z) nrow(z$junctions), integer(1))
  expect_true(all(diff(nj) > 0))
})

test_that("an evolved amplicon is called chromothriptic, not episomal", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T9", copies = 40, seed = 7)
  ev <- sim_evolve(ec, rounds = 2, n_breaks = 20, seed = 7)
  inp <- sim_to_epitracer(ev, seed = 7)

  ct <- q(call_chromothripsis(inp$ecdna_gr, inp$breakpoints_gr,
                              inp$cnv_gr, inp$cancer_genes_gr))
  se <- q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                               inp$cnv_gr, inp$cancer_genes_gr))
  expect_true(any(ct$chromothripsis == "TRUE"))
  expect_false(any(se$episomal == "TRUE"))
})

test_that("fusing two episomes gives a chimera with amplified translocations", {
  a <- sim_episome(seed_locus("EGFR"), sample = "T10", copies = 40, seed = 8)
  b <- sim_episome(seed_locus("CDK4"), sample = "T10", copies = 30, seed = 9)
  ch <- sim_fuse_episomes(a, b, n_breaks = 12, seed = 10)

  expect_length(unique(as.character(GenomeInfoDb::seqnames(ch$fragments))), 2L)
  expect_gt(sum(ch$junctions$svclass == "TRA"), 0L)
  ## both parents contribute fragments, and both oncogenes ride the circle
  expect_length(unique(ch$fragments$origin), 2L)
  expect_setequal(ch$preserved$gene, c("EGFR", "CDK4"))

  inp <- sim_to_epitracer(ch, seed = 10)
  mn <- q(call_micronucleation(inp$ecdna_gr, inp$breakpoints_gr,
                               inp$cnv_gr, inp$cancer_genes_gr))
  expect_true(any(mn$micronucleation == "TRUE"))
})

test_that("sim_noise scales the technical model and zero noise is exact", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T11", copies = 50, seed = 11)
  exact <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(0), seed = 1)$cnv_gr)
  noisy <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(3), seed = 1)$cnv_gr)

  expect_equal(round(max(exact$copyNumber)), 51)
  expect_gt(stats::sd(noisy$copyNumber - exact$copyNumber), 0)

  ## low purity and shallow depth cut read support
  deep <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(0, depth = 100, purity = 0.9),
                                 seed = 1)$breakpoints_gr)
  shallow <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(0, depth = 20, purity = 0.3),
                                    seed = 1)$breakpoints_gr)
  expect_gt(max(deep$VF), max(shallow$VF))
})

test_that("sim_cohort labels every sample and emits usable caller inputs", {
  co <- sim_cohort(n = 6, seed = 12)

  expect_equal(nrow(co$truth), 6L)
  expect_true(all(co$truth$class %in% c("episomal", "chromothriptic", "chimeric")))
  expect_length(unique(co$truth$WGS_ID), 6L)
  expect_true(all(c("ecdna_gr", "breakpoints_gr", "cnv_gr", "cancer_genes_gr")
                  %in% names(co)))

  ## the inputs satisfy the callers' column contract
  expect_true(all(c("WGS_ID", "event", "svclass", "PURPLE_AF", "PURPLE_JCN",
                    "VF", "PURPLE_CN", "insLen", "HOMLEN")
                  %in% names(S4Vectors::mcols(co$breakpoints_gr))))
  expect_true(all(c("sample", "copyNumber", "ploidy", "majorAlleleCopyNumber",
                    "minorAlleleCopyNumber")
                  %in% names(S4Vectors::mcols(co$cnv_gr))))
})

test_that("sim_benchmark recovers the trajectory class each caller is built for", {
  co <- sim_cohort(n = 12, seed = 13)
  bm <- sim_benchmark(co)

  expect_equal(nrow(bm$calls), 12L)
  s <- bm$summary
  ## simple excision fires on born-and-unshattered amplicons, and on those only
  expect_equal(s[class == "episomal"]$simple_excision, 1)
  expect_equal(s[class == "chromothriptic"]$simple_excision, 0)
  ## shattered amplicons are called chromothriptic; chimeras micronucleated
  expect_equal(s[class == "chromothriptic"]$chromothripsis, 1)
  expect_equal(s[class == "chimeric"]$micronucleation, 1)
})

test_that("simulated amplicons render as plotter inputs", {
  ec <- sim_evolve(sim_episome(seed_locus("EGFR"), sample = "T12", copies = 40,
                               seed = 14),
                   rounds = 1, n_breaks = 12, seed = 14)
  pin <- sim_to_plot_inputs(ec, seed = 14)

  expect_named(pin, c("cnv_data", "sv_data", "wgd_data"))
  expect_gt(nrow(pin$sv_data), 0L)
  expect_true(all(c("chrom1", "start1", "chrom2", "start2", "strand1", "strand2",
                    "svclass", "VF", "JCN", "sample") %in% names(pin$sv_data)))
  p <- q(plot_sv_linear(sample = "T12", cnv_data = pin$cnv_data,
                        sv_data = pin$sv_data, wgd_data = pin$wgd_data,
                        chromosome = "chr7"))
  expect_s3_class(p, "ggplot")
})

test_that("sim_replicate raises the copy level without touching the structure", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T13", copies = 5, seed = 20)
  rep5 <- sim_replicate(ec, rounds = 5, fold = 2)

  ## structure is byte-identical; only the level moved
  expect_equal(rep5$copies, 5 * 2^5)
  expect_equal(rep5$fragments, ec$fragments)
  expect_equal(rep5$junctions, ec$junctions)
  expect_equal(rep5$round, ec$round)
  expect_equal(summary(rep5)$n_junctions, summary(ec)$n_junctions)
})

test_that("a replicating episome gains boundary read support and no new junctions", {
  ## sim_noise(0) is the exact, deterministic rendering: no copy-number error, no
  ## read-support dispersion, no dropout and no false-positive junctions
  clean <- sim_noise(0)
  ec <- sim_episome(seed_locus("EGFR"), sample = "T14", copies = 8, seed = 21)
  traj <- sim_replicate(ec, rounds = 4, fold = 2, keep_all = TRUE)

  boundary <- t(vapply(traj, function(x) {
    bp <- gr2dt(sim_to_epitracer(x, noise = clean, seed = 21)$breakpoints_gr)
    dup <- bp[svclass == "DUP"][1]
    c(n_event = length(unique(bp$event)), vf = dup$VF, jcn = dup$PURPLE_JCN)
  }, numeric(3)))

  ## exactly two junctions throughout: the boundary DUP and the excision scar
  expect_true(all(boundary[, "n_event"] == 2))
  ## read support and junction copy number both climb monotonically
  expect_true(all(diff(boundary[, "vf"]) > 0))
  expect_true(all(diff(boundary[, "jcn"]) > 0))
  ## and support tracks the copy level exactly: four doublings multiply VF by 16
  expect_equal(unname(boundary[5, "vf"] / boundary[1, "vf"]), 16)
})

test_that("a clean episome carries no internal SVs when false positives are off", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T15", copies = 40, seed = 22)
  bp <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(fp_rate = 0, dropout = 0),
                               seed = 22)$breakpoints_gr)

  expect_setequal(unique(bp$svclass), c("DUP", "DEL"))
  amp <- ec$fragments
  ## no breakend strictly inside the circle -- every one sits at a boundary
  inside <- bp[start > GenomicRanges::start(amp)[1] &
                 start < GenomicRanges::end(amp)[1]]
  expect_equal(nrow(inside), 0L)
})

test_that("replication alone never makes an episome look chromothriptic", {
  clean <- sim_noise(fp_rate = 0, dropout = 0)
  ec <- sim_episome(seed_locus("EGFR"), sample = "T16", copies = 10, seed = 23)
  big <- sim_replicate(ec, rounds = 5, fold = 2)
  inp <- sim_to_epitracer(big, noise = clean, seed = 23)

  se <- q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                               inp$cnv_gr, inp$cancer_genes_gr))
  ct <- q(call_chromothripsis(inp$ecdna_gr, inp$breakpoints_gr,
                              inp$cnv_gr, inp$cancer_genes_gr))
  expect_true(any(se$episomal == "TRUE"))
  expect_false(any(ct$chromothripsis == "TRUE"))
})

test_that("an internal deletion splits the circle and adds one DEL junction", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T17", copies = 2, seed = 30)
  ec <- sim_replicate(ec, rounds = 4)                      # 32 intact circles
  expect_equal(ec$copies, 32)

  before_span <- sum(as.numeric(GenomicRanges::width(ec$fragments)))
  ec <- sim_internal_deletion(ec, start = 55300000, end = 55350000, copies = 1)

  pop <- EpiTracer:::.sim_population(ec)
  expect_length(pop, 2L)
  ## the altered molecule is drawn OUT of the intact pool: total is unchanged
  expect_equal(ec$copies, 32)
  expect_setequal(vapply(pop, function(v) v$copies, numeric(1)), c(31, 1))

  ## the deleted species is one circle broken into two fragments, shorter by
  ## exactly the deleted interval -- this is the regression guard for GRanges
  ## `setdiff` dispatching to element-wise set difference instead of interval
  ## arithmetic, which silently returned the fragment untouched
  delv <- pop[[which(vapply(pop, function(v) v$label, character(1)) == "internal_del")]]
  expect_length(delv$fragments, 2L)
  expect_equal(sum(as.numeric(GenomicRanges::width(delv$fragments))),
               before_span - 50001)

  ## exactly one new junction, deletion-orientation, at the cut
  newj <- ec$junctions[mechanism == "internal_deletion"]
  expect_equal(nrow(newj), 1L)
  expect_equal(newj$svclass, "DEL")
  expect_equal(newj$start1, 55299999)
  expect_equal(newj$start2, 55350001)
})

test_that("a late internal deletion is emitted below the boundary duplication", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T18", copies = 2, seed = 31)
  ec <- sim_replicate(ec, rounds = 4)
  ec <- sim_internal_deletion(ec, start = 55300000, end = 55350000, copies = 1)
  inp <- sim_to_epitracer(ec, noise = sim_noise(0), seed = 31)
  bp <- gr2dt(inp$breakpoints_gr)

  dup  <- bp[svclass == "DUP"][1]
  idel <- bp[svclass == "DEL" & start == 55299999][1]

  ## the boundary rides every circle; the deletion rides only the altered one
  expect_equal(dup$PURPLE_JCN, 32)
  expect_equal(idel$PURPLE_JCN, 1)
  expect_lt(idel$VF, dup$VF)

  ## copy number dips by exactly the deleted species' abundance, and only there
  cn <- gr2dt(inp$cnv_gr)
  expect_equal(cn[start == 55300000]$copyNumber, 32)      # 31 circles + 1 allele
  expect_equal(cn[start == 55350001]$copyNumber, 33)      # 32 circles + 1 allele
})

test_that("replication grows every species and preserves their proportions", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T19", copies = 2, seed = 32)
  ec <- sim_replicate(ec, rounds = 4)
  ec <- sim_internal_deletion(ec, start = 55300000, end = 55350000, copies = 1)
  grown <- sim_replicate(ec, rounds = 3, fold = 2)

  pop <- EpiTracer:::.sim_population(grown)
  expect_equal(grown$copies, 256)
  expect_setequal(vapply(pop, function(v) v$copies, numeric(1)), c(248, 8))

  ## the late deletion has grown but still sits well below the boundary junction
  bp <- gr2dt(sim_to_epitracer(grown, noise = sim_noise(0), seed = 32)$breakpoints_gr)
  expect_equal(bp[svclass == "DUP"][1]$PURPLE_JCN, 256)
  expect_equal(bp[start == 55299999][1]$PURPLE_JCN, 8)
})

test_that("an internal deletion may not remove a preserved oncogene", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T20", copies = 10, seed = 33)
  ## EGFR spans 55019021-55211628; a deletion over it must be refused
  expect_error(sim_internal_deletion(ec, start = 55100000, end = 55150000),
               "preserved oncogene")
  ## and a randomly placed one never lands on it
  for (i in 1:10) {
    d <- sim_internal_deletion(ec, width = 5e4, seed = 200 + i)
    delv <- EpiTracer:::.sim_population(d)
    delv <- delv[[which(vapply(delv, function(v) v$label, character(1)) == "internal_del")]]
    expect_gt(sum(GenomicRanges::countOverlaps(delv$fragments, ec$preserved)), 0L)
  }
})

test_that("a circle carrying a late internal deletion is still called episomal", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T21", copies = 4, seed = 34)
  ec <- sim_replicate(ec, rounds = 4)
  ec <- sim_internal_deletion(ec, width = 5e4, copies = 2, seed = 34)
  ec <- sim_replicate(ec, rounds = 2)
  inp <- sim_to_epitracer(ec, noise = sim_noise(fp_rate = 0, dropout = 0), seed = 34)

  se <- q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                               inp$cnv_gr, inp$cancer_genes_gr))
  ## the boundary DUP still outranks the internal deletion on read support, so
  ## the amplicon reads as a simple excision that has since lost an internal piece
  expect_true(any(se$episomal == "TRUE"))
  expect_true(any(se$duplication_at_boundary_has_highest_VF == "TRUE"))
})

test_that("junction read support is proportionate to junction copy number", {
  clean <- sim_noise(fp_rate = 0, dropout = 0)
  ec <- sim_episome(seed_locus("EGFR"), sample = "T22", copies = 2, seed = 40)
  ec <- sim_replicate(ec, rounds = 4)                                  # 32 copies
  ec <- sim_internal_deletion(ec, start = 55300000, end = 55350000, copies = 1)

  ## expected support is exactly proportionate: 32x the copies, 32x the support
  exact <- gr2dt(sim_to_epitracer(ec, noise = sim_noise(0), seed = 40)$breakpoints_gr)
  dup0 <- exact[svclass == "DUP"][1]; del0 <- exact[start == 55299999][1]
  expect_equal(dup0$VF / dup0$PURPLE_JCN, del0$VF / del0$PURPLE_JCN)
  expect_equal(dup0$VF / del0$VF, dup0$PURPLE_JCN / del0$PURPLE_JCN)

  ## and it stays proportionate on average once noise is on
  draws <- vapply(1:120, function(i) {
    bp <- gr2dt(sim_to_epitracer(ec, noise = clean, seed = 500 + i)$breakpoints_gr)
    c(dup = bp[svclass == "DUP"][1]$VF, del = bp[start == 55299999][1]$VF)
  }, numeric(2))
  per_copy_dup <- mean(draws["dup", ], na.rm = TRUE) / 32
  per_copy_del <- mean(draws["del", ], na.rm = TRUE) / 1
  expect_equal(per_copy_del, per_copy_dup, tolerance = 0.15)

  ## support is a COUNT, so relative noise shrinks as support grows -- the flat
  ## log-normal this replaced gave both junctions the same CV, which hid the
  ## proportionality behind implausible scatter on high-support junctions
  cv <- function(v) stats::sd(v, na.rm = TRUE) / mean(v, na.rm = TRUE)
  expect_lt(cv(draws["dup", ]), cv(draws["del", ]))
  expect_lt(cv(draws["dup", ]), 0.2)
})

test_that("a neutral minority species keeps a frozen share, a fitter one expands", {
  base <- sim_replicate(sim_episome(seed_locus("EGFR"), sample = "T23",
                                    copies = 2, seed = 41), rounds = 4)
  share <- function(x) {
    p <- EpiTracer:::.sim_population(x)
    p[[2]]$copies / x$copies
  }
  mk <- function(fit) sim_internal_deletion(base, start = 55300000, end = 55350000,
                                            copies = 1, fitness = fit)

  ## neutral: both species double together, so the ratio never moves. This is
  ## why a single-molecule deletion sits ~32-fold below the boundary forever.
  neutral <- mk(1)
  expect_equal(share(neutral), 1 / 32)
  expect_equal(share(sim_replicate(neutral, rounds = 6)), 1 / 32)

  ## fitter: the shorter circle outgrows the intact one and takes share
  fitter <- mk(1.3)
  expect_equal(share(fitter), 1 / 32)
  s6 <- share(sim_replicate(fitter, rounds = 6))
  expect_gt(s6, share(fitter))
  expect_equal(s6, share(sim_replicate(fitter, rounds = 6)))   # deterministic

  ## and the gap between the two junctions narrows accordingly
  gap <- function(x) {
    bp <- gr2dt(sim_to_epitracer(x, noise = sim_noise(0), seed = 41)$breakpoints_gr)
    bp[svclass == "DUP"][1]$PURPLE_JCN / bp[start == 55299999][1]$PURPLE_JCN
  }
  expect_gt(gap(sim_replicate(neutral, rounds = 6)),
            gap(sim_replicate(fitter,  rounds = 6)))
})

test_that("the copy-number ceiling rescales species without changing their ratio", {
  base <- sim_replicate(sim_episome(seed_locus("EGFR"), sample = "T24",
                                    copies = 2, seed = 42), rounds = 4)
  ec <- sim_internal_deletion(base, start = 55300000, end = 55350000,
                              copies = 1, fitness = 1.3)
  capped <- sim_replicate(ec, rounds = 8, fold = 2, max_cn = 250)

  expect_lte(summary(capped)$max_cn, 250 + 1e-6)
  ## at the plateau the total is fixed but the fitter species keeps taking share,
  ## so its junction keeps gaining support while the boundary's flattens
  p <- EpiTracer:::.sim_population(capped)
  expect_gt(p[[2]]$copies / capped$copies, 1 / 32)
})

test_that("fission religates the fragments into several separate episomes", {
  ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T25",
                    copies = 90, seed = 50)
  fis <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 3L, seed = 50)

  pop <- EpiTracer:::.sim_population(fis)
  expect_gt(length(pop), 1L)
  expect_true(fis$circular)
  expect_equal(fis$round, 1L)
  expect_true("micronuclear_fission" %in% fis$history)

  ## the parent's copies are split across the surviving daughters
  expect_equal(fis$copies, ec$copies)
  ## every daughter is smaller than the parent, and no fragment is shared
  parent_span <- sum(as.numeric(GenomicRanges::width(ec$fragments)))
  for (v in pop) {
    expect_lt(sum(as.numeric(GenomicRanges::width(v$fragments))), parent_span)
    ## and every surviving daughter carries the oncogene
    expect_gt(sum(GenomicRanges::countOverlaps(v$fragments, ec$preserved)), 0L)
  }
})

test_that("fission destroys the founder boundary duplication", {
  founder_kept <- vapply(1:15, function(i) {
    ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T26",
                      copies = 90, seed = 60 + i)
    fis <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 3L,
                                   seed = 600 + i)
    f <- ec$junctions[origin == "amplicon"]
    nrow(merge(fis$junctions, f[, .(chrom1, start1, chrom2, start2)],
               by = c("chrom1", "start1", "chrom2", "start2"))) > 0L
  }, logical(1))
  ## it survives only if a daughter rejoins those exact two ends, which after
  ## twenty cuts essentially never happens
  expect_false(any(founder_kept))
})

test_that("fission products are no longer called episomal, though they are", {
  ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T27",
                    copies = 90, host_ploidy = 3, seed = 70)
  fis <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 3L, seed = 70)
  clean <- sim_noise(fp_rate = 0, dropout = 0)

  before <- sim_to_epitracer(ec,  noise = clean, seed = 70)
  after  <- sim_to_epitracer(fis, noise = clean, seed = 70)

  se_b <- q(call_simple_excision(before$ecdna_gr, before$breakpoints_gr,
                                 before$cnv_gr, before$cancer_genes_gr))
  se_a <- q(call_simple_excision(after$ecdna_gr, after$breakpoints_gr,
                                 after$cnv_gr, after$cancer_genes_gr))

  ## the parent is recognised; its products are not, even though every surviving
  ## daughter is a circle carrying the oncogene -- the spanning duplication the
  ## caller keys on did not survive the shattering
  expect_true(any(se_b$episomal == "TRUE"))
  expect_false(any(se_a$episomal == "TRUE"))
  expect_true(all(vapply(EpiTracer:::.sim_population(fis),
                         function(v) length(v$fragments) > 0L, logical(1))))
})

test_that("keep_parent leaves the founder circle alongside its daughters", {
  ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T28",
                    copies = 90, seed = 80)
  fis <- sim_shatter_to_episomes(ec, n_breaks = 15L, n_circles = 2L,
                                 keep_parent = TRUE, copies = 3, seed = 80)

  pop <- EpiTracer:::.sim_population(fis)
  labels <- vapply(pop, function(v) v$label, character(1))
  expect_true("founder" %in% labels)
  ## the parent keeps its copies; the daughters arise as a minority
  expect_equal(pop[[which(labels == "founder")]]$copies, 90)
  expect_gt(fis$copies, ec$copies)

  ## and because the parent survives, so does its boundary duplication
  bp <- gr2dt(sim_to_epitracer(fis, noise = sim_noise(0), seed = 80)$breakpoints_gr)
  expect_gt(max(bp[svclass == "DUP"]$PURPLE_JCN), 50)
})

test_that("fission requires enough fragments to make the requested circles", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T29", copies = 40, seed = 90)
  expect_error(sim_shatter_to_episomes(ec, n_breaks = 1L, n_circles = 6L, seed = 90),
               "too few fragments")
})

test_that("neutral division does not raise ecDNA copy number, however long it runs", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T30", copies = 20, seed = 100)

  ## random segregation is copy-number-neutral in expectation: replication
  ## doubles, division halves. Eighty generations move the mean nowhere.
  mean_over <- function(g) mean(vapply(1:5, function(i)
    sim_segregate(ec, generations = g, selection = 0, addicted = FALSE,
                  seed = i)$copies, numeric(1)))
  expect_equal(mean_over(10L), 20, tolerance = 0.15)
  expect_equal(mean_over(80L), 20, tolerance = 0.20)

  ## what division DOES generate is cell-to-cell spread, which bulk averages away
  short <- attr(sim_segregate(ec, generations = 5L,  selection = 0,
                              addicted = FALSE, seed = 1), "cells")
  long  <- attr(sim_segregate(ec, generations = 40L, selection = 0,
                              addicted = FALSE, seed = 1), "cells")
  cv <- function(v) stats::sd(v) / mean(v)
  expect_gt(cv(long), cv(short))
  expect_gt(mean(long == 0), mean(short == 0))     # a tail loses the amplicon
})

test_that("only selection raises copy number, and read support follows it", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T31", copies = 20, seed = 101)
  clean <- sim_noise(fp_rate = 0, dropout = 0)

  neutral <- sim_segregate(ec, generations = 25L, selection = 0,
                           addicted = FALSE, seed = 2)
  selected <- sim_segregate(ec, generations = 25L, selection = 0.1, seed = 2)
  expect_gt(selected$copies, 2 * neutral$copies)

  ## boundary read support tracks the copy level, not the number of divisions
  vf <- function(x) gr2dt(sim_to_epitracer(x, noise = clean,
                                           seed = 2)$breakpoints_gr)[svclass == "DUP"][1]$VF
  expect_gt(vf(selected), vf(neutral))
  expect_equal(vf(selected) / vf(neutral), selected$copies / neutral$copies,
               tolerance = 0.15)

  ## and longer under selection means more copies
  expect_gt(sim_segregate(ec, generations = 40L, selection = 0.05, seed = 3)$copies,
            sim_segregate(ec, generations = 10L, selection = 0.05, seed = 3)$copies)
})

test_that("sim_segregate keeps species proportions and respects max_cn", {
  ec <- sim_replicate(sim_episome(seed_locus("EGFR"), sample = "T32",
                                  copies = 4, seed = 102), rounds = 3)
  ec <- sim_internal_deletion(ec, start = 55300000, end = 55350000, copies = 4)
  before <- vapply(EpiTracer:::.sim_population(ec), function(v) v$copies, numeric(1))

  seg <- sim_segregate(ec, generations = 15L, selection = 0.1, max_cn = 120, seed = 4)
  after <- vapply(EpiTracer:::.sim_population(seg), function(v) v$copies, numeric(1))

  expect_equal(after / sum(after), before / sum(before))   # proportions preserved
  expect_lte(summary(seg)$max_cn, 120 + 1e-6)
  expect_type(attr(seg, "cells"), "integer")
})

test_that("reconfiguration preserves the founder junction as the closing junction", {
  ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T33",
                    copies = 90, seed = 110)
  founder <- ec$junctions[origin == "amplicon"][1]

  rc <- sim_reconfigure(ec, n_breaks = 12L, copies = NULL, seed = 110)
  expect_true(attr(rc, "preserved_founder"))

  ## the circle is now a mosaic, but the founder junction is still traversed
  pop <- EpiTracer:::.sim_population(rc)
  expect_gt(length(pop[[1]]$fragments), 1L)
  j <- EpiTracer:::.sim_junctions_from_fragments(pop[[1]]$fragments, TRUE, "", 0L)
  expect_true(any(j$start1 == founder$start1 & j$start2 == founder$start2))

  ## and preserve = "none" does not protect it
  free <- sim_reconfigure(ec, n_breaks = 12L, copies = NULL,
                          preserve = "none", seed = 110)
  expect_false(isTRUE(attr(free, "preserved_founder")))
})

test_that("the founder junction stays clonal with NO unshattered molecule left", {
  ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "T34",
                    copies = 90, host_ploidy = 3, seed = 111)
  founder <- ec$junctions[origin == "amplicon"][1]

  ## round 1 consumes the entire founder species; later rounds split lineages off
  x <- sim_reconfigure(ec, n_breaks = 12L, copies = NULL, seed = 21, label = "L1")
  for (i in 2:5) {
    dom <- max(vapply(EpiTracer:::.sim_population(x), function(v) v$copies, numeric(1)))
    x <- sim_reconfigure(x, n_breaks = 12L, copies = dom * 0.5,
                         seed = 20 + i, label = paste0("L", i))
  }

  ## nothing in the population is an unshattered circle any more
  pop <- EpiTracer:::.sim_population(x)
  expect_true(all(vapply(pop, function(v) length(v$fragments) > 1L, logical(1))))
  expect_gt(length(pop), 1L)

  ## yet the founder junction rides EVERY molecule, so it is emitted at the full
  ## amplicon copy number and tops the read-support hierarchy
  bp <- gr2dt(sim_to_epitracer(x, noise = sim_noise(fp_rate = 0, dropout = 0),
                               seed = 99)$breakpoints_gr)
  bp <- bp[PURPLE_CN > 9][!duplicated(event)]
  isf <- bp$start %in% c(founder$start1, founder$start2) & bp$svclass == "DUP"
  expect_true(any(isf))
  expect_equal(bp[isf]$PURPLE_JCN[1], x$copies)
  expect_equal(max(bp$VF), bp[isf]$VF[1])
})

test_that("reconfiguration falls back gracefully when the junction is gone", {
  ec <- sim_episome(seed_locus("EGFR"), sample = "T35", copies = 40, seed = 112)
  ## a junction whose ends are not fragment boundaries cannot be preserved
  bogus <- data.frame(start1 = 1L, start2 = 2L)
  rc <- sim_reconfigure(ec, n_breaks = 10L, copies = NULL,
                        preserve = bogus, seed = 112)
  expect_false(attr(rc, "preserved_founder"))
  expect_gt(length(EpiTracer:::.sim_population(rc)[[1]]$fragments), 1L)
})
