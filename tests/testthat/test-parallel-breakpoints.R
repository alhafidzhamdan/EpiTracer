test_that("find_parallel_breakpoints pairs same-orientation, near, distinct-junction breakends", {
  # J1@100 and J2@110 are same strand (+), 10 bp apart, different junctions -> a pair.
  # J3 is the opposite strand and far enough away that it neither forms a +/-
  # insertion adjacency with them nor makes independent origin plausible.
  p <- find_parallel_breakpoints(pos = c(100, 110, 5e5),
                                 strand = c("+", "+", "-"),
                                 event = c("J1", "J2", "J3"))
  expect_equal(nrow(p), 1L)
  expect_equal(p$pos1, 100); expect_equal(p$pos2, 110)
  expect_equal(p$strand, "+"); expect_equal(p$dist, 10)
})

## The three criteria the source paper applies (Zhang, Mendez-Dorantes, Burns &
## Pellman, Nat Genet 58, 88-99, 2026). Only the second was implemented before,
## which made call_brf() fire on essentially any junction-dense amplicon.

test_that("breakends in a +/- insertion adjacency are excluded first", {
  # 100(+) and 400(-) are a +/- adjacency within 20 kb: an insertion or an
  # overlapping junction, a different process, so both breakends are dropped --
  # taking the 100/110 parallel pair with them.
  p <- find_parallel_breakpoints(c(100, 110, 400), c("+", "+", "-"),
                                 c("J1", "J2", "J3"))
  expect_equal(nrow(p), 0L)
  # switching the exclusion off recovers the pair
  p2 <- find_parallel_breakpoints(c(100, 110, 400), c("+", "+", "-"),
                                  c("J1", "J2", "J3"),
                                  exclude_insertion_adjacency = FALSE,
                                  max_indep_p = Inf)
  expect_equal(nrow(p2), 1L)
})

test_that("identical breakends are not a parallel pair", {
  # one DNA end joined to two partners reports twice at the SAME coordinate; the
  # paper's sister ends are adjacent but NON-identical, separated by resection
  p <- find_parallel_breakpoints(c(100, 100, 5e5), c("+", "+", "-"),
                                 c("J1", "J2", "J3"))
  expect_equal(nrow(p), 0L)
  # a resection tract of even 1 bp qualifies
  expect_equal(nrow(find_parallel_breakpoints(c(100, 101, 5e5), c("+", "+", "-"),
                                              c("J1", "J2", "J3"))), 1L)
})

test_that("the independence test rejects chance pairs in dense regions", {
  # a +/+ pair 2 kb apart with the nearest (-) breakend 200 kb to the right:
  # independent origin is unlikely (ratio 0.01), so it is kept
  far <- find_parallel_breakpoints(c(100, 2100, 2e5), c("+", "+", "-"),
                                   c("J1", "J2", "J3"))
  expect_equal(nrow(far), 1L)
  expect_lt(far$indep_p, 0.05)

  # The same pair with a (-) breakend only 10 kb to the right: the ancestral
  # segment can be no longer than that, so the pair could easily be coincidence.
  # The insertion-adjacency exclusion is switched off here so that the ONLY thing
  # separating the two calls below is the independence test itself.
  near <- find_parallel_breakpoints(c(100, 2100, 12100), c("+", "+", "-"),
                                    c("J1", "J2", "J3"),
                                    exclude_insertion_adjacency = FALSE)
  expect_equal(nrow(near), 0L)
  # ... and it reappears when the independence test is switched off
  loose <- find_parallel_breakpoints(c(100, 2100, 12100), c("+", "+", "-"),
                                     c("J1", "J2", "J3"),
                                     exclude_insertion_adjacency = FALSE,
                                     max_indep_p = Inf)
  expect_equal(nrow(loose), 1L)
  expect_equal(loose$indep_p, 0.2)   # 2 kb apart / 10 kb to the nearest (-)
})

test_that("-/- pairs look left for the nearest opposite breakend", {
  # for -/- the relevant neighbour is the first (+) breakpoint to the LEFT
  p <- find_parallel_breakpoints(c(100, 2e5, 2.02e5), c("+", "-", "-"),
                                 c("J1", "J2", "J3"))
  expect_equal(nrow(p), 1L)
  expect_equal(p$strand, "-")
  expect_equal(p$opp_dist, 2e5 - 100)
})

test_that("find_parallel_breakpoints excludes opposite strand, too-far, and same-junction pairs", {
  # opposite strand within range -> not parallel
  expect_equal(nrow(find_parallel_breakpoints(c(100, 110), c("+", "-"), c("A", "B"))), 0L)
  # same strand but > max_dist apart -> not adjacent
  expect_equal(nrow(find_parallel_breakpoints(c(100, 40000), c("+", "+"), c("A", "B"))), 0L)
  # same strand, near, but SAME junction (both breakends of one SV) -> not counted
  expect_equal(nrow(find_parallel_breakpoints(c(100, 110), c("+", "+"), c("A", "A"))), 0L)
  # max_dist is honoured
  expect_equal(nrow(find_parallel_breakpoints(c(100, 30000), c("+", "+"), c("A", "B"),
                                              max_dist = 5e4)), 1L)
})

test_that("find_parallel_breakpoints ignores NA strand and handles < 2 breakends", {
  expect_equal(nrow(find_parallel_breakpoints(numeric(0), character(0), character(0))), 0L)
  expect_equal(nrow(find_parallel_breakpoints(100, "+", "A")), 0L)
  expect_equal(nrow(find_parallel_breakpoints(c(100, 110), c("+", NA), c("A", "B"))), 0L)
})
