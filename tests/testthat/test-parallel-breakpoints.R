test_that("find_parallel_breakpoints pairs same-orientation, near, distinct-junction breakends", {
  # J1@100 and J2@110 are same strand (+), 10 bp apart, different junctions -> a pair.
  # J3@500 is the opposite strand and far -> no pair.
  p <- find_parallel_breakpoints(pos = c(100, 110, 500),
                                 strand = c("+", "+", "-"),
                                 event = c("J1", "J2", "J3"))
  expect_equal(nrow(p), 1L)
  expect_equal(p$pos1, 100); expect_equal(p$pos2, 110)
  expect_equal(p$strand, "+"); expect_equal(p$dist, 10)
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
