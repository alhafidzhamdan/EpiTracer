#' Adjacent parallel breakpoints (breakage-replication/fusion hallmark)
#'
#' Finds pairs of structural-variant breakends that are the SAME orientation,
#' within `max_dist` bp of each other, and from DISTINCT junctions -- "adjacent
#' parallel breakpoints", the hallmark footprint of the breakage-replication/
#' fusion process, in which a broken DNA end is replicated before repair to yield
#' two same-orientation sister ends a short distance apart (Mendez-Dorantes,
#' Zhang, Burns & Pellman, *Nat Genet* 2026; doi:10.1038/s41588-025-02434-5).
#'
#' Used both by [call_simple_excision()] (to set the `brf` annotation and
#' `n_parallel_pairs` count) and by plotting code to mark the pairs.
#'
#' @param pos Numeric vector of breakend positions (bp).
#' @param strand Character vector of breakend orientations (`"+"`/`"-"`), the same
#'   length as `pos`. Breakends with `NA` strand are ignored.
#' @param event Character vector of junction identifiers, the same length as
#'   `pos`; two breakends of the same junction never form a pair.
#' @param max_dist Integer; the maximum distance (bp) between the two breakends of
#'   a pair (default `2e4`, the 20 kb heuristic from the source paper).
#'
#' @return A data.frame with one row per adjacent parallel breakpoint pair and
#'   columns `pos1`, `pos2` (with `pos1 <= pos2`), `strand`, and `dist`
#'   (`pos2 - pos1`). Zero rows if no pair qualifies.
#' @examples
#' find_parallel_breakpoints(pos = c(100, 110, 500), strand = c("+", "+", "-"),
#'                           event = c("J1", "J2", "J3"))
#' @export
find_parallel_breakpoints <- function(pos, strand, event, max_dist = 2e4) {
  stopifnot(length(pos) == length(strand), length(pos) == length(event))
  empty <- data.frame(pos1 = numeric(0), pos2 = numeric(0),
                      strand = character(0), dist = numeric(0),
                      stringsAsFactors = FALSE)
  ok <- !is.na(pos) & !is.na(strand)
  pos <- as.numeric(pos)[ok]; strand <- as.character(strand)[ok]
  event <- as.character(event)[ok]
  n <- length(pos)
  if (n < 2L) return(empty)
  o <- order(pos); pos <- pos[o]; strand <- strand[o]; event <- event[o]
  rows <- list()
  for (a in seq_len(n - 1L))
    for (b in (a + 1L):n) {
      if (pos[b] - pos[a] > max_dist) break
      if (strand[b] == strand[a] && event[b] != event[a])
        rows[[length(rows) + 1L]] <- data.frame(
          pos1 = pos[a], pos2 = pos[b], strand = strand[a],
          dist = pos[b] - pos[a], stringsAsFactors = FALSE)
    }
  if (length(rows)) do.call(rbind, rows) else empty
}
