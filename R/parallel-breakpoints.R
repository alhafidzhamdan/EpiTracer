#' Adjacent parallel breakpoints (breakage-replication/fusion hallmark)
#'
#' Finds **adjacent parallel breakpoints** --- the hallmark footprint of the
#' breakage-replication/fusion (BRF) process, in which a staggered DNA end is
#' replicated before repair, so that a single ancestral end yields two
#' same-orientation sister ends a short distance apart (Zhang, Mendez-Dorantes,
#' Burns & Pellman, *Nat Genet* **58**, 88-99, 2026;
#' doi:10.1038/s41588-025-02434-5).
#'
#' The source paper's procedure has three parts, all applied here:
#'
#' 1. **Insertion / overlap exclusion.** Breakends that participate in a `+/-`
#'    adjacency within `max_dist` (an insertion or an overlapping-segment
#'    junction) are removed *first*. The paper is explicit that this is necessary,
#'    because the breakends of adjacent insertions also fall inside the distance
#'    threshold but arise by a different process.
#' 2. **Adjacency and orientation.** From what remains, pairs of breakends with
#'    the SAME orientation (`+/+` or `-/-`), between `min_dist` and `max_dist` bp
#'    apart, and from DISTINCT junctions are collected. The paper's sister ends
#'    are "adjacent but **non-identical**", separated by the resection tract
#'    (429 bp and 2,059 bp in its worked examples), so breakends reported at the
#'    *same* coordinate are one DNA end joined to two partners, not two
#'    replicated ends, and are excluded by `min_dist`.
#' 3. **Independence test.** Two same-orientation breakends can also land near
#'    each other by chance. The paper bounds that probability by the ratio of the
#'    distance between the parallel breakends to the size of the ancestral
#'    segment they would have come from; that segment must be at least as long as
#'    the distance out to the nearest OPPOSITE-orientation breakend --- to the
#'    right of a `+/+` pair, to the left of a `-/-` pair. A pair is kept when
#'    `dist / opp_dist <= max_indep_p`.
#'
#' Step 3 is what separates a genuine BRF footprint from coincidence, and it
#' matters most exactly where it is easiest to get wrong: inside a junction-dense
#' amplicon, same-orientation breakends fall within 20 kb of each other
#' constantly by chance, but so does the nearest opposite-orientation breakend,
#' so the ratio stays large and the pair is correctly rejected. Applying steps 1
#' and 2 alone makes BRF fire on essentially any heavily rearranged amplicon.
#'
#' @param pos Numeric vector of breakend positions (bp).
#' @param strand Character vector of breakend orientations (`"+"`/`"-"`), the same
#'   length as `pos`. Breakends with `NA` strand are ignored.
#' @param event Character vector of junction identifiers, the same length as
#'   `pos`; two breakends of the same junction never form a pair.
#' @param max_dist Integer; the maximum distance (bp) between the two breakends of
#'   a pair (default `2e4`, the 20 kb heuristic from the source paper).
#' @param min_dist Integer; the minimum distance (bp) between them (default `1`,
#'   i.e. the two breakends must be non-identical, as the source paper requires).
#' @param max_indep_p Numeric; keep a pair only when its estimated probability of
#'   having arisen independently, `dist / opp_dist`, is at most this (default
#'   `0.05`, the threshold used in the source paper). Set to `NULL` or `Inf` to
#'   skip the independence test and return every adjacent same-orientation pair.
#' @param exclude_insertion_adjacency Logical; drop breakends involved in a `+/-`
#'   adjacency within `max_dist` before pairing (default `TRUE`, as in the source
#'   paper). Set `FALSE` for the geometric pair-finding alone.
#'
#' @return A data.frame with one row per adjacent parallel breakpoint pair and
#'   columns `pos1`, `pos2` (with `pos1 <= pos2`), `strand`, `dist`
#'   (`pos2 - pos1`), `opp_dist` (distance to the nearest opposite-orientation
#'   breakend on the relevant side; `Inf` if there is none) and `indep_p`
#'   (`dist / opp_dist`, the bound on the chance of independent origin). Zero
#'   rows if no pair qualifies.
#' @seealso [call_brf()], [call_simple_excision()]
#' @examples
#' ## two + breakends 10 bp apart from distinct junctions, with the nearest
#' ## - breakend far away: a confident BRF footprint
#' find_parallel_breakpoints(pos = c(100, 110, 5e5), strand = c("+", "+", "-"),
#'                           event = c("J1", "J2", "J3"))
#'
#' ## the same pair inside a junction-dense region is rejected: the nearest
#' ## opposite-orientation breakend is close, so independent origin is likely
#' find_parallel_breakpoints(pos = c(100, 110, 300), strand = c("+", "+", "-"),
#'                           event = c("J1", "J2", "J3"))
#' @export
find_parallel_breakpoints <- function(pos, strand, event, max_dist = 2e4,
                                      min_dist = 1, max_indep_p = 0.05,
                                      exclude_insertion_adjacency = TRUE) {
  stopifnot(length(pos) == length(strand), length(pos) == length(event))
  empty <- data.frame(pos1 = numeric(0), pos2 = numeric(0),
                      strand = character(0), dist = numeric(0),
                      opp_dist = numeric(0), indep_p = numeric(0),
                      stringsAsFactors = FALSE)
  if (is.null(max_indep_p)) max_indep_p <- Inf

  ok <- !is.na(pos) & !is.na(strand)
  pos <- as.numeric(pos)[ok]; strand <- as.character(strand)[ok]
  event <- as.character(event)[ok]
  if (length(pos) < 2L) return(empty)
  o <- order(pos); pos <- pos[o]; strand <- strand[o]; event <- event[o]

  ## (1) drop breakends that form a +/- adjacency -- insertions and overlapping
  ## segments, which fall inside the same distance threshold but are a different
  ## process from the replication of a single end.
  if (isTRUE(exclude_insertion_adjacency)) {
    drop <- logical(length(pos))
    for (a in seq_len(length(pos) - 1L))
      for (b in (a + 1L):length(pos)) {
        if (pos[b] - pos[a] > max_dist) break
        if (strand[a] == "+" && strand[b] == "-") drop[c(a, b)] <- TRUE
      }
    pos <- pos[!drop]; strand <- strand[!drop]; event <- event[!drop]
    if (length(pos) < 2L) return(empty)
  }

  ## (2) same orientation, within max_dist, distinct junctions
  n <- length(pos)
  rows <- list()
  for (a in seq_len(n - 1L))
    for (b in (a + 1L):n) {
      if (pos[b] - pos[a] > max_dist) break
      if (strand[b] == strand[a] && event[b] != event[a] &&
          (pos[b] - pos[a]) >= min_dist) {
        ## (3) distance out to the nearest opposite-orientation breakend: to the
        ## right for a +/+ pair, to the left for a -/- pair. Measured from the
        ## near edge of the pair, the conservative choice (smallest denominator).
        opp <- if (strand[a] == "+") {
          cand <- pos[strand == "-" & pos > pos[b]]
          if (length(cand)) min(cand) - pos[b] else Inf
        } else {
          cand <- pos[strand == "+" & pos < pos[a]]
          if (length(cand)) pos[a] - max(cand) else Inf
        }
        d <- pos[b] - pos[a]
        p <- if (is.infinite(opp)) 0 else d / opp
        if (p <= max_indep_p)
          rows[[length(rows) + 1L]] <- data.frame(
            pos1 = pos[a], pos2 = pos[b], strand = strand[a], dist = d,
            opp_dist = opp, indep_p = p, stringsAsFactors = FALSE)
      }
    }
  if (length(rows)) do.call(rbind, rows) else empty
}
