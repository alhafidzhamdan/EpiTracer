## Small, dependency-free GRanges helpers replacing the two gUtils functions
## EpiTracer used (gUtils::gr2dt and the `%$%` operator), so the package installs
## from CRAN/Bioconductor alone with no GitHub-only remote.

#' Convert a GRanges to a data.table
#'
#' Internal replacement for `gUtils::gr2dt()` for the columns EpiTracer needs:
#' `seqnames`, `start`, `end`, `width`, `strand` (character `seqnames`/`strand`)
#' plus every metadata column.
#'
#' @param x A [GenomicRanges::GRanges].
#' @return A [data.table::data.table].
#' @keywords internal
#' @importFrom data.table as.data.table :=
gr2dt <- function(x) {
  dt <- data.table::as.data.table(as.data.frame(x, row.names = NULL))
  if ("seqnames" %in% names(dt)) dt[, seqnames := as.character(seqnames)]
  if ("strand"   %in% names(dt)) dt[, strand   := as.character(strand)]
  dt[]
}

#' Convert a data.frame/data.table to GRanges
#'
#' Internal replacement for `regioneR::toGRanges()`: builds a
#' [GenomicRanges::GRanges] from a table with `seqnames`, `start`, `end`
#' columns, keeping every other column as metadata. Dependency-free (uses only
#' GenomicRanges), so EpiTracer needs neither regioneR nor a GitHub remote.
#'
#' @param x A `data.frame`/`data.table` with `seqnames`, `start`, `end` columns.
#' @return A [GenomicRanges::GRanges].
#' @keywords internal
#' @importFrom GenomicRanges makeGRangesFromDataFrame
to_granges <- function(x) {
  x <- as.data.frame(x)
  ## makeGRangesFromDataFrame() errors on a 0-row table; return an empty GRanges
  ## that still carries the metadata columns, so callers can keep piping/selecting
  ## (e.g. an amplicon with no in-region breakpoints classifies as not-episomal
  ## rather than crashing the run).
  if (nrow(x) == 0L) {
    gr <- GenomicRanges::GRanges()
    extra <- setdiff(names(x), c("seqnames", "start", "end", "width", "strand"))
    if (length(extra)) S4Vectors::mcols(gr) <- S4Vectors::DataFrame(x[, extra, drop = FALSE])
    return(gr)
  }
  GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
}
