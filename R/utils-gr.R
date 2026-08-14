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
