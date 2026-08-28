#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' GenomicRanges overlap-annotation operator
#'
#' `x \%$\% y` returns the GRanges `x` with the metadata columns of the first
#' range in `y` that each range of `x` overlaps (`NA` where there is no overlap).
#' A small, dependency-free replacement for the gUtils operator of the same name,
#' covering the annotate-by-overlap use inside EpiTracer (where the ranges of `y`
#' — amplicons, genes, copy-number segments — do not overlap one another).
#'
#' @name %$%
#' @rdname gr-annotate
#' @keywords internal
#' @param x A GRanges to be annotated.
#' @param y A GRanges whose metadata is transferred onto `x`.
#' @return `x` with metadata columns added from overlapping ranges in `y`.
`%$%` <- function(x, y) {
  hit   <- GenomicRanges::findOverlaps(x, y, select = "first")
  ycols <- S4Vectors::mcols(y)
  if (!is.null(ycols) && ncol(ycols) > 0L) {
    add <- ycols[hit, , drop = FALSE]     # NA rows where there is no overlap
    mc  <- S4Vectors::mcols(x)
    if (is.null(mc)) mc <- S4Vectors::DataFrame(row.names = seq_along(x))
    for (nm in colnames(add)) mc[[nm]] <- add[[nm]]
    S4Vectors::mcols(x) <- mc
  }
  x
}

## Null-coalescing default, used by the simulation functions to fall back to a
## derived value when an argument is left NULL.
`%||%` <- function(x, y) if (is.null(x)) y else x
