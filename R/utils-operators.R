#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' GenomicRanges overlap-annotation operator
#'
#' Re-exported from \pkg{gUtils}. `x \%$\% y` annotates the ranges in `x`
#' with the metadata columns of `y` by genomic overlap. See
#' \code{gUtils::\link[gUtils]{\%$\%}} for details.
#'
#' @name %$%
#' @rdname gr-annotate
#' @keywords internal
#' @importFrom gUtils %$%
#' @usage x \%$\% y
#' @param x A GRanges to be annotated.
#' @param y A GRanges whose metadata is transferred onto `x`.
#' @return `x` with metadata columns added from overlapping ranges in `y`.
NULL
