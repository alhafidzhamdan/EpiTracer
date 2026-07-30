#' Round a value up to a "nice" axis top
#'
#' Rounds `x` up so an axis top lands on a round number carrying `digits`
#' significant figures (e.g. `nice_ceiling(103) == 110`, `nice_ceiling(1.005) ==
#' 1.1`). Non-finite or non-positive `x` returns `1`. Shared by the copy-number /
#' read-support axes of [plot_sv_linear()] and [plot_sv_reconstruction()].
#'
#' @param x Numeric value to round up.
#' @param digits Number of significant figures in the result (default `2`).
#' @return A single numeric axis top.
#' @keywords internal
#' @noRd
nice_ceiling <- function(x, digits = 2) {
  if (!is.finite(x) || x <= 0) return(1)
  mag <- 10^(floor(log10(x)) - (digits - 1)); ceiling(x / mag) * mag
}
