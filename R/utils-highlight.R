## Shared SV-highlight helper used by plot_sv_linear(), plot_sv_reconstruction()
## and plot_sv_circos(). Given a per-sample SV table and a set of event ids, it
## returns a logical vector marking which SV rows to emphasise. The event id is
## matched against `highlight_id_col`; when that is NULL it auto-detects a `name`
## then `event` column (the usual BEDPE / breakpoint identifiers). This is what
## lets any event set -- a chromoplexy cycle's junctions, a TBA boundary, a set
## of fold-backs -- be drawn bold on top of the full rearrangement background.

#' Which SV rows match a highlight set
#'
#' @param sv A per-sample SV `data.frame`/`data.table`.
#' @param highlight_events Character vector of event ids to emphasise (or `NULL`).
#' @param highlight_id_col Name of the id column in `sv`; `NULL` auto-detects
#'   `name` then `event`.
#' @return A logical vector, one per row of `sv` (all `FALSE` when
#'   `highlight_events` is `NULL`/empty).
#' @keywords internal
.resolve_highlight <- function(sv, highlight_events, highlight_id_col = NULL) {
  n <- nrow(sv)
  if (is.null(highlight_events) || !length(highlight_events)) return(rep(FALSE, n))
  col <- highlight_id_col
  if (is.null(col))
    col <- if ("name" %in% names(sv)) "name" else if ("event" %in% names(sv)) "event" else NA_character_
  if (is.na(col) || !col %in% names(sv))
    stop("`highlight_events` was supplied but no id column was found in sv_data; ",
         "set `highlight_id_col` to a column holding the SV identifiers.", call. = FALSE)
  as.character(sv[[col]]) %in% as.character(highlight_events)
}
