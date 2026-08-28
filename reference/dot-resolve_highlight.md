# Which SV rows match a highlight set

Which SV rows match a highlight set

## Usage

``` r
.resolve_highlight(sv, highlight_events, highlight_id_col = NULL)
```

## Arguments

- sv:

  A per-sample SV `data.frame`/`data.table`.

- highlight_events:

  Character vector of event ids to emphasise (or `NULL`).

- highlight_id_col:

  Name of the id column in `sv`; `NULL` auto-detects `name` then
  `event`.

## Value

A logical vector, one per row of `sv` (all `FALSE` when
`highlight_events` is `NULL`/empty).
