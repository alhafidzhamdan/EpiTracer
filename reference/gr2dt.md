# Convert a GRanges to a data.table

Internal replacement for `gUtils::gr2dt()` for the columns EpiTracer
needs: `seqnames`, `start`, `end`, `width`, `strand` (character
`seqnames`/`strand`) plus every metadata column.

## Usage

``` r
gr2dt(x)
```

## Arguments

- x:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html).
