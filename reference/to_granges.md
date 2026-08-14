# Convert a data.frame/data.table to GRanges

Internal replacement for `regioneR::toGRanges()`: builds a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
from a table with `seqnames`, `start`, `end` columns, keeping every
other column as metadata. Dependency-free (uses only GenomicRanges), so
EpiTracer needs neither regioneR nor a GitHub remote.

## Usage

``` r
to_granges(x)
```

## Arguments

- x:

  A `data.frame`/`data.table` with `seqnames`, `start`, `end` columns.

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html).
