# GenomicRanges overlap-annotation operator

`x \%$\% y` returns the GRanges `x` with the metadata columns of the
first range in `y` that each range of `x` overlaps (`NA` where there is
no overlap). A small, dependency-free replacement for the gUtils
operator of the same name, covering the annotate-by-overlap use inside
EpiTracer (where the ranges of `y` — amplicons, genes, copy-number
segments — do not overlap one another).

## Usage

``` r
x %$% y
```

## Arguments

- x:

  A GRanges to be annotated.

- y:

  A GRanges whose metadata is transferred onto `x`.

## Value

`x` with metadata columns added from overlapping ranges in `y`.
