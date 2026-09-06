# Derive ecDNA intervals from an AmpliconArchitect cycles file

Generalises
[`call_episomal_from_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_from_ampliconsuite.md)
to **any** AmpliconArchitect output, including bare AA runs without
AmpliconClassifier. The ecDNA extent is taken from the reconstructed
cycle: the segments of the dominant (highest `Copy_count`) cycle, or of
all cycles. The AA "source" vertex (segment id `0`) is dropped. Returned
coordinates use the exclusive-end convention of the AmpliconClassifier
BED, so the result plugs straight into the `intervals` argument of
[`call_episomal_from_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_from_ampliconsuite.md).

## Usage

``` r
aa_cycle_intervals(cycles_file, cycle = c("dominant", "all"))
```

## Arguments

- cycles_file:

  Path to an AmpliconArchitect `*_cycles.txt`.

- cycle:

  Which cycle(s) to use: `"dominant"` (default, highest copy count) or
  `"all"`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of `chr`, `start`, `end` (exclusive end), or `NULL` if no usable
cycle/segments are found.

## See also

[`call_episomal_from_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_from_ampliconsuite.md),
[`read_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_ampliconsuite.md)
