# Build EpiTracer caller inputs from copy-number and SV tables

Assembles the `cnv_gr` and `breakpoints_gr`
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
that the amplicon-mechanism callers
([`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md),
[`call_brf()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_brf.md),
etc.) take, from the tables returned by
[`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md)
/
[`read_sv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_sv.md).
If the CN and SV tables use different sample labels (common with PURPLE
file naming), pass `sample` to force a single name; the SVs are then
taken from the whole SV table.

## Usage

``` r
prepare_amplicon_inputs(
  cnv,
  sv,
  sample = NULL,
  ploidy_mode = c("per_chr", "global")
)
```

## Arguments

- cnv, sv:

  Copy-number and SV tables (from
  [`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md)
  /
  [`read_sv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_sv.md)).

- sample:

  Sample to build for. Defaults to the first CN sample.

- ploidy_mode:

  `"per_chr"` (default) sets each breakend/segment ploidy to the local
  per-chromosome baseline (recovers amplicons on polysomic chromosomes);
  `"global"` uses the CN table's own ploidy.

## Value

A list with `cnv_gr` and `breakpoints_gr`.

## See also

[`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md),
[`read_sv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_sv.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
