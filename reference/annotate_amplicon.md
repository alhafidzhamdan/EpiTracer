# Per-amplicon breakpoint annotation (shared set-up)

Internal helper shared by
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
and the standalone mechanism callers
([`call_brf()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_brf.md),
[`call_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_micronucleation.md),
[`call_bfb()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.md)).
For one amplicon it builds the structural-variant breakpoint table
falling within (or just outside, by `ext`) the amplicon, annotated with
oncogene overlap and allele-specific copy number, and returns it
together with the sample's full breakpoint set, copy-number segments,
ploidy and the amplicon ranges.

## Usage

``` r
annotate_amplicon(
  this_amplicon_id,
  ecdna_gr,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07
)
```

## Arguments

- this_amplicon_id:

  Character scalar; one value of `ecdna_gr$ID`.

- ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr, ext:

  See
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md).

## Value

A list with `bp` (annotated per-breakend data.table for the amplicon),
`sample_bp` (all sample breakends, a GRanges), `cnv` (sample CN
GRanges), `ploidy` (numeric) and `amplicon_gr` (the amplicon ranges); or
`NULL` if the amplicon has no breakpoints in range.
