# Detect focal-amplicon seeds from copy number alone

Derives candidate focal-amplicon regions ("seeds") directly from
allele-specific copy-number segments, so
[`call_episomal_ecdna()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_ecdna.md)
can run **without** an AmpliconArchitect amplicon catalogue. Per sample,
segments with `copyNumber > min_cn_ratio * ploidy` are merged across
gaps up to `gap` bp and regions narrower than `min_width` are dropped;
each surviving region is labelled with an `ID` and `WGS_ID`, matching
the `ecdna_gr` contract.

## Usage

``` r
detect_amplicon_seeds(cnv_gr, min_cn_ratio = 3, gap = 1e+06, min_width = 1e+05)
```

## Arguments

- cnv_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of allele-specific copy-number segments with metadata columns
  `sample`, `copyNumber`, `ploidy`.

- min_cn_ratio:

  Numeric; a segment is amplified where
  `copyNumber > min_cn_ratio * ploidy` (default `3`).

- gap:

  Integer; merge amplified segments separated by at most this many bp
  into one seed (default `1e6`).

- min_width:

  Integer; drop seeds narrower than this (default `1e5`).

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
of amplicon seeds with metadata columns `ID` and `WGS_ID` (empty if no
amplification is present).

## See also

[`call_episomal_ecdna()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_ecdna.md)

## Examples

``` r
seeds <- detect_amplicon_seeds(ex_caller_inputs$cnv_gr)
seeds
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames            ranges strand |      WGS_ID             ID
#>          <Rle>         <IRanges>  <Rle> | <character>    <character>
#>   [1]     chr7 55000000-55500000      * |   EXAMPLE01 EXAMPLE01_amp1
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
