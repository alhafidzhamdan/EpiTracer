# A seed locus for simulation

Convenience wrapper that resolves an oncogene symbol to the amplicon
seed and preserved-region pair that
[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md)
expects, using the oncogene panel bundled with EpiTracer. The returned
region is the gene span padded by `flank` – a plausible excised
neighbourhood – and carries the gene itself as the locus that must be
preserved through shattering.

## Usage

``` r
seed_locus(gene, genome = c("hg38", "hg19", "mm10"), flank = 250000)
```

## Arguments

- gene:

  Gene symbol(s), e.g. `"EGFR"` or `c("CDK4", "MDM2")`.

- genome:

  One of `"hg38"`, `"hg19"`, `"mm10"`.

- flank:

  Bp of neighbouring sequence co-excised with the gene (default `2.5e5`,
  giving a sub-megabase episome typical of a simple excision).

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
of the seed region(s), with a `preserved` attribute holding the gene
spans and a `gene` metadata column.

## See also

[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`gene_locus()`](https://alhafidzhamdan.github.io/EpiTracer/reference/gene_locus.md)

## Examples

``` r
seed_locus("EGFR")
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames            ranges strand |        gene
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]     chr7 54769021-55461628      * |        EGFR
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
attr(seed_locus("EGFR"), "preserved")
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames            ranges strand |        gene
#>          <Rle>         <IRanges>  <Rle> | <character>
#>   [1]     chr7 55019021-55211628      * |        EGFR
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
