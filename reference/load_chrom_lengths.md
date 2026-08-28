# Chromosome lengths for a bundled genome

Returns the length (bp) of each chromosome from the bundled cytoband
data (the maximum band coordinate per chromosome). Used by
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
to test telomere proximity for the breakage-fusion-bridge (BFB)
annotation: a p-telomere sits at coordinate ~0 and a q-telomere at ~the
chromosome length.

## Usage

``` r
load_chrom_lengths(genome = c("hg38", "hg19", "mm10"))
```

## Arguments

- genome:

  One of `"hg38"`, `"hg19"`, `"mm10"`.

## Value

A named numeric vector of chromosome lengths (names are `seqnames`, e.g.
`"chr7"`). Empty if the genome's cytoband data is unavailable.

## Examples

``` r
load_chrom_lengths("hg38")["chr7"]
#>      chr7 
#> 159345973 
```
