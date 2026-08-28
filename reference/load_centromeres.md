# Centromere regions for a bundled genome

Returns the centromere (cytoband `acen`) span per chromosome from the
bundled cytoband data. Used by
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
to flag and exclude "chromosomal bridge" events – amplicons whose
amplified edge terminates *within* a centromere (a dicentric / bridge
signature) rather than forming a self-contained episome.

## Usage

``` r
load_centromeres(genome = c("hg38", "hg19", "mm10"))
```

## Arguments

- genome:

  One of `"hg38"`, `"hg19"`, `"mm10"`.

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
with one range per chromosome (the centromere span). Empty if the
genome's cytobands carry no `acen` bands.

## Examples

``` r
load_centromeres("hg38")
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> GRanges object with 24 ranges and 0 metadata columns:
#>        seqnames              ranges strand
#>           <Rle>           <IRanges>  <Rle>
#>    [1]     chr1 121700000-125100000      *
#>    [2]     chr2   91800000-96000000      *
#>    [3]     chr3   87800000-94000000      *
#>    [4]     chr4   48200000-51800000      *
#>    [5]     chr5   46100000-51400000      *
#>    ...      ...                 ...    ...
#>   [20]    chr20   25700000-30400000      *
#>   [21]    chr21   10900000-13000000      *
#>   [22]    chr22   13700000-17400000      *
#>   [23]     chrX   58100000-63800000      *
#>   [24]     chrY   10300000-10600000      *
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```
