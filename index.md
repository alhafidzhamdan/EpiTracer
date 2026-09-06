# EpiTracer

> **EpiTracer** — calling and visualising episomal extrachromosomal DNA
> (ecDNA) amplicons formed by simple excision, and telling them apart
> from the other mechanisms of focal amplification, from WGS copy-number
> and structural-variant calls.

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("alhafidzhamdan/EpiTracer")
```

## Quick start

``` r

library(EpiTracer)

# bundled example: a textbook episomal EGFR ecDNA
res <- call_simple_excision(
  ecdna_gr        = ex_caller_inputs$ecdna_gr,      # or NULL to detect seeds from CN
  breakpoints_gr  = ex_caller_inputs$breakpoints_gr,
  cnv_gr          = ex_caller_inputs$cnv_gr,
  cancer_genes_gr = ex_caller_inputs$cancer_genes_gr
)
unique(res$episomal)          #> "TRUE"
```

Already have [AmpliconSuite](https://github.com/AmpliconSuite) output
(e.g. from [AmpliconRepository](https://ampliconrepository.org))? Apply
the same episome model directly to an AmpliconArchitect graph +
AmpliconClassifier interval BED:

``` r

call_episomal_from_ampliconsuite(
  graph_file = "SAMPLE_amplicon1_graph.txt",
  intervals  = "SAMPLE_amplicon1_ecDNA_1_intervals.bed"
)
```

## Documentation

Everything — input formats, the founder-boundary and mechanism callers,
and the plotters — is on the **[package
website](https://alhafidzhamdan.github.io/EpiTracer/)** and its
vignettes:

- [Get
  started](https://alhafidzhamdan.github.io/EpiTracer/articles/epitracer.html)
- [Calling episomal
  ecDNA](https://alhafidzhamdan.github.io/EpiTracer/articles/calling.html)
- [Calling episomal ecDNA from
  AmpliconSuite](https://alhafidzhamdan.github.io/EpiTracer/articles/ampliconsuite.html)
- [Plotting copy number, structural variants and
  kataegis](https://alhafidzhamdan.github.io/EpiTracer/articles/plotting.html)
- [Reconstructing a focal amplicon wave by
  wave](https://alhafidzhamdan.github.io/EpiTracer/articles/reconstruction.html)

## Status

**v0.1.0 — alpha.** Under active development; the API may still change.
Released under the [MIT
License](https://alhafidzhamdan.github.io/EpiTracer/LICENSE). © 2026
Alhafidz Hamdan.
