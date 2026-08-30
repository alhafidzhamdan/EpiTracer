# EpiTracer

<!-- badges: start -->
[![R-CMD-check](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/pkgdown.yaml/badge.svg)](https://alhafidzhamdan.github.io/EpiTracer/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

> **EpiTracer** — calling and visualising episomal extrachromosomal DNA (ecDNA)
> amplicons formed by simple excision, and telling them apart from the other
> mechanisms of focal amplification, from WGS copy-number and structural-variant
> calls.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("alhafidzhamdan/EpiTracer")
```

## Quick start

```r
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

## Documentation

Everything — input formats, the founder-boundary and mechanism callers, and the
plotters — is on the **[package website](https://alhafidzhamdan.github.io/EpiTracer/)**
and its vignettes:

- [Get started](https://alhafidzhamdan.github.io/EpiTracer/articles/epitracer.html)
- [Calling episomal ecDNA](https://alhafidzhamdan.github.io/EpiTracer/articles/calling.html)
- [Plotting copy number, structural variants and kataegis](https://alhafidzhamdan.github.io/EpiTracer/articles/plotting.html)
- [Reconstructing a focal amplicon wave by wave](https://alhafidzhamdan.github.io/EpiTracer/articles/reconstruction.html)

## Status

**v0.0.1 — alpha.** Under active development; the API may still change.
Released under the [MIT License](LICENSE). © 2026 Alhafidz Hamdan.
