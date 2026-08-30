# EpiTracer

<!-- badges: start -->
[![R-CMD-check](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/alhafidzhamdan/EpiTracer/actions/workflows/pkgdown.yaml/badge.svg)](https://alhafidzhamdan.github.io/EpiTracer/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

> **EpiTracer** — *calling and visualising extrachromosomal circular DNA
> (ecDNA) amplicons formed by simple excision, and telling them apart from the
> other mechanisms that drive focal amplification.*

EpiTracer works from whole-genome sequencing (WGS) copy-number and
structural-variant calls. It **nominates episomal ecDNA** — circles born from a
simple excision-and-recircularisation event — by finding the **founder
circularisation junction** that bounds the amplicon, and it provides a matching
set of **mechanism callers**, **rearrangement plotters**, and an **ecDNA
simulation suite**. Click any function for its documentation.

### Episomal ecDNA

| Function | Purpose |
| :--- | :--- |
| [**`call_simple_excision()`**](docs/call_simple_excision.md) | The episomal ecDNA caller (formerly `call_episomal_ecdna()`) — flags amplicons bounded by a copy-gaining **boundary duplication** with non-gained flanks and a deletion **"excision scar"**. |
| [**`call_founder_boundary()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_founder_boundary.html) | Nominates the **founder circularisation junction** via VF-stratified reconstruction — the founder must be a *copy-gaining* junction (a boundary DUP = simple excision; an inversion = inverted duplication, not episomal; a deletion is ancestral). |
| [**`detect_amplicon_seeds()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/detect_amplicon_seeds.html) | Detects focal-amplicon "seeds" from copy number alone, so the caller can run **without AmpliconArchitect**. |

### Mechanism callers — *is the amplification really episomal, or something else?*

Each takes the same WGS inputs and reports one mechanism's structural signature,
so a focal amplification can be assigned to its likely route of formation.

| Function | Mechanism |
| :--- | :--- |
| [**`call_chromothripsis()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.html) | Chromothripsis **within** an (episomal) amplicon — clustered breakpoints, random fragment joins, copy-number oscillation. |
| [**`call_micronucleation()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_micronucleation.html) | Micronucleation followed by chromothripsis. |
| [**`call_bfb()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.html) | Breakage–fusion–bridge (BFB). |
| [**`call_brf()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_brf.html) | Breakage–replication/fusion (BRF). |
| [**`call_translocation_bridge_amp()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_translocation_bridge_amp.html) | Translocation-bridge amplification (TBA). |
| [**`call_chromoplexy()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromoplexy.html) | Chromoplexy — closed, balanced rearrangement chains. |
| [**`call_wgd()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/call_wgd.html) | Whole-genome doubling (WGD) from allele-specific copy number. |

### Visualisation

| Function | Purpose |
| :--- | :--- |
| [**`plot_sv_linear()`**](docs/plot_sv_linear.md) | General-purpose linear plotter for structural variants and copy number, with an optional stacked SNV panel. |
| [**`plot_sv_reconstruction()`**](docs/plot_sv_reconstruction.md) | Reconstructs a focal amplicon wave by wave — stacks SV + CN panels ordered by read support (`VF`). |
| [**`plot_sv_circos()`**](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_circos.html) | Whole-genome structural-variant circos plot for one sample. |

### Simulation

An `sim_*` suite generates ecDNA amplicons with **known mechanisms** — a circle
at birth, after rounds of micronucleation and chromothripsis, or fused from two
episomes — and renders them as caller or plotter inputs, for benchmarking the
callers above. See the
[**simulation reference**](https://alhafidzhamdan.github.io/EpiTracer/reference/index.html#simulation).

📖 **Full documentation & function reference:** [**alhafidzhamdan.github.io/EpiTracer**](https://alhafidzhamdan.github.io/EpiTracer/) · get started with `vignette("epitracer")`.

## Vignettes

| Article | What it covers |
| :--- | :--- |
| [**Get started**](https://alhafidzhamdan.github.io/EpiTracer/articles/epitracer.html) | The end-to-end walk-through, from inputs to a first episomal call and figure. |
| [**Calling episomal ecDNA**](https://alhafidzhamdan.github.io/EpiTracer/articles/calling.html) | `call_simple_excision()` in depth — the episome heuristic, the excision scar, and running with or without AmpliconArchitect. |
| [**Calling amplicon-formation mechanisms**](https://alhafidzhamdan.github.io/EpiTracer/articles/mechanisms.html) | Runnable examples for `call_chromothripsis()`, `call_micronucleation()`, `call_bfb()`, `call_brf()` and `call_chromoplexy()` — each mechanism's structural signature. |
| [**Plotting copy number, structural variants and kataegis**](https://alhafidzhamdan.github.io/EpiTracer/articles/plotting.html) | Building `plot_sv_linear()` figures, including the stacked SNV/kataegis panel. |
| [**Reconstructing a focal amplicon wave by wave**](https://alhafidzhamdan.github.io/EpiTracer/articles/reconstruction.html) | Using `plot_sv_reconstruction()` to order rearrangements by read support and read off the founder junction. |

> [!TIP]
> **[Ongoing work] Interactive ecDNA simulations** accompany EpiTracer:
> **[alhafidzhamdan.github.io/simulations](https://alhafidzhamdan.github.io/simulations/)** —
> browser-based models of ecDNA birth and evolution over 100 cell divisions (neutral drift,
> positive selection with structural evolution, a micronucleus/chromothripsis origin, and
> spatial selection across tumour niches), with genome profiles drawn in the EpiTracer style.

> [!IMPORTANT]
> `call_simple_excision()` works best with focal amplicon calls from
> [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect)
> together with allele-specific copy-number (CNV) and structural-variant (SV)
> calls from
> [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple), part
> of the Hartwig Medical Foundation (HMF) pipeline.
> **AmpliconArchitect is optional**: pass `ecdna_gr = NULL` and EpiTracer
> detects focal-amplicon seeds directly from the copy-number segments
> (`detect_amplicon_seeds()`), so it can run from SV + CN alone. Accurate
> allele-specific copy-number segmentation remains important for reliably
> nominating amplicons formed via simple excision.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("alhafidzhamdan/EpiTracer")
```

## Quick start

A small example dataset ships with the package, so you can confirm your
installation works in seconds:

```r
library(EpiTracer)

# Bundled example: a textbook episomal EGFR ecDNA
res <- call_simple_excision(
  ecdna_gr        = ex_caller_inputs$ecdna_gr,
  breakpoints_gr  = ex_caller_inputs$breakpoints_gr,
  cnv_gr          = ex_caller_inputs$cnv_gr,
  cancer_genes_gr = ex_caller_inputs$cancer_genes_gr
)
unique(res$episomal)          #> "TRUE"

# ...or run it without AmpliconArchitect — detect seeds from copy number:
call_simple_excision(
  ecdna_gr        = NULL,
  breakpoints_gr  = ex_caller_inputs$breakpoints_gr,
  cnv_gr          = ex_caller_inputs$cnv_gr,
  cancer_genes_gr = ex_caller_inputs$cancer_genes_gr
)
```

See the **[get-started vignette](https://alhafidzhamdan.github.io/EpiTracer/articles/epitracer.html)**
for the full walk-through, including the `plot_sv_linear()` figure.

## Status

**v0.0.1 — alpha.** Under active development; the API may still change.

## License

Released under the [MIT License](LICENSE). © 2026 Alhafidz Hamdan.
