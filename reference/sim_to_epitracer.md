# Render a simulated amplicon as EpiTracer caller inputs

Turns an
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
into the four GRanges
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
and the mechanism callers take, applying a model of what a short-read
pipeline would actually report: finite purity and depth, copy-number
segmentation error and resolution, junction dropout and false positives,
and repair-pathway microhomology. With `noise = sim_noise(0)` the
rendering is exact.

## Usage

``` r
sim_to_epitracer(x, noise = sim_noise(), min_cn_ratio = 3, seed = NULL)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
  or a list of them (one per sample) to render as one cohort.

- noise:

  A
  [`sim_noise()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_noise.md)
  specification, or a numeric scale passed to it.

- min_cn_ratio:

  Amplification threshold used to define the reported amplicon
  footprint, matching
  [`detect_amplicon_seeds()`](https://alhafidzhamdan.github.io/EpiTracer/reference/detect_amplicon_seeds.md)
  (default `3`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A named list of four
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
— `ecdna_gr`, `breakpoints_gr`, `cnv_gr`, `cancer_genes_gr` — plus
`truth`, a
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of the ground-truth mechanism labels per sample. The first four can be
passed straight to any EpiTracer caller.

## See also

[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`sim_to_plot_inputs()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_plot_inputs.md),
[`sim_cohort()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_cohort.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01", copies = 50)
inp <- sim_to_epitracer(ec, seed = 1)
#> Error in S4Vectors:::load_package_gracefully("XVector", "by the range() method ",     "for CompressedIRangesList objects"): Could not load package XVector. Is it installed?
#> 
#>   Note that the XVector package is required by the range() method for
#>   CompressedIRangesList objects. Please install it with:
#> 
#>     BiocManager::install("XVector")
res <- call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                            inp$cnv_gr, inp$cancer_genes_gr)
#> Error: object 'inp' not found
unique(res$episomal)
#> Error: object 'res' not found
```
