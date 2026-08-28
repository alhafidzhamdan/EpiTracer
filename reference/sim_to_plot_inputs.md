# Render a simulated amplicon as plotter inputs

The same simulated amplicon in the data-frame format
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)
and
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)
expect, so a simulated trajectory can be inspected with the package's
own viewers — the quickest way to see whether a parameter set reproduces
the look of a patient amplicon.

## Usage

``` r
sim_to_plot_inputs(x, noise = sim_noise(), seed = NULL)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
  or a list of them.

- noise:

  A
  [`sim_noise()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_noise.md)
  specification, or a numeric scale.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A named list of three `data.frame`s — `cnv_data`, `sv_data`, `wgd_data`
— matching
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)'s
arguments.

## See also

[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md),
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md),
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)

## Examples

``` r
ec <- sim_evolve(sim_episome(seed_locus("EGFR"), sample = "SIM04"),
                 rounds = 2, n_breaks = 15, seed = 1)
pin <- sim_to_plot_inputs(ec, seed = 1)
#> Error in S4Vectors:::load_package_gracefully("XVector", "by the range() method ",     "for CompressedIRangesList objects"): Could not load package XVector. Is it installed?
#> 
#>   Note that the XVector package is required by the range() method for
#>   CompressedIRangesList objects. Please install it with:
#> 
#>     BiocManager::install("XVector")
p <- plot_sv_linear(sample = "SIM04", cnv_data = pin$cnv_data,
                    sv_data = pin$sv_data, wgd_data = pin$wgd_data,
                    chromosome = "chr7")
#> Error: object 'pin' not found
```
