# Score EpiTracer's mechanism callers against simulated ground truth

Runs the mechanism callers over a simulated cohort and joins their
per-amplicon verdicts to the known truth, returning both the per-sample
calls and a summary of how each trajectory class was classified. This is
how the simulation suite is meant to be consumed: build a cohort whose
mechanisms are known by construction, then ask what EpiTracer says about
it.

## Usage

``` r
sim_benchmark(
  cohort,
  callers = c("simple_excision", "chromothripsis", "micronucleation"),
  mc.cores = 1L
)
```

## Arguments

- cohort:

  The list returned by
  [`sim_cohort()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_cohort.md).

- callers:

  Character vector naming which callers to run; any of
  `"simple_excision"`, `"chromothripsis"`, `"micronucleation"`, `"brf"`,
  `"bfb"`. Defaults to the first three, the mechanisms this suite
  generates.

- mc.cores:

  Cores passed to the callers (default `1`).

## Value

A list with `calls` (one row per amplicon: truth columns plus one
logical per caller) and `summary` (per-class firing rate of each
caller).

## Details

Note the classes are **not** mutually exclusive by design — a chimeric
amplicon has genuinely been through micronucleation *and*
chromothripsis, so a `TRUE` from both callers is correct, not a false
positive. The summary therefore reports the rate at which each caller
fires per class rather than a single confusion matrix.

## See also

[`sim_cohort()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_cohort.md)

## Examples

``` r
co <- sim_cohort(n = 6, seed = 1)
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
#> Error in S4Vectors:::load_package_gracefully("XVector", "by the range() method ",     "for CompressedIRangesList objects"): Could not load package XVector. Is it installed?
#> 
#>   Note that the XVector package is required by the range() method for
#>   CompressedIRangesList objects. Please install it with:
#> 
#>     BiocManager::install("XVector")
bm <- sim_benchmark(co)
#> Error: object 'co' not found
bm$summary
#> Error: object 'bm' not found
```
