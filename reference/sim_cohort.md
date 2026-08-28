# Simulate a cohort of ecDNA amplicons with known mechanisms

Draws `n` independent samples from a mixture of ecDNA trajectories and
renders them as one set of EpiTracer inputs with a ground-truth label
per sample. The three trajectory classes correspond to the questions the
suite was built to answer: what a **freshly born** episome looks like,
what it looks like after **one or several rounds** of micronucleation
and chromothripsis, and what a **chimera** of two co-encapsulated
episomes looks like.

## Usage

``` r
sim_cohort(
  n = 30,
  mixture = c(episomal = 1, chromothriptic = 1, chimeric = 1),
  genes = .sim_default_genes,
  rounds = 1:3,
  n_breaks = 8:30,
  copies = c(15, 80),
  host_ploidy = c(2, 2, 2, 3),
  homology = c("nhej", "nhej", "mmej"),
  noise = sim_noise(),
  prefix = "SIM",
  seed = NULL,
  mc.cores = 1L
)
```

## Arguments

- n:

  Number of samples (default `30`).

- mixture:

  Named numeric vector of class weights, normalised internally.
  Recognised names are `episomal` (born and unshattered),
  `chromothriptic` (one or more micronucleation rounds) and `chimeric`
  (two episomes fused). Defaults to an even split.

- genes:

  Character vector of oncogene symbols to draw seed loci from (default:
  a panel of recurrently ecDNA-amplified oncogenes).

- rounds:

  Integer vector of micronucleation rounds to sample from for the
  `chromothriptic` class (default `1:3`).

- n_breaks:

  Integer vector of chromothriptic breakpoint counts to sample from
  (default `8:30`).

- copies:

  Numeric length-2 range of ecDNA copies per cell at birth (default
  `c(15, 80)`).

- host_ploidy:

  Numeric vector of host baseline ploidies to sample from (default
  `c(2, 2, 2, 3)`, so a quarter of samples sit on a polysomic chromosome
  – the case that defeats a ploidy-based flank test).

- homology:

  Repair pathways to sample the circularisation junction from (default
  `c("nhej", "nhej", "mmej")`).

- noise:

  A
  [`sim_noise()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_noise.md)
  specification, or a numeric scale.

- prefix:

  Sample-name prefix (default `"SIM"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- mc.cores:

  Cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`).

## Value

A named list: the four caller inputs (`ecdna_gr`, `breakpoints_gr`,
`cnv_gr`, `cancer_genes_gr`), a `truth`
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with one row per simulated amplicon, and `sims`, the list of underlying
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
objects (so any sample can be re-rendered at different noise, or
plotted).

## See also

[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`sim_evolve()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_evolve.md),
[`sim_fuse_episomes()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_fuse_episomes.md),
[`sim_benchmark()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_benchmark.md)

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
co$truth[, .(WGS_ID, class, rounds, chimeric, truth_max_cn)]
#>     WGS_ID          class rounds chimeric truth_max_cn
#>     <char>         <char>  <int>   <lgcl>        <num>
#> 1: SIM0001 chromothriptic      3    FALSE    200.00000
#> 2: SIM0002       chimeric      1     TRUE    146.60006
#> 3: SIM0003       chimeric      1     TRUE     77.72181
#> 4: SIM0004       episomal      0    FALSE     15.51498
#> 5: SIM0005 chromothriptic      2    FALSE    200.00000
#> 6: SIM0006       episomal      0    FALSE     67.70956
```
