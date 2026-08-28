# Evolve a simulated ecDNA through several rounds

Applies
[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md)
`rounds` times, so one call produces the whole trajectory from a clean
episome to a heavily rearranged amplicon. Returns the final amplicon by
default, or the full trajectory (`keep_all = TRUE`) so the rounds can be
compared — the natural way to ask at what point EpiTracer stops calling
an amplicon episomal and starts calling it chromothriptic.

## Usage

``` r
sim_evolve(x, rounds = 1L, keep_all = FALSE, ..., seed = NULL)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).

- rounds:

  Integer number of micronucleation + chromothripsis rounds.

- keep_all:

  Logical; return a list of
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
  objects, one per round including the input as element 1 (default
  `FALSE`, returning only the final amplicon).

- ...:

  Passed to
  [`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md)
  (`n_breaks`, `del_p`, `inv_p`, `dup_p`, `amplify`, ...). A
  vector-valued `n_breaks` is recycled across rounds, so
  `n_breaks = c(20, 5)` models a severe first shattering followed by a
  milder one.

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) once, before the
  first round.

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
or a list of them when `keep_all = TRUE`.

## See also

[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SIM03")
traj <- sim_evolve(ec, rounds = 3, n_breaks = 12, keep_all = TRUE, seed = 1)
do.call(rbind, lapply(traj, summary))[, c("rounds", "n_junctions", "circle_mb")]
#>    rounds n_junctions circle_mb
#>     <int>       <int>     <num>
#> 1:      0           2  0.692608
#> 2:      1          13  0.816901
#> 3:      2          35  0.954532
#> 4:      3          63  1.366450
```
