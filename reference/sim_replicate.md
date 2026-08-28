# Simulate rounds of ecDNA amplification

Amplifies a circle **without changing its structure**: the fragment list
and the junction ledger are untouched, and only `copies` (ecDNA
molecules per tumour cell) rises, by `fold` per round.

## Usage

``` r
sim_replicate(
  x,
  rounds = 1L,
  fold = 2,
  jitter = 0,
  max_cn = Inf,
  keep_all = FALSE,
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
  typically straight from
  [`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md).

- rounds:

  Integer number of replication rounds (default `1`).

- fold:

  Multiplicative gain in copies per round (default `2`, one doubling).
  Non-integer values model partial selective sweeps. When the amplicon
  is a mixture of circle species each grows by `fold * fitness`, so a
  species with a replication advantage expands its share of the
  population.

- jitter:

  Log-normal dispersion applied to each round's fold change (default
  `0`, i.e. deterministic). Set e.g. `0.15` to model uneven segregation
  between rounds.

- max_cn:

  Ceiling on peak copy number (default `Inf`; see
  [`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md)).

- keep_all:

  Logical; return a list of
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
  objects, one per round including the input as element 1 (default
  `FALSE`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
or a list of them when `keep_all = TRUE`.

## Details

`fold` is a **net growth factor of copy number under selection, not a
cell division**. Cell divisions do not raise ecDNA copy number: with no
centromere the copies partition at random between daughters, so
replication doubles and division halves, and the expected copy number
per cell is unchanged. The neutral value of `fold` is therefore `1`, and
`fold = 2` is an aggressive selective sweep. Use
[`sim_segregate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_segregate.md)
to let the level follow from a number of generations and a selection
coefficient instead of setting it by hand.

Because junction copy number is `multiplicity * copies`, and simulated
read support is proportional to junction copy number, the effect on the
emitted data is that the **boundary duplication gains read support round
after round** while no new junction appears. That is the signature of a
clean episome amplifying: one junction, rising `VF`, `PURPLE_JCN` and
`PURPLE_CN`, flanks unchanged at the host baseline.

Contrast
[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
which changes the structure (new junctions, oscillating copy number)
rather than only the level.

## See also

[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "REP01", copies = 4)
traj <- sim_replicate(ec, rounds = 5, fold = 2, keep_all = TRUE)

## structure is fixed; only the level rises
do.call(rbind, lapply(traj, summary))[, c("copies", "max_cn", "n_junctions")]
#>    copies max_cn n_junctions
#>     <num>  <num>       <int>
#> 1:      4      4           2
#> 2:      8      8           2
#> 3:     16     16           2
#> 4:     32     32           2
#> 5:     64     64           2
#> 6:    128    128           2
```
