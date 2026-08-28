# Simulate one round of micronucleation and chromothripsis

Encapsulates the ecDNA in a micronucleus, shatters it, and religates the
surviving fragments in random order and orientation back into a circle —
the ecDNA -\> micronucleus -\> chromothripsis route by which a clean
episome becomes a complex amplicon.

## Usage

``` r
sim_micronucleation(
  x,
  n_breaks = 10L,
  del_p = 0.2,
  inv_p = 0.5,
  dup_p = 0.25,
  max_dup = 3L,
  amplify = 1.5,
  max_cn = 200,
  homology = c("nhej", "mmej", "nahr"),
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md),
  from
  [`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md)
  or a previous round.

- n_breaks:

  Integer number of chromothriptic breakpoints (default `10`). Patient
  chromothripsis footprints typically carry 10-50.

- del_p:

  Probability a non-preserved fragment is lost (default `0.2`).

- inv_p:

  Probability a surviving fragment is inverted (default `0.5`).

- dup_p:

  Per-trial probability of an extra incorporation of a fragment (default
  `0.25`), giving the copy-number oscillation described above. Set `0`
  for strict shatter-and-rejoin with no reduplication.

- max_dup:

  Maximum times one fragment may appear in the religated circle (default
  `3`).

- amplify:

  Multiplicative change in copies per cell across the round (default
  `1.5`), modelling selection for the fitter, oncogene-retaining circle.
  Set `1` for no amplification.

- max_cn:

  Ceiling on the amplicon's peak copy number (default `200`). Copies per
  cell are scaled back so `copies * max(coverage)` stays at or below it,
  modelling the plateau beyond which further oncogene dosage buys no
  selective advantage (the `N_fitness` saturation of the Bernhard *et
  al.* model). Without it, compounding `amplify` and fragment
  reduplication drive copy number far past anything observed in
  patients.

- homology:

  Repair pathway sealing the new junctions (default `"nhej"`;
  micronuclear religation is canonically blunt-ended NHEJ).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
with the religated fragment list, the accumulated junction ledger and
`round` incremented.

## Details

Three things happen, in order:

1.  **Shattering.** The circle is cut at `n_breaks` positions drawn
    uniformly over its sequence (so wide fragments are cut more often
    than narrow ones).

2.  **Fragment loss and reduplication.** Each fragment is lost with
    probability `del_p`, except fragments carrying a preserved oncogene,
    which always survive — a circle that loses its oncogene is not
    selected. Each surviving fragment is incorporated
    `1 + rbinom(1, max_dup - 1, dup_p)` times. This reduplication is
    what generates the **oscillating copy-number profile** that
    distinguishes a chromothriptic amplicon from a clean episome (which
    would otherwise be a flat two-level profile), and it is the step
    [`call_chromothripsis()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.md)
    scores as `cn_oscillations`.

3.  **Random rejoining.** Fragments are shuffled and each is inverted
    with probability `inv_p`, then religated into a circle. Because
    order and orientation are random, the four intrachromosomal junction
    orientations (`DEL`, `DUP`, `h2hINV`, `t2tINV`) come out close to
    equally represented — the **random fragment joins** hallmark, which
    [`call_chromothripsis()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.md)
    tests with a chi-squared goodness-of-fit and which separates
    chromothripsis from orientation-biased mechanisms such as BFB.

Applying this repeatedly (see
[`sim_evolve()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_evolve.md))
compounds the rearrangement: junction counts grow roughly linearly in
the number of rounds while the circle contracts, reproducing the small,
junction-dense amplicons seen in patients.

## Circle populations

A shattering event happens inside ONE micronucleus, to ONE molecule, and
the surviving product is what sweeps. So when the amplicon is a mixture
of circle species (see
[`sim_internal_deletion()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_internal_deletion.md)),
this operator acts on the most abundant species and the population
collapses back to one — the shattered descendant. Minority species
present before the event are not carried through.

## References

Zhang, C.-Z. *et al.* Chromothripsis from DNA damage in micronuclei.
*Nature* **522**, 179-184 (2015).

## See also

[`sim_evolve()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_evolve.md),
[`sim_fuse_episomes()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_fuse_episomes.md),
[`call_chromothripsis()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01")
mn <- sim_micronucleation(ec, n_breaks = 12, seed = 1)
summary(mn)
#>    sample circular n_fragments n_chr circle_mb copies n_variants max_cn
#>    <char>   <lgcl>       <int> <int>     <num>  <num>      <int>  <num>
#> 1:  SIM01     TRUE          13     1  0.816901     60          1    120
#>    n_junctions n_del n_dup n_h2h n_t2t n_tra rounds
#>          <int> <int> <int> <int> <int> <int>  <int>
#> 1:          13     5     4     2     2     0      1
#>                               history
#>                                <char>
#> 1: simple_excision -> micronucleation
```
