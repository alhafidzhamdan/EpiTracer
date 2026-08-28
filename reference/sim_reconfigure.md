# Reconfigure an ecDNA while preserving its founder junction

Shatters a circle and religates it into a NEW circle, but keeps a
nominated junction — by default the oldest amplicon junction, i.e. the
founder circularisation duplication — intact as the closing junction.
Everything between its two ends is cut, reordered and reoriented at
random.

## Usage

``` r
sim_reconfigure(
  x,
  n_breaks = 10L,
  copies = NULL,
  preserve = "founder",
  del_p = 0.05,
  inv_p = 0.5,
  fitness = 1,
  label = "reconfigured",
  homology = c("nhej", "mmej", "nahr"),
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).

- n_breaks:

  Integer shattering breakpoints applied to the interior (default `10`).

- copies:

  Copies per cell of the reconfigured lineage (default `NULL`, meaning
  it takes the whole source species — a clean sweep). Give a number to
  have it arise as a subclone alongside its parent.

- preserve:

  Which junction to keep: `"founder"` (the oldest amplicon junction, the
  default), `"none"`, or a one-row `data.frame`/`data.table` with
  `start1` and `start2` naming the junction's two breakends.

- del_p, inv_p:

  Fragment loss and inversion probabilities for the interior (defaults
  `0.05` and `0.5`). Fragments carrying a preserved oncogene are never
  lost.

- fitness:

  Replication rate of the new lineage relative to its parent (default
  `1`); see
  [`sim_internal_deletion()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_internal_deletion.md).

- label:

  Name for the new species (default `"reconfigured"`).

- homology:

  Repair pathway sealing the new junctions (default `"nhej"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).
When the nominated junction's ends have been cut away by an earlier
round it cannot be preserved; the reconfiguration still happens and the
result carries the attribute `"preserved_founder" = FALSE`.

## Details

The preservation is **imposed, not derived**. This operator exists to
explore the counterfactual — what an amplicon would look like if one
junction were carried through repeated rebuilding — and the resulting
population shows a junction hierarchy ordered by age, with the preserved
junction clonal and each later round a step below, on a population in
which nothing is unshattered.

## There is no known selection for a specific junction

Do not read the preservation as a mechanism. A religated circle is a
circle whichever junction closes it: circularity is topological, so a
descendant that loses the founder duplication is no less viable than one
that keeps it, and nothing distinguishes an arbitrary duplication
between two arbitrary ends as a target of selection. Under free
reshuffling the founder survives a round with probability about
`(1 - del_p) / (2 * n_breaks)`, so its chance of remaining clonal across
several rounds is negligible.

A **clonal** founder duplication in real data is therefore evidence that
the molecules carrying it were never extensively reshuffled — the
parsimonious reading being that chromothripsis struck a minority of the
population and the rest was untouched. The one indirect route to genuine
selection is on the whole CONFIGURATION rather than the junction:
retaining the founder junction amounts to retaining the excised segment
contiguous and in reference orientation, which preserves the oncogene's
regulatory neighbourhood, and ecDNA is known to carry co-amplified
enhancers whose spacing to the promoter matters. But that argument
selects for the intact arrangement as a whole, and so would preserve the
internal junctions too — not the founder alone, which is what this
operator does. Treat `preserve` as a hypothesis to test, not as biology
already established.

## See also

[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_shatter_to_episomes()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_shatter_to_episomes.md),
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "RC01", copies = 90)
## three rounds of reconfiguration, each a subclone, nothing left unshattered
for (i in 1:3) ec <- sim_reconfigure(ec, n_breaks = 12, copies = 12, seed = i)
ec
#> <ecdna_sim> RC01 
#>   structure   : circular (episome) 
#>   fragments   : 1 on chr7 
#>   circle size : 2.59 Mb
#>   copies/cell : 90 
#>   population  : 4 circle species
#>       founder           54.0 copies (60.0%)  fitness 1.00    1 frag    2593 kb
#>       reconfigured      12.0 copies (13.3%)  fitness 1.00   12 frag    2114 kb
#>       reconfigured      12.0 copies (13.3%)  fitness 1.00   11 frag    2019 kb
#>       reconfigured      12.0 copies (13.3%)  fitness 1.00   13 frag    2593 kb
#>   junctions   : 35 (DEL=4, DUP=9, h2hINV=11, t2tINV=11) 
#>   rounds      : 3 
#>   history     : simple_excision -> reconfiguration x3 
```
