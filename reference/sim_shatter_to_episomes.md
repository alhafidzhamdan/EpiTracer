# Simulate shattering into several new episomes

The ecDNA is encapsulated in a micronucleus and shattered, and the
surviving fragments religate into **several separate circles** rather
than one. Each daughter is a new episome with its own circularisation
junction, and the population that results is a mixture of them.

## Usage

``` r
sim_shatter_to_episomes(
  x,
  n_breaks = 20L,
  n_circles = 3L,
  del_p = 0.1,
  inv_p = 0.5,
  require_oncogene = TRUE,
  keep_parent = FALSE,
  copies = 1,
  fitness = 1,
  homology = c("nhej", "mmej", "nahr"),
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).
  The most abundant circle species is the one captured in the
  micronucleus.

- n_breaks:

  Integer number of shattering breakpoints (default `20`).

- n_circles:

  Number of daughter circles the fragments religate into (default `3`).

- del_p:

  Probability a fragment is lost outright during religation (default
  `0.1`).

- inv_p:

  Probability a fragment is inverted (default `0.5`).

- require_oncogene:

  Logical; drop daughters that carry no preserved locus (default
  `TRUE`). A circle without the oncogene confers no advantage and is
  diluted out. Set `FALSE` to keep every daughter, modelling the moment
  before selection has acted.

- keep_parent:

  Logical; does the intact parent circle survive alongside the daughters
  (default `FALSE`)? `FALSE` models the whole ecDNA population passing
  through the micronucleus, so the daughters inherit and split its
  copies. `TRUE` models a single molecule being shattered while the rest
  of the population carries on, giving a parent that keeps its founder
  junction plus daughters arising as minority species.

- copies:

  Copies per cell of each daughter when `keep_parent = TRUE` (default
  `1`, a single molecule). Ignored when `keep_parent = FALSE`, where the
  parent's copies are split evenly across the surviving daughters.

- fitness:

  Replication rate of the daughters relative to the parent (default
  `1`). Daughters are shorter, so a value above `1` is the expected
  direction; see
  [`sim_internal_deletion()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_internal_deletion.md).

- homology:

  Repair pathway sealing the new junctions (default `"nhej"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
whose population is the surviving daughters (plus the parent when
`keep_parent = TRUE`). Its `history` gains `"micronuclear_fission"`.

## Details

This is the branch
[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md)
does not cover: that operator religates the fragments back into a single
circle, whereas here the fragments partition, so one amplicon becomes a
family of smaller ones.

Two consequences fall out that are worth simulating for:

- **The founder boundary duplication is usually destroyed.** It survives
  only if some daughter happens to rejoin the two ends it joined, which
  is unlikely once the circle has been cut in many places. An amplicon
  that has been through fission therefore keeps a high, flat copy number
  over whatever its surviving daughters retain, but has no
  duplication-orientation junction spanning it — so
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
  no longer recognises it, even though it is still, genuinely, a set of
  episomes.

- **Copy number becomes patchy.** Sequence that ends up in a daughter
  that was lost, or in one that failed selection, drops back to the host
  baseline, so the flat amplicon of the parent breaks into blocks at
  different levels.

## References

Shoshani, O. *et al.* Chromothripsis drives the evolution of gene
amplification in cancer. *Nature* **591**, 137-141 (2021).

## See also

[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "FIS01", copies = 60)
fis <- sim_shatter_to_episomes(ec, n_breaks = 20, n_circles = 3, seed = 1)
fis                                   # several daughter species
#> <ecdna_sim> FIS01 
#>   structure   : circular (episome) 
#>   fragments   : 9 on chr7 
#>   circle size : 0.25 Mb
#>   copies/cell : 60 
#>   population  : 2 circle species
#>       episome1          30.0 copies (50.0%)  fitness 1.00    9 frag     246 kb
#>       episome2          30.0 copies (50.0%)  fitness 1.00    6 frag     250 kb
#>   junctions   : 16 (DEL=4, DUP=4, h2hINV=4, t2tINV=4) 
#>   rounds      : 1 
#>   history     : simple_excision -> micronuclear_fission 
## the founder boundary duplication is gone
fis$junctions[origin == "amplicon" & mechanism == "simple_excision"]
#> Empty data.table (0 rows and 12 cols): chrom1,start1,strand1,chrom2,start2,strand2...
```
