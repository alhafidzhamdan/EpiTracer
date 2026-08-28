# Simulate an internal deletion arising within an ecDNA

A single ecDNA molecule (or `copies` of them) loses an interval from
inside the circle and re-seals. The altered molecule is drawn **out of**
the existing pool, so the cell now carries a mixture: intact circles
plus deleted ones, each replicating from then on.

## Usage

``` r
sim_internal_deletion(
  x,
  start = NULL,
  end = NULL,
  width = 50000,
  copies = 1,
  fitness = 1,
  label = "internal_del",
  homology = c("nhej", "mmej", "nahr"),
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).

- start, end:

  Explicit deletion boundaries. When `NULL` (the default) an interval of
  `width` bp is drawn at random from inside the circle, avoiding the
  preserved loci.

- width:

  Deletion size in bp when drawing at random (default `5e4`).

- copies:

  Number of ecDNA molecules that acquire the deletion (default `1` — the
  minimal event, a single molecule). These are taken out of the intact
  pool, so the total copy number is unchanged at the moment it happens.

- fitness:

  Replication rate of the deleted circle relative to the intact one
  (default `1`, neutral). This matters more than it looks: a NEUTRAL
  minority species keeps a **fixed share** of the population forever,
  because both species double together — one molecule in 32 stays at 3%
  however long it replicates, and its junction is emitted ~32-fold below
  the boundary. Real amplicons show internal junctions at roughly a half
  to a twentieth of the founder's read support, which needs the deleted
  circle to expand. A value above `1` does that, and is the biologically
  expected direction: a circle that has shed non-essential sequence but
  kept its oncogene is shorter, replicates faster and packages more
  copies per cell.

- label:

  Name for the new species in the population (default `"internal_del"`).

- homology:

  Repair pathway sealing the deletion (default `"nhej"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
carrying two circle species.
[`summary()`](https://rdrr.io/r/base/summary.html) reports the
population, and
[`sim_replicate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_replicate.md)
grows both.

## Details

This is what makes the *order* of events observable. The deleted circle
still carries the original boundary duplication, so that junction's copy
number is the whole population's; the new internal deletion is carried
only by the altered species, so its copy number — and therefore its read
support — is just that species' abundance. An internal deletion that
arose late sits far below the boundary junction it lies inside, which is
exactly the read-support stratification
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)
separates into waves. A deletion introduced at birth and allowed to
sweep would instead sit level with the boundary junction, and the two
events would be indistinguishable.

The deleted interval never overlaps a preserved oncogene: a circle that
loses its oncogene confers no advantage and would not persist.

## See also

[`sim_replicate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_replicate.md),
[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "DEL01", copies = 2)
ec <- sim_replicate(ec, rounds = 4)          # 32 intact circles
ec <- sim_internal_deletion(ec, width = 5e4, fitness = 1.3, seed = 1)
ec <- sim_replicate(ec, rounds = 4)          # the shorter circle expands
ec                                           # two species, DEL share rising
#> <ecdna_sim> DEL01 
#>   structure   : circular (episome) 
#>   fragments   : 1 on chr7 
#>   circle size : 0.69 Mb
#>   copies/cell : 541.7 
#>   population  : 2 circle species
#>       founder          496.0 copies (91.6%)  fitness 1.00    1 frag     693 kb
#>       internal_del      45.7 copies ( 8.4%)  fitness 1.30    2 frag     643 kb
#>   junctions   : 3 (DEL=2, DUP=1) 
#>   rounds      : 0 
#>   history     : simple_excision -> replication x4 -> internal_deletion -> replication x4 
```
