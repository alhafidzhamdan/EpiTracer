# Simulate chimeric ecDNA formation by co-encapsulation of two episomes

Two independently born episomes are captured in the **same**
micronucleus, shattered together, and religated into a single chimeric
circle carrying fragments of both — the route by which an ecDNA acquires
sequence from a second, non-homologous locus.

## Usage

``` r
sim_fuse_episomes(
  x,
  y,
  n_breaks = 10L,
  del_p = 0.15,
  inv_p = 0.5,
  dup_p = 0.2,
  max_dup = 3L,
  copies = NULL,
  sample = NULL,
  homology = c("nhej", "mmej", "nahr"),
  seed = NULL
)
```

## Arguments

- x, y:

  Two
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
  objects to co-encapsulate. `x` supplies the sample name, host
  background and copy level unless overridden.

- n_breaks:

  Integer breakpoints applied to the pooled fragment set (default `10`).

- del_p, inv_p, dup_p, max_dup:

  Fragment loss, inversion and reduplication, as in
  [`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md).

- copies:

  Copies per cell of the chimera; defaults to the mean of the two
  parents, the chimera replacing both.

- sample:

  Sample name for the chimera (defaults to `x`'s).

- homology:

  Repair pathway sealing the new junctions (default `"nhej"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
whose fragments carry an `origin` column naming the parental episome
each came from, so the chimeric contribution is traceable.

## Details

When the two episomes come from different chromosomes the chimera
necessarily carries **inter-chromosomal junctions with both breakends
inside amplified copy number**, which is precisely the signature
[`call_chimeric_amplicon()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chimeric_amplicon.md)
detects. It also puts two oncogenes on one circle, the co-amplification
(e.g. *EGFR* with *CDK4*, or *MYC* with a distant enhancer) that is
common in patient amplicons and cannot arise from a single excision.

## Circle populations

A shattering event happens inside ONE micronucleus, to ONE molecule, and
the surviving product is what sweeps. So when the amplicon is a mixture
of circle species (see
[`sim_internal_deletion()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_internal_deletion.md)),
this operator acts on the most abundant species and the population
collapses back to one — the shattered descendant. Minority species
present before the event are not carried through.

## References

Hung, K. L. *et al.* ecDNA hubs drive cooperative intermolecular
oncogene expression. *Nature* **600**, 731-736 (2021).

## See also

[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`call_chimeric_amplicon()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chimeric_amplicon.md)

## Examples

``` r
a <- sim_episome(seed_locus("EGFR"), sample = "SIM02", copies = 40)
b <- sim_episome(seed_locus("CDK4"), sample = "SIM02", copies = 30)
ch <- sim_fuse_episomes(a, b, seed = 1)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
summary(ch)   # n_chr == 2, n_tra > 0
#> Warning: The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
#>    sample circular n_fragments n_chr circle_mb copies n_variants max_cn
#>    <char>   <lgcl>       <int> <int>     <num>  <num>      <int>  <num>
#> 1:  SIM02     TRUE          15     2   1.58216     35          1     70
#>    n_junctions n_del n_dup n_h2h n_t2t n_tra rounds
#>          <int> <int> <int> <int> <int> <int>  <int>
#> 1:          16     3     0     3     0    10      1
#>                                                        history
#>                                                         <char>
#> 1: simple_excision -> fuse(simple_excision) -> chimeric_fusion
```
