# Summarise a simulated amplicon

Ground-truth summary statistics of an
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
object – the quantities the mechanism callers try to recover, so a
simulation can be checked against what EpiTracer reports for it.

## Usage

``` r
# S3 method for class 'ecdna_sim'
summary(object, ...)
```

## Arguments

- object:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
  from
  [`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md)
  and friends.

- ...:

  Ignored.

## Value

A one-row
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with the sample, structure, fragment and junction counts, per-class
junction counts, circle size, copies per cell, peak copy number, number
of chromosomes involved and the mechanism history.

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01")
summary(ec)
#>    sample circular n_fragments n_chr circle_mb copies n_variants max_cn
#>    <char>   <lgcl>       <int> <int>     <num>  <num>      <int>  <num>
#> 1:  SIM01     TRUE           1     1  0.692608     40          1     40
#>    n_junctions n_del n_dup n_h2h n_t2t n_tra rounds         history
#>          <int> <int> <int> <int> <int> <int>  <int>          <char>
#> 1:           2     1     1     0     0     0      0 simple_excision
```
