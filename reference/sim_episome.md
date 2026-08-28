# Simulate the birth of an episomal ecDNA

Creates a simple-excision episome: one or more genomic regions are
excised from the chromosome and circularised. This is the mechanism
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
is built to detect, and the simulated amplicon carries its three
signatures by construction — a **duplication-orientation junction at the
amplicon boundary** (the circularisation junction, `-/+`, emitted as
`svclass == "DUP"`), **non-gained flanks** (the host chromosome stays at
its baseline ploidy immediately outside the amplicon), and an **excision
scar** (a low-copy deletion junction joining the flanks on the
chromosome the circle left behind).

## Usage

``` r
sim_episome(
  regions,
  sample = "SIM01",
  copies = 40,
  preserved = NULL,
  host_ploidy = 2,
  host_flank = 1.5e+07,
  scar = TRUE,
  retain_source = TRUE,
  homology = c("nhej", "mmej", "nahr"),
  genome = "hg38",
  seed = NULL
)
```

## Arguments

- regions:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of the region(s) to excise. Several regions produce a circle assembled
  from several segments, as when an episome forms from a locus already
  broken by internal deletions. Use
  [`seed_locus()`](https://alhafidzhamdan.github.io/EpiTracer/reference/seed_locus.md)
  to build one from a gene symbol.

- sample:

  Character sample identifier carried into the emitted inputs.

- copies:

  Numeric ecDNA copies per tumour cell (default `40`). Copy number is
  `copies * coverage(fragments)` over the host baseline, so this sets
  the amplification level; patient ecDNA amplicons typically sit at
  20-100.

- preserved:

  Oncogene loci that later shattering must not delete. Defaults to the
  `preserved` attribute of `regions` when
  [`seed_locus()`](https://alhafidzhamdan.github.io/EpiTracer/reference/seed_locus.md)
  made it, else to `regions` themselves.

- host_ploidy:

  Numeric baseline copy number of the host chromosome (default `2`). Set
  higher to model an episome on a polysomic chromosome – e.g. `3` for
  the chr7 gain that accompanies EGFR ecDNA in glioblastoma.

- host_flank:

  Bp of host chromosome to represent on each side of the amplicon
  (default `1.5e7`), giving the flanks the callers test.

- scar:

  Logical; emit the chromosomal excision scar (default `TRUE`). Set
  `FALSE` to model an episome whose source chromosome was lost, so the
  scar is unobservable.

- retain_source:

  Logical; does the source chromosome keep the excised region on its
  other allele (default `TRUE`)? When `TRUE` the excised interval sits
  one copy below the flanks on the chromosome, the usual single-allele
  excision.

- homology:

  Repair pathway that sealed the circularisation junction, one of
  `"nhej"` (blunt), `"mmej"` (2-10 bp microhomology) or `"nahr"` (\>=14
  bp homology). Sets `HOMLEN` on the emitted breakends, the
  breakpoint-homology signature
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
  reports.

- genome:

  Genome build for chromosome lengths (default `"hg38"`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for
  reproducibility.

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
object.

## References

Shoshani, O. *et al.* Chromothripsis drives the evolution of gene
amplification in cancer. *Nature* **591**, 137-141 (2021).

## See also

[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01", copies = 50)
ec
#> <ecdna_sim> SIM01 
#>   structure   : circular (episome) 
#>   fragments   : 1 on chr7 
#>   circle size : 0.69 Mb
#>   copies/cell : 50 
#>   junctions   : 2 (DEL=1, DUP=1) 
#>   rounds      : 0 
#>   history     : simple_excision 
summary(ec)
#>    sample circular n_fragments n_chr circle_mb copies n_variants max_cn
#>    <char>   <lgcl>       <int> <int>     <num>  <num>      <int>  <num>
#> 1:  SIM01     TRUE           1     1  0.692608     50          1     50
#>    n_junctions n_del n_dup n_h2h n_t2t n_tra rounds         history
#>          <int> <int> <int> <int> <int> <int>  <int>          <char>
#> 1:           2     1     1     0     0     0      0 simple_excision
```
