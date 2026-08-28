# Grow an ecDNA population by segregation and selection

Evolves the amplicon's copy number over `generations` of cell division
under an explicit model of how ecDNA actually propagates: each S phase
every circle replicates, and at mitosis the `2k` copies partition **at
random** between the two daughters (`Binom(2k, 1/2)`), because ecDNA has
no centromere. Cells are then sampled in proportion to a fitness that
rises with oncogene dosage.

## Usage

``` r
sim_segregate(
  x,
  generations = 20L,
  selection = 0.05,
  plateau = 40,
  n_cells = 2000L,
  addicted = TRUE,
  max_cn = Inf,
  seed = NULL
)
```

## Arguments

- x:

  An
  [ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md).

- generations:

  Integer number of cell divisions to run (default `20`).

- selection:

  Dosage-selection coefficient (default `0.05`). `0` is neutral; the
  population mean copy number then does not rise, however many
  generations elapse.

- plateau:

  Copy number above which extra dosage buys no further advantage
  (default `40`), the saturation seen in ecDNA fitness models.

- n_cells:

  Number of cells simulated (default `2000`). Larger is smoother and
  slower; the returned copy number is a population mean, so a few
  thousand cells is ample.

- addicted:

  Logical; do cells that lose the amplicon entirely die (default `TRUE`,
  oncogene addiction)?

- max_cn:

  Ceiling on peak copy number (default `Inf`).

- seed:

  Optional integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

An
[ecdna_sim](https://alhafidzhamdan.github.io/EpiTracer/reference/ecdna_sim.md)
with `copies` set to the evolved population mean. The per-cell
copy-number distribution is attached as the `"cells"` attribute, for
single-cell or heterogeneity work; bulk rendering ignores it. Species
proportions are rescaled together (co-segregation of distinct circle
species within a cell is not modelled).

## Details

Use it in place of
[`sim_replicate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_replicate.md)
when you want the copy level to *follow from* how long the tumour has
been evolving and how strongly the amplicon is selected, rather than
being set by hand.

## Why division rate does not raise read support

Random segregation is copy-number-neutral in expectation — replication
doubles, division halves, so `E[copies in daughter] = copies in parent`.
In a growing tumour with both daughters retained the ecDNA pool grows in
step with the cell count, but so does the diploid background it is
measured against, and a sequencing library is loaded to a fixed depth
rather than a fixed number of cells. Read support is therefore a
**ratio** — roughly `depth * (junction copies / genome copies) * purity`
— and dividing faster leaves it unchanged. With `selection = 0` this
function reproduces that: copy number stays flat however many
generations are run, while cell-to-cell variance grows and a tail of
cells loses the amplicon entirely. Only `selection > 0` raises the
population mean, at a rate of about
`d log(copies)/d generation ~ selection`. Since bulk WGS sees only that
mean, the per-cell spread is invisible to it — which is why
[`sim_replicate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_replicate.md)'s
deterministic growth is an adequate model for bulk data, and this
function's value is in deriving the level rather than asserting it.

## See also

[`sim_replicate()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_replicate.md),
[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md)

## Examples

``` r
ec <- sim_episome(seed_locus("EGFR"), sample = "SEG01", copies = 20)

## neutral: twenty generations of division leave the level where it started
summary(sim_segregate(ec, generations = 20, selection = 0, seed = 1))$copies
#> [1] 20.821

## under selection it climbs
summary(sim_segregate(ec, generations = 20, selection = 0.1, seed = 1))$copies
#> [1] 56.3165
```
