# Resolve gene symbol(s) to a plotting locus

Looks a gene symbol up in the bundled oncogene panel for `genome` and
returns a window suitable for
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)'s
`loci` argument, so a user can ask to plot e.g. `"EGFR"` or `"CDKN2A"`
by name. When `cnv` is supplied the window is grown to the copy-number
event (amplification or homozygous deletion) overlapping the gene, so
the zoom fits the event; otherwise the gene span plus `flank` is used.

## Usage

``` r
gene_locus(
  genes,
  genome = "hg38",
  cnv = NULL,
  sample = NULL,
  flank = 5e+05,
  target = c("any", "amp", "homdel"),
  min_cn_ratio = 3,
  homdel_thresh = 0.5
)
```

## Arguments

- genes:

  Character vector of gene symbols (case-insensitive).

- genome:

  One of `"hg38"`, `"hg19"`, `"mm10"`.

- cnv:

  Optional copy-number table (from
  [`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md))
  to size the window to the overlapping event.

- sample:

  Optional sample to subset `cnv` by.

- flank:

  Bp added on each side of the resolved region (default `5e5`).

- target:

  `"amp"` sizes to an amplified segment, `"homdel"` to a homozygous
  deletion, `"any"` (default) to the gene span.

- min_cn_ratio, homdel_thresh:

  Amplification (`copyNumber > min_cn_ratio * ploidy`) and
  homozygous-deletion (`copyNumber < homdel_thresh`) thresholds.

## Value

A data.frame with `chr`, `start`, `end`, `gene` (one row per found
gene).

## See also

[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md),
[`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md)
