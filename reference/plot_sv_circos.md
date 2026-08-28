# Whole-genome SV circos plot for one sample

Draws a per-sample circos: chromosome-label + banded ideogram, a
copy-number status ring (amplification / gain / LOH / homozygous
deletion / neutral, relative to ploidy), separate intra- and
inter-chromosomal breakpoint-density rings, and structural-variant links
coloured by class (DUP / TRA / h2hINV / t2tINV / other). Titled with the
sample and its whole-genome-doubling status.

## Usage

``` r
plot_sv_circos(
  sample,
  sv_data,
  cnv_data,
  wgd_data = NULL,
  genome = c("hg38", "hg19", "mm10"),
  bin_size = 1e+06,
  highlight_events = NULL,
  highlight_id_col = NULL,
  highlight_colour = "#d95f0e",
  dim_unhighlighted = FALSE,
  outdir = NULL,
  overwrite = TRUE
)
```

## Arguments

- sample:

  Character scalar; the sample to plot (matched against the `sample`
  column of `sv_data` / `cnv_data`).

- sv_data:

  A `data.frame`/`data.table` of structural variants in BEDPE layout:
  `sample`, `chrom1`, `start1`, `chrom2`, `start2`, `svclass` (and
  optionally `end1`/`end2`; when absent the start positions are used).

- cnv_data:

  A `data.frame`/`data.table` of copy-number segments: `sample`,
  `seqnames`, `start`, `end`, `copyNumber`, `ploidy`,
  `minorAlleleCopyNumber`.

- wgd_data:

  Optional `data.frame` with `sample` and `Polyploidy` (`"wgd"` vs
  anything else) for the title; `NULL` labels the sample only.

- genome:

  Reference build; one of `"hg38"` (default), `"hg19"`, `"mm10"`.

- bin_size:

  Integer; genome tile width (bp) for the density and copy-number rings
  (default `1e6`).

- highlight_events:

  Optional character vector of SV identifiers to draw bold on top of the
  faint full-genome links (e.g. the junctions of one chromoplexy cycle).
  Matched against `highlight_id_col`.

- highlight_id_col:

  Optional name of the column in `sv_data` holding the identifiers
  matched by `highlight_events`; `NULL` (default) auto-detects `name`
  then `event`.

- highlight_colour:

  Colour for highlighted links (default `"#d95f0e"`).

- dim_unhighlighted:

  Logical; when `TRUE`, non-highlighted links are greyed so the
  highlighted set stands out (default `FALSE`).

- outdir:

  Optional directory; when supplied the plot is written to
  `<outdir>/<sample>_circos.pdf`. When `NULL` the plot is drawn on the
  current graphics device (open one yourself, e.g. with
  [`grDevices::pdf()`](https://rdrr.io/r/grDevices/pdf.html)).

- overwrite:

  Logical; when `FALSE` and the output PDF already exists, the sample is
  skipped (default `TRUE`).

## Value

Invisibly `NULL`; called for its plotting side effect.

## Details

A port of a bespoke circos plotter to EpiTracer conventions: hg38 by
default, the ideogram and genome tiling built from the bundled
`chr_info`, and the gUtils/regioneR helpers replaced by the package's
own
[`to_granges()`](https://alhafidzhamdan.github.io/EpiTracer/reference/to_granges.md)
/ `\%$\%`. `chrY` and `chrM` are dropped. Requires the circlize package.

## See also

[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)
