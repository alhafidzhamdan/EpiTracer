# Linear copy-number and structural rearrangement plot

Draws allele-specific copy number, structural variant arcs, karyotype
ideograms and gene labels for one or more loci laid out side-by-side on
a single concatenated x-axis. Because all loci share one coordinate
system, the structural variants that interconnect separate amplicons
(e.g. the junctions of a multi-fragment / hub ecDNA) are drawn as arcs
spanning the loci.

## Usage

``` r
plot_sv_linear(
  sample,
  cnv_data,
  sv_data,
  wgd_data = NULL,
  genome = c("hg38", "hg19", "mm10"),
  karyotype = NULL,
  gene_coord = NULL,
  chromosome = NULL,
  chromosome_range = NULL,
  loci = NULL,
  events = "amp",
  cluster_gap = 5e+06,
  flank_pct = 10,
  min_cn_ratio = 3,
  gain_ratio = 1.4,
  loh_thresh = 0.5,
  homdel_thresh = 0.5,
  min_amp_width = 1e+05,
  gap_frac = 0.06,
  genes_to_highlight = NULL,
  gene_label_angle = NULL,
  repel_labels = TRUE,
  cn_max = NULL,
  displayExon = FALSE,
  cds_gr = NULL,
  offset_gene = 1.15,
  ymax_highlight_ratio = 1.08,
  karyotype_rel_size = 0.048,
  loh_position_ratio = 0.5,
  highlight_amp = TRUE,
  highlight_hom_del = TRUE,
  amplicons = NULL,
  parallel_breakpoints = NULL,
  highlight_events = NULL,
  highlight_id_col = NULL,
  highlight_colour = "#d95f0e",
  dim_unhighlighted = FALSE,
  wgd_sample_col = NULL,
  snv_data = NULL,
  snv_sample_col = NULL,
  snv_y = c("imd", "vaf", "cn"),
  snv_type_col = NULL,
  vaf_col = "allelic_freq",
  snv_cn_col = "variant_cn",
  snv_timing = FALSE,
  snv_timing_pre_frac = 0.5,
  snv_timing_post_mcn = 1.5,
  snv_timing_colours = c(`Pre-amplification` = "#d1495b", `Post-amplification` =
    "#1d3557", Unknown = "grey70"),
  snv_rel_height = 0.4,
  snv_point_size = 0.7,
  snv_alpha = 0.6,
  snv_colour = "#1d3557",
  vaf_max = 1,
  outdir = NULL,
  save = TRUE,
  plot_width_custom = NULL,
  plot_height_custom = NULL,
  verbose = FALSE
)
```

## Arguments

- sample:

  Character scalar; sample identifier (matched in `cnv_data`, `sv_data`,
  `wgd_data`).

- cnv_data:

  A data.frame of copy-number segments with columns `sample`,
  `seqnames`, `start`, `end`, `copyNumber`, `ploidy`,
  `majorAlleleCopyNumber`, `minorAlleleCopyNumber`.

- sv_data:

  A data.frame of SVs with columns `chrom1`, `start1`, `chrom2`,
  `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`.

- wgd_data:

  Optional data.frame with a sample-identifier column (see
  `wgd_sample_col`) and a `Polyploidy` column (`"No"` = diploid,
  otherwise WGD). If supplied, the sample's WGD status is annotated in
  the plot title; if `NULL` (default) the title shows just the sample
  name.

- genome:

  Genome build for the bundled references: `"hg38"` (default), `"hg19"`
  or `"mm10"`. Selects the karyotype ideogram and oncogene panel used
  when `karyotype` / `gene_coord` are not given. Ignored for either
  reference that is supplied explicitly.

- karyotype:

  Ideogram bands (UCSC `cytoBand` / `chr_info` style), as a data.frame
  (first column = chromosome; needs `start`, `end`, `gieStain`) or an
  `.rds` path. `NULL` (default) uses the bundled reference for `genome`.

- gene_coord:

  Gene coordinates as a data.frame with columns
  `chr`,`start`,`end`,`strand`,`gene`, or a path to a headerless
  tab-separated BED-like file with those columns. `NULL` (default) uses
  the bundled oncogene panel for `genome`; supply your own to label
  other genes.

- chromosome:

  Optional character vector of chromosomes to display, e.g.
  `c("chr7", "chr12")`.

- chromosome_range:

  Optional two-column matrix/data.frame of `start`,`end` window limits,
  one row per entry in `chromosome`.

- loci:

  Optional explicit loci: a `data.frame` with columns
  `chr`,`start`,`end`, or a character vector of `"chr:start-end"`
  strings. Takes precedence over `chromosome`/`chromosome_range`.

- events:

  Character vector of copy-number event types to target when
  auto-detecting loci (used only when neither `loci` nor
  `chromosome_range` is given). Any of: `"amp"` (amplification,
  `copyNumber > min_cn_ratio * ploidy`), `"gain"`
  (`> gain_ratio * ploidy` but not amplified), `"loh"`
  (`minorAlleleCopyNumber < loh_thresh`), `"homdel"` (homozygous
  deletion, `copyNumber < homdel_thresh`). Default `"amp"`. When
  explicit loci are supplied the whole region is plotted regardless of
  event type – the function is a general CN/SV viewer, not
  amplification-only.

- cluster_gap:

  Numeric; when auto-detecting, consecutive event segments more than
  this many bp apart start a new locus, so scattered focal events become
  separate panels instead of one whole-chromosome span (default `5e6`).

- flank_pct:

  Numeric percentage by which each auto-detected (CN-status) region is
  extended on both sides – i.e. the flanking window shown around an
  amplicon / deletion, as a percent of its width (default `10` =
  +/-10%). Only applies to auto-detected loci, not to explicit
  `chromosome`/`loci`.

- min_cn_ratio:

  Numeric; amplification threshold as a multiple of ploidy
  (`copyNumber > min_cn_ratio * ploidy`, default `3`).

- gain_ratio:

  Numeric; gain threshold as a multiple of ploidy (default `1.4`); a
  "gain" is `> gain_ratio * ploidy` but not amplified.

- loh_thresh:

  Numeric; LOH minor-allele threshold
  (`minorAlleleCopyNumber < loh_thresh`, default `0.5`).

- homdel_thresh:

  Numeric; homozygous-deletion copy-number threshold
  (`copyNumber < homdel_thresh`, default `0.5`).

- min_amp_width:

  Numeric; drop auto-detected loci whose total event span is below this
  many bp (default `1e5`).

- gap_frac:

  Numeric; gap between loci as a fraction of the total plotted width
  (default `0.06`).

- genes_to_highlight:

  Optional character vector of gene symbols. If `NULL`, a default
  oncogene panel is used; only genes falling inside a locus window are
  drawn.

- gene_label_angle:

  Optional numeric label rotation in degrees. If `NULL` (default) labels
  are horizontal, switching to 45 degrees automatically when genes are
  crowded.

- repel_labels:

  Logical; use ggrepel to de-collide gene labels (default `TRUE`; falls
  back to plain labels if ggrepel is absent).

- cn_max:

  Optional numeric override for the copy-number axis top; if `NULL` it
  is rounded up from the data to two significant figures, with a floor
  of 2 so near-diploid regions still use a full 0-2 axis.

- displayExon:

  Logical; if `TRUE`, draw exon models (from `cds_gr`) for in-window
  genes instead of a point + label (default `FALSE`).

- cds_gr:

  Optional
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of CDS/exon ranges (metadata column `gene_name`); required when
  `displayExon = TRUE`.

- offset_gene, ymax_highlight_ratio, karyotype_rel_size,
  loh_position_ratio:

  Layout tuning parameters.

- highlight_amp, highlight_hom_del:

  Logical; shade amplified / homozygously deleted segments.

- amplicons:

  Optional distinct-amplicon overlay: a
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  or data.frame of amplicon regions (columns `seqnames`/`chr`, `start`,
  `end`, and an optional `ID`/`label`). Each amplicon that falls in the
  plotted window is drawn as a short horizontal bar across the TOP of
  the plot spanning its extent, in its own colour and with its label
  above, so several distinct amplicons (e.g. a focal episome sitting
  inside a larger, separately detected amplified span) are visually
  separable without washing over the data. Genomic coordinates are
  mapped through the same per-locus transform as the rest of the plot,
  so the bars align with the copy-number track.

- parallel_breakpoints:

  Optional data.frame of *adjacent parallel breakpoint* pairs to
  highlight (the breakage-replication/fusion hallmark): columns
  `chr`/`seqnames`, `pos1`, `pos2` (the two same-orientation breakends
  of a pair) and an optional `strand`. Each pair is drawn near the
  baseline as a bracket joining the two breakends, with a caret marking
  each, so BRF breakpoints stand out among the other junctions.

- highlight_events:

  Optional character vector of SV identifiers to draw bold on top of the
  rest (e.g. the junctions of one chromoplexy cycle, a TBA boundary, a
  fold-back set). Matched against `highlight_id_col`.

- highlight_id_col:

  Optional name of the column in `sv_data` holding the identifiers
  matched by `highlight_events`; `NULL` (default) auto-detects `name`
  then `event`.

- highlight_colour:

  Colour for highlighted SVs (default `"#d95f0e"`).

- dim_unhighlighted:

  Logical; when `TRUE`, non-highlighted SVs are greyed so the
  highlighted set stands out (default `FALSE`).

- wgd_sample_col:

  Optional name of the sample column in `wgd_data` (default: `sample`,
  falling back to `WGS_ID`).

- snv_data:

  Optional SNV/SSM table as a `data.frame` or
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  with (at least) `seqnames`/`start` position columns, a
  sample-identifier column (see `snv_sample_col`) and, for
  `snv_y = "vaf"`, a variant-allele-frequency column (see `vaf_col`).
  When supplied, a second SNV panel (intermutation-distance rainfall by
  default, or VAF; see `snv_y`) is drawn directly beneath the
  copy-number / SV panel, sharing the same concatenated genomic x-axis,
  and the two are returned/saved as one stacked figure (requires the
  patchwork package). When `NULL` (default) only the CN/SV panel is
  drawn and behaviour is unchanged.

- snv_sample_col:

  Optional name of the sample column in `snv_data` (default: `sampleID`,
  falling back to `sample`).

- snv_y:

  What the SNV panel's y-axis shows: `"imd"` (default) plots the
  intermutation distance – the bp distance to the previous SNV on the
  same chromosome, a rainfall plot, on a log10 axis; `"vaf"` plots the
  variant-allele frequency; `"cn"` plots the SNV copy number (mutation
  copy number, see `snv_cn_col`) – in an amplicon this times each SNV
  against the amplification. Only single-nucleotide variants are used in
  every case (indels / MNVs are excluded, see `snv_type_col`);
  intermutation distances are computed across all of the sample's SNVs
  per chromosome before restricting to the plotted loci, so window-edge
  mutations keep their true neighbour distance.

- snv_type_col:

  Optional name of a mutation-type column in `snv_data` used to keep
  SNVs only (rows where the value is `"SNV"`). Defaults to `type` when
  present; if no such column exists, SNVs are inferred from single-base
  `ref`/`mut` columns.

- vaf_col:

  Name of the VAF column in `snv_data` (default `allelic_freq`), used
  when `snv_y = "vaf"`. Values outside `[0, vaf_max]` are treated as
  artefacts and dropped.

- snv_cn_col:

  Name of the SNV copy-number column in `snv_data` (default
  `variant_cn`), used when `snv_y = "cn"` and for `snv_timing`.

- snv_timing:

  Logical; if `TRUE`, classify each SNV by its timing relative to the
  focal amplification and colour the points accordingly, with a legend
  collected to the right of the figure. Uses the mutation copy number
  (`snv_cn_col`) against the amplified-allele copy number (`major_cn`)
  at the site (also needs `minor_cn`): a site is amplified when
  `major_cn + minor_cn > min_cn_ratio * ploidy`; within it, an SNV is
  "Pre-amplification" when its copy number is
  `>= snv_timing_pre_frac * major_cn` (and `>= 2`), "Post-amplification"
  when it is `<= snv_timing_post_mcn`, and "Unknown" otherwise or when
  the site is not amplified. Works with any `snv_y`.

- snv_timing_pre_frac:

  Numeric; fraction of the amplified-allele copy number an SNV's copy
  number must reach to be called pre-amplification (default `0.5`).

- snv_timing_post_mcn:

  Numeric; SNV copy-number at or below which a variant in an amplified
  site is called post-amplification (default `1.5`).

- snv_timing_colours:

  Named character vector of point colours for the `"Pre-amplification"`,
  `"Post-amplification"` and `"Unknown"` classes.

- snv_rel_height:

  Numeric; height of the SNV panel relative to the CN/SV panel (default
  `0.4`).

- snv_point_size, snv_alpha, snv_colour:

  Point size, alpha and colour for the SNV VAF scatter (defaults `0.7`,
  `0.6`, `"#1d3557"`).

- vaf_max:

  Numeric top of the VAF axis; also the upper bound of the plausible VAF
  window used to drop artefactual values (default `1`).

- outdir:

  Optional directory in which to write the plot. If `NULL` no file is
  written and only the plot object is returned.

- save:

  Logical; if `TRUE` (default) and `outdir` is supplied, write the plot
  as a PDF.

- plot_width_custom, plot_height_custom:

  Optional numeric overrides for the output dimensions (inches).

- verbose:

  Logical; print progress/diagnostic messages.

## Value

The assembled
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object (invisibly when a file is written); written paths are attached as
attribute `"path"`.

## Details

Optionally, supplying `snv_data` adds a second panel of the sample's
small mutations directly beneath the copy-number / SV panel, sharing the
same concatenated genomic x-axis so each mutation lines up with the copy
number and rearrangements it sits within. Its y-axis shows the
intermutation distance (a rainfall plot, the default), the
variant-allele frequency, or the SNV (mutation) copy number (see
`snv_y`), and the points can be coloured by their timing relative to the
focal amplification (`snv_timing`). The two panels are returned and
saved as one stacked figure (requires patchwork); see the `snv_*`
parameters for the full set of controls.

The loci to display are resolved in this order:

1.  `loci` if supplied (explicit windows);

2.  `chromosome` + `chromosome_range` (explicit windows, one per
    chromosome);

3.  `chromosome` alone – the amplified region on each named chromosome
    is auto-detected and padded by `flank_pct`; a chromosome with no
    amplification falls back to its whole length;

4.  nothing – every amplified locus in the sample is detected
    automatically.

A locus is "amplified" where `copyNumber > min_cn_ratio * ploidy`; its
extent is the min-max of amplified segments on the chromosome, padded by
`flank_pct`, with single-segment artefacts below `min_amp_width`
dropped.

## See also

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single focused locus:
plot_sv_linear("DO11441T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
               chromosome = "chr7",
               chromosome_range = matrix(c(52e6, 56e6), nrow = 1),
               outdir = "plots")

# All amplified loci, interconnected (multi-fragment ecDNA hub):
plot_sv_linear("DUMC12T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
               outdir = "plots")

# Explicit loci:
plot_sv_linear("S1", cnv, sv, wgd, karyotype = K, gene_coord = G,
               loci = c("chr4:50e6-64e6", "chr12:57e6-59e6"))
} # }
```
