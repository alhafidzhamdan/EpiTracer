# Temporally-stratified rearrangement plot (SV+CN grouped by read support)

A companion to
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)
that decomposes a locus's structural rearrangements into read-support
(variant-fraction, `VF`) strata and draws one stacked panel per stratum.
Junctions with the highest read support – taken as a proxy for the
earliest / highest-copy events – are plotted at the top, then
successively lower-`VF` strata below, so the figure reads as a rough
temporal decomposition of amplicon evolution.

## Usage

``` r
plot_sv_reconstruction(
  sample,
  cnv_data,
  sv_data,
  wgd_data = NULL,
  karyotype = system.file("extdata", "chr_info_hg38.rds", package = "EpiTracer"),
  gene_coord = system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
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
  vf_col = "VF",
  k = "auto",
  max_k = 4,
  vf_breaks = NULL,
  min_vf = 1,
  isolate_founder = TRUE,
  founder_offscale = FALSE,
  vf_scale = c("shared", "per_panel"),
  cn_near_flank = 1e+05,
  cn_display = c("reconstruct", "actual"),
  prior_sv_alpha = 0.15,
  founder_alpha = 0.5,
  cn_border_lines = TRUE,
  cn_border_min_step = NULL,
  drop_empty_strata = TRUE,
  panel_rel_height = 1.4,
  genes_to_highlight = NULL,
  gene_label_angle = NULL,
  repel_labels = TRUE,
  displayExon = FALSE,
  cds_gr = NULL,
  cn_max = NULL,
  offset_gene = 1.15,
  ymax_highlight_ratio = 1.08,
  karyotype_rel_size = 0.048,
  loh_position_ratio = 0.5,
  highlight_amp = TRUE,
  highlight_hom_del = TRUE,
  wgd_sample_col = NULL,
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

- karyotype:

  Ideogram bands (UCSC `cytoBand` / `chr_info` style), as a data.frame
  (first column = chromosome; needs `start`, `end`, `gieStain`) or an
  `.rds` path. Defaults to the bundled hg38 reference.

- gene_coord:

  Gene coordinates as a data.frame with columns
  `chr`,`start`,`end`,`strand`,`gene`, or a path to a headerless
  tab-separated BED-like file with those columns. Defaults to a small
  bundled hg38 table covering the default oncogene panel; supply your
  own for other genes or genome builds.

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

- vf_col:

  Name of the read-support column in `sv_data` (default `"VF"`).

- k:

  Number of `VF` strata. `"auto"` (default) chooses the number of
  clusters automatically (via Ckmeans.1d.dp BIC when available,
  otherwise a within-cluster-sum-of-squares elbow capped at `max_k`); an
  integer forces that many. Reduced automatically when a sample has
  fewer distinct `VF` values than `k`.

- max_k:

  Upper bound on the number of strata when `k = "auto"` (default `4`).

- vf_breaks:

  Optional numeric vector of explicit `VF` cut points (upper edges),
  e.g. `c(200, 80)` for strata `>200`, `80-200`, `<80`. Overrides `k`.

- min_vf:

  Junctions with `VF < min_vf` are dropped before clustering (default
  `1`).

- isolate_founder:

  Logical (default `TRUE`); pull the single highest-`VF` junction (the
  defining excision) out of the clustering and give it its own top
  panel, headed "Max VF". Its arc is drawn at the same weight as every
  other junction. It is excluded from the k-means / breaks step (so one
  dominant outlier no longer skews the clusters), and the remaining
  strata are headed "Cluster 1..X" (highest `VF` first). Ties at the
  maximum `VF` are treated as one group. See `founder_offscale` for how
  it relates to the axis.

- founder_offscale:

  Logical (default `FALSE`); controls the shared read-support axis when
  a founder is isolated. `FALSE`: the axis reflects the true maximum
  `VF` (founder included), so the founder arc's height is literal.
  `TRUE`: the axis is driven by the highest non-founder junction (so the
  lower strata are not compressed by a large founder outlier) and the
  founder arc is clamped to the axis top – its true `VF` is still
  labelled, and every panel's read-support axis title is marked
  "(scaled)" to flag that the axis no longer maps 1:1 to `VF`. Preferred
  when the founder `VF` greatly exceeds the next highest.

- vf_scale:

  How arc heights (and the right-hand read-support axis) are scaled.
  `"shared"` (default) puts every panel on one global
  read-support-to-copy-number scale, so arc heights are directly
  comparable across strata – best for multi-locus / hub amplicons.
  `"per_panel"` scales each panel to its own stratum's `VF` range, so a
  single very-high-`VF` junction does not crush the dynamic range of the
  lower strata – usually the better choice for a dense single-chromosome
  amplicon. The left copy-number axis and the x-axis stay shared either
  way, so panels remain aligned.

- cn_near_flank:

  Numeric; in `cn_display = "actual"`, a copy-number segment is drawn in
  a panel if it lies within this many bp of any of that stratum's
  breakpoints (and, for intra-locus junctions, if it lies between the
  two breakpoints). Default `1e5`.

- cn_display:

  How the major-allele copy number is drawn per panel. `"reconstruct"`
  (default) shows an *iterative* reconstruction of how the copy-number
  profile is built up wave by wave. The breakpoints of the junctions
  introduced so far (strata \\1..k\\) partition each locus into
  intervals, and every interval is drawn flat at the minimum of the two
  copy-number segments immediately adjacent to (just inside) its
  bounding breakpoints. The founder panel therefore shows a single flat
  baseline over the founder span; each lower panel introduces more
  breakpoints, so the intervals subdivide and copy-number structure
  emerges; and the bottom panel is the full observed bulk copy number
  that the reconstruction evolves toward. This is a visual model, not a
  formal amplicon deconvolution. `"actual"` instead draws the real
  observed CN segments near each stratum's breakpoints in every panel
  (see `cn_near_flank`). The minor-allele CN is always drawn in full
  regardless of this setting. In `"reconstruct"` the allele-loss bars
  (LOH where minor CN `< loh_thresh`, and homozygous deletions where
  total CN `< homdel_thresh`) are placed on the same timeline: a loss
  segment is revealed in the wave of the earliest intra-locus deletion
  junction whose deleted span covers it, or – if no deletion covers it –
  from the founder wave down (treated as pre-existing).

- prior_sv_alpha:

  Numeric in `[0, 1]`; opacity at which the junctions of earlier
  (higher-`VF`) strata are re-drawn faintly in each lower panel, so the
  accumulation of rearrangements is visible going down. The current
  stratum is always drawn at full opacity. Set to `0` to show only the
  current stratum's junctions in each panel (default `0.15`).

- founder_alpha:

  Numeric in `[0, 1]`; opacity at which the founder junction (the
  defining, highest-`VF` event) is re-drawn as a prior in the lower
  panels. Kept higher than `prior_sv_alpha` so the founder stays
  trackable all the way down the reconstruction, while later waves fade
  more. The founder arc is also drawn with a slightly bolder line than
  the other junctions. Has no effect on the founder's own (top) panel,
  where it is always at full opacity (default `0.5`).

- cn_border_lines:

  Logical (default `TRUE`); draw vertical dashed guide lines at the
  borders of the final (fully reconstructed) copy-number segments.
  Because the panels share the x-axis, the guides line up across the
  whole stack, marking where copy-number transitions occur.

- cn_border_min_step:

  Numeric minimum copy-number change for an *internal* border line to be
  drawn (each locus's outer edges are always drawn). If `NULL` (default)
  a proportional threshold (`max(0.25, 0.015 * cn_axis_max)`) is used,
  so a guide is drawn at essentially every perceptible copy-number step
  (only flat, equal-level segment splits are skipped). Set larger to
  show only the biggest transitions, smaller to also mark sub-copy
  noise.

- drop_empty_strata:

  Logical; drop strata that end up with no drawable in-locus junction
  (default `TRUE`).

- panel_rel_height:

  Numeric relative height of the (label-bearing) bottom panel versus the
  others (default `1.4`).

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

- displayExon:

  Logical; if `TRUE`, draw exon models (from `cds_gr`) for in-window
  genes instead of a point + label (default `FALSE`).

- cds_gr:

  Optional
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of CDS/exon ranges (metadata column `gene_name`); required when
  `displayExon = TRUE`.

- cn_max:

  Optional numeric override for the copy-number axis top; if `NULL` it
  is rounded up from the data to two significant figures, with a floor
  of 2 so near-diploid regions still use a full 0-2 axis.

- offset_gene, ymax_highlight_ratio, karyotype_rel_size,
  loh_position_ratio:

  Layout tuning parameters.

- highlight_amp, highlight_hom_del:

  Logical; shade amplified / homozygously deleted segments.

- wgd_sample_col:

  Optional name of the sample column in `wgd_data` (default: `sample`,
  falling back to `WGS_ID`).

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

A patchwork object stacking the per-stratum panels (invisibly when a
file is written); the written path is attached as attribute `"path"`,
and the per-junction stratum assignment as attribute `"strata"`.

## Details

Strata are found by one-dimensional clustering of `log(VF)` (k-means;
the optimal 1-D solver Ckmeans.1d.dp is used when installed, otherwise
[stats::kmeans](https://rdrr.io/r/stats/kmeans.html)). `VF` here is the
junction read-support count (as emitted by PURPLE), which is strongly
right-skewed, so clustering is done on the log scale – equal-width `VF`
bins would collapse almost every junction into the lowest bracket. All
panels share the same concatenated-locus x-axis and the same copy-number
/ read-support y-scaling, so arc heights are comparable across strata;
each panel draws only the copy-number segments lying near that stratum's
breakpoints, emphasising the CN change each set of events bounds.

Loci are resolved exactly as in
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)
(explicit `loci`; `chromosome` + `chromosome_range`; `chromosome` alone;
or auto-detection of all amplified loci).

## See also

[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md),
[`call_episomal_ecdna()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_ecdna.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plot_sv_reconstruction("DUMC12T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
              k = "auto", outdir = "plots")

# Fixed read-support brackets:
plot_sv_reconstruction("DUMC12T1", cnv, sv, karyotype = K, gene_coord = G,
              vf_breaks = c(200, 80))
} # }
```
