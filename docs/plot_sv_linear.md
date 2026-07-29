[← EpiTracer](../README.md)

# `plot_sv_linear()` — general-purpose flexible SV + CNV plotter

A generic, all-purpose plotter for **structural variants and copy number**. It
draws one or more loci side-by-side on a single concatenated x-axis with
allele-specific copy-number tracks, SV arcs, karyotype ideograms, LOH /
homozygous-deletion bars and gene labels (saved as PDF). Point it at a single
focused locus (`chromosome` + `chromosome_range`), give explicit `loci`, or let
it **auto-detect** copy-number events (amplification, gain, LOH, homozygous
deletion) — the structural variants interconnecting separate amplicons
(multi-fragment / hub ecDNA junctions) are drawn as arcs spanning the loci.

Although developed alongside the caller, it is **not ecDNA-specific**: use it to
visualise any CNV / SV landscape at a single locus, across whole chromosomes, or
genome-wide.

Optionally, supplying `snv_data` adds a **second SNV panel** directly beneath
the CN / SV panel, sharing the same concatenated genomic x-axis so each mutation
lines up with the copy number and rearrangements it sits within. Its y-axis
shows the intermutation distance (a rainfall plot, the default), the
variant-allele frequency, or the SNV (mutation) copy number (`snv_y`), and the
points can be coloured by their timing relative to the focal amplification
(`snv_timing`). The two panels are saved as one stacked figure (requires the
[patchwork](https://patchwork.data-imaginist.com) package).

## Example plots

**Single locus** — an episomal *EGFR* ecDNA in `DO11441T1`: a boundary
duplication (magenta arc reaching the top) with its **excision-scar** deletion
(light-blue arc below the baseline), over an ~100-copy amplicon.

![EGFR episomal ecDNA in DO11441T1](../man/figures/README-DO11441T1-EGFR.png)

**Multiple loci** — a multi-fragment co-amplification hub in `DUMC12T1`:
*PDGFRA* (chr4), a chr9p amplicon and *CDK4* / *MDM2* (chr12), auto-detected and
drawn side-by-side on one axis with the structural variants interconnecting them.

![Multi-chromosome amplicon hub in DUMC12T1](../man/figures/README-DUMC12T1-multichrom.png)

## Input data

Uses PURPLE CNV segments (`cnv_data`) and PURPLE SV / BEDPE breakpoints
(`sv_data`) for one sample. Karyotype ideogram and gene coordinates default to
small bundled hg38 references, so a minimal call needs only `sample`, `cnv_data`
and `sv_data`.

## Usage

```r
library(EpiTracer)

# A single focused locus (karyotype / gene_coord default to bundled hg38 refs;
# wgd_data is optional — supply it to annotate WGD status in the title):
plot_sv_linear(
  sample     = "DO11441T1",
  cnv_data   = cnv_df, sv_data = sv_df,
  chromosome = "chr7",
  chromosome_range = matrix(c(52e6, 56e6), nrow = 1),
  outdir     = "plots"
)

# Every amplified locus, interconnected (multi-fragment ecDNA hub) — auto-detected:
plot_sv_linear(
  sample    = "DUMC12T1",
  cnv_data  = cnv_df, sv_data = sv_df, wgd_data = wgd_df,
  outdir    = "plots",   # loci auto-detected; or events = "homdel", loci = c("chr4:...", ...)
  flank_pct = 10         # extend each auto-detected window by ±10%
)

# With a stacked SNV panel beneath the CN/SV plot (rainfall by default;
# snv_y = "vaf" or "cn"; add snv_timing = TRUE to colour by amplification timing):
plot_sv_linear(
  sample     = "DO11441T1",
  cnv_data   = cnv_df, sv_data = sv_df,
  chromosome = "chr7",
  chromosome_range = matrix(c(52e6, 56e6), nrow = 1),
  snv_data   = snv_df,   # needs seqnames/start, a sample column, and a VAF column for snv_y = "vaf"
  snv_y      = "imd",
  outdir     = "plots"
)
```

## Arguments

Only `sample`, `cnv_data` and `sv_data` are **required**; everything else is
optional.

**Required**

| Argument | Description |
| :--- | :--- |
| `sample`   | Sample identifier (matched in `cnv_data` / `sv_data`). |
| `cnv_data` | PURPLE CN segments: `sample`, `seqnames`, `start`, `end`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`. |
| `sv_data`  | PURPLE SVs / BEDPE: `chrom1`, `start1`, `chrom2`, `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`. |

**Optional — references & sample metadata**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `wgd_data`       | `NULL` | WGD table (`Polyploidy` column); annotates WGD status in the title when supplied. |
| `karyotype`      | bundled hg38 | Ideogram bands (data.frame or `.rds` path). |
| `gene_coord`     | bundled hg38 oncogenes | Gene coordinates (data.frame or BED path). |
| `wgd_sample_col` | `NULL` | Sample-id column in `wgd_data` (`sample`, else `WGS_ID`). |

**Optional — which loci to plot**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `chromosome`       | `NULL`  | Chromosome(s) to display. |
| `chromosome_range` | `NULL`  | Two-column start/end window(s), one row per `chromosome`. |
| `loci`             | `NULL`  | Explicit loci: data.frame (`chr`,`start`,`end`) or `"chr:start-end"` strings; overrides `chromosome`. |
| `events`           | `"amp"` | Auto-detect target(s): `"amp"`, `"gain"`, `"loh"`, `"homdel"`. |
| `cluster_gap`      | `5e6`   | Events more than this many bp apart become separate loci. |
| `flank_pct`        | `10`    | Extend each auto-detected window by ±X% of its width. |
| `min_amp_width`    | `1e5`   | Drop auto-detected loci with total event span below this (bp). |

**Optional — event thresholds**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `min_cn_ratio`  | `3`   | Amplification: `copyNumber > min_cn_ratio × ploidy`. |
| `gain_ratio`    | `1.4` | Gain: `> gain_ratio × ploidy`, but not amplified. |
| `loh_thresh`    | `0.5` | LOH: `minorAlleleCopyNumber < loh_thresh`. |
| `homdel_thresh` | `0.5` | Homozygous deletion: `copyNumber < homdel_thresh`. |

**Optional — gene labels**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `genes_to_highlight` | `NULL`  | Gene symbols to label (default oncogene panel if `NULL`). |
| `gene_label_angle`   | `NULL`  | Label rotation in degrees; auto-angles when crowded if `NULL`. |
| `repel_labels`       | `TRUE`  | De-collide labels with ggrepel. |
| `displayExon`        | `FALSE` | Draw exon models (requires `cds_gr`). |
| `cds_gr`             | `NULL`  | GRanges of CDS/exon ranges (`gene_name`), for `displayExon`. |

**Optional — SNV panel**

Supply `snv_data` to add the stacked SNV panel; all other `snv_*` arguments are
ignored when `snv_data` is `NULL`.

| Argument | Default | Description |
| :--- | :--- | :--- |
| `snv_data`            | `NULL`          | SNV/SSM table (data.frame or GRanges) with `seqnames`/`start`, a sample column, and (for `snv_y = "vaf"`) a VAF column. Enables the panel. |
| `snv_sample_col`      | `NULL`          | Sample-id column in `snv_data` (`sampleID`, else `sample`). |
| `snv_y`               | `"imd"`         | Panel y-axis: `"imd"` (intermutation distance / rainfall), `"vaf"`, or `"cn"` (mutation copy number). |
| `snv_type_col`        | `NULL`          | Mutation-type column, used to keep SNVs only (`type`; else inferred from single-base `ref`/`mut`). |
| `vaf_col`             | `"allelic_freq"`| VAF column, used when `snv_y = "vaf"`. |
| `snv_cn_col`          | `"variant_cn"`  | SNV copy-number column, used when `snv_y = "cn"` and for `snv_timing`. |
| `vaf_max`             | `1`             | Top of the VAF axis; also the upper bound for dropping artefactual VAFs. |
| `snv_timing`          | `FALSE`         | Colour SNVs by timing relative to the amplification (pre- / post- / unknown); needs `snv_cn_col`, `major_cn`, `minor_cn`. |
| `snv_timing_pre_frac` | `0.5`           | Fraction of the amplified-allele CN an SNV must reach to be called pre-amplification. |
| `snv_timing_post_mcn` | `1.5`           | SNV copy-number at or below which a variant in an amplified site is called post-amplification. |
| `snv_timing_colours`  | 3 colours       | Point colours for the `Pre-amplification` / `Post-amplification` / `Unknown` classes. |
| `snv_rel_height`      | `0.4`           | SNV panel height relative to the CN/SV panel. |
| `snv_point_size`      | `0.7`           | SNV point size. |
| `snv_alpha`           | `0.6`           | SNV point alpha. |
| `snv_colour`          | `"#1d3557"`     | SNV point colour (when not timing-coloured). |

**Optional — layout & appearance**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `cn_max`               | `NULL`  | Copy-number axis top (auto from data if `NULL`). |
| `gap_frac`             | `0.06`  | Gap between loci as a fraction of total width. |
| `offset_gene`          | `1.15`  | Gene-label height relative to max CN. |
| `ymax_highlight_ratio` | `1.08`  | Height of amp / homdel shading relative to max CN. |
| `karyotype_rel_size`   | `0.048` | Ideogram height relative to the CN axis. |
| `loh_position_ratio`   | `0.5`   | LOH / homdel bar position within the ideogram gap. |
| `highlight_amp`        | `TRUE`  | Shade amplified segments. |
| `highlight_hom_del`    | `TRUE`  | Shade homozygously deleted segments. |
| `plot_width_custom`    | `NULL`  | Override output width (inches). |
| `plot_height_custom`   | `NULL`  | Override output height (inches). |

**Optional — output**

| Argument | Default | Description |
| :--- | :--- | :--- |
| `outdir`  | `NULL`  | Directory to write the PDF; if `NULL`, only the plot object is returned. |
| `save`    | `TRUE`  | Write the PDF when `outdir` is supplied. |
| `verbose` | `FALSE` | Print progress / diagnostics. |

Returns the `ggplot` object (invisibly when a file is written); any written PDF
path is attached as `attr(p, "path")`.
