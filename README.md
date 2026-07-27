# EpiTracer

<!-- badges: start -->
<!-- badges: end -->

**EpiTracer: Calling and visualising extrachromosomal circular DNA amplicons
likely generated from simple excision events.**

It detects such amplicons from whole-genome sequencing (WGS) data and visualises
the structural rearrangements underlying focal amplifications.

The package currently provides two functions:

- **`call_episomal_ecdna()`** — an episomal ecDNA caller. For each ecDNA
  amplicon it locates the structural variant breakpoints at the amplicon
  boundaries and flags amplicons whose structure is consistent with the
  *episome* model of formation: a circular amplicon bounded by a duplication
  (DUP) breakpoint, arising from an otherwise non-amplified chromosomal region,
  often leaving a deletion "excision scar" at the origin locus.
- **`plot_sv_linear()`** — a linear allele-specific copy-number and structural
  rearrangement plotter. It draws one or more loci side-by-side on a single
  concatenated x-axis with CN tracks, SV arcs, karyotype ideograms, LOH /
  homozygous-deletion bars, and gene labels (saved as PDF). Point it at a
  single focused locus (`chromosome` + `chromosome_range`), give explicit
  `loci`, or let it **auto-detect every amplified locus** in the sample — the
  structural variants interconnecting separate amplicons (multi-fragment / hub
  ecDNA junctions) are drawn as arcs spanning the loci.

## Installation

Install with `BiocManager`, which resolves the Bioconductor
(`GenomicRanges`, `regioneR`), GitHub ([`gUtils`](https://github.com/mskilab/gUtils))
and `Remotes:` dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("alhafidzhamdan/EpiTracer")
```

## The episome heuristic

For each amplicon, `call_episomal_ecdna()`:

1. finds DUP breakpoints at the amplicon boundaries that are themselves
   amplified (`PURPLE_CN > 3 × ploidy`);
2. requires the boundary DUP to carry the **highest variant fraction (VF)** of
   any DUP in the amplicon;
3. requires the chromosomal segments immediately **flanking both boundaries to
   be non-gained** (consistent with a circle excised from a diploid region);
4. flags a shared flanking **deletion** as a candidate **excision scar**.

## Usage

```r
library(EpiTracer)

# Inputs are GRanges objects (e.g. derived from AmpliconArchitect + PURPLE):
episomal <- call_episomal_ecdna(
  ecdna_gr        = ecdna_gr,        # amplicon regions; needs $ID, $WGS_ID
  breakpoints_gr  = breakpoints_gr,  # SV breakpoints (PURPLE-style columns)
  cnv_gr          = cnv_gr,          # allele-specific CN segments
  cancer_genes_gr = cancer_genes_gr, # cancer gene loci for annotation
  ext             = 1e7,
  mc.cores        = 4
)

# A single focused locus (karyotype/gene_coord default to bundled hg38 refs;
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
  sample     = "DUMC12T1",
  cnv_data   = cnv_df, sv_data = sv_df, wgd_data = wgd_df,
  outdir     = "plots",   # loci auto-detected; or events = "homdel", loci = c("chr4:...", ...)
  flank_pct  = 10         # extend each auto-detected window by +/-10%
)
```

## Data requirement

**Copy-number (CNV) and structural variant (SV) data must come from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple), part of
the Hartwig Medical Foundation (HMF) pipeline.** EpiTracer relies on PURPLE's
allele-specific copy-number segments and its SV/breakpoint output (variant
fraction, junction copy number, PURPLE copy number, etc.); other callers'
outputs are not supported unless coerced to the columns below.

## Input format

`call_episomal_ecdna()` expects these PURPLE/HMF metadata columns:

| Object            | Required metadata |
|-------------------|-------------------|
| `ecdna_gr`        | `ID`, `WGS_ID` |
| `breakpoints_gr`  | `WGS_ID`, `event`, `svclass`, `PURPLE_AF`, `PURPLE_JCN`, `VF`, `PURPLE_CN`, `insLen`, `HOMLEN` |
| `cnv_gr`          | `sample`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber` |
| `cancer_genes_gr` | any (used for overlap annotation) |

`plot_sv_linear()` uses PURPLE CNV segments (`cnv_data`) and PURPLE SV/BEDPE
breakpoints (`sv_data`) for the same sample.

## Arguments

### `call_episomal_ecdna()`

Detects episomal ecDNA amplicons. **Required** arguments have no default.

| Argument | Required | Default | Description |
|---|:---:|---|---|
| `ecdna_gr`        | **yes** | — | GRanges of ecDNA amplicon regions; needs `ID`, `WGS_ID`. |
| `breakpoints_gr`  | **yes** | — | GRanges of PURPLE SV breakpoints (see *Input format*). |
| `cnv_gr`          | **yes** | — | GRanges of PURPLE allele-specific CN segments. |
| `cancer_genes_gr` | **yes** | — | GRanges of cancer-gene loci, for oncogene annotation. |
| `ext`             | no | `1e7`   | bp to extend each amplicon by when searching for boundary SVs. |
| `mc.cores`        | no | `1`     | Cores for parallel processing across amplicons. |
| `verbose`         | no | `FALSE` | Print per-amplicon progress. |

### `plot_sv_linear()`

Only `sample`, `cnv_data` and `sv_data` are **required**; everything else is
optional.

**Required**

| Argument | Description |
|---|---|
| `sample`   | Sample identifier (matched in `cnv_data` / `sv_data`). |
| `cnv_data` | PURPLE CN segments: `sample`, `seqnames`, `start`, `end`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`. |
| `sv_data`  | PURPLE SVs/BEDPE: `chrom1`, `start1`, `chrom2`, `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`. |

**Optional — references & sample metadata**

| Argument | Default | Description |
|---|---|---|
| `wgd_data`       | `NULL` | WGD table (`Polyploidy` column); annotates WGD status in the title when supplied. |
| `karyotype`      | bundled hg38 | Ideogram bands (data.frame or `.rds` path). |
| `gene_coord`     | bundled hg38 oncogenes | Gene coordinates (data.frame or BED path). |
| `wgd_sample_col` | `NULL` | Sample-id column in `wgd_data` (`sample`, else `WGS_ID`). |

**Optional — which loci to plot**

| Argument | Default | Description |
|---|---|---|
| `chromosome`       | `NULL`  | Chromosome(s) to display. |
| `chromosome_range` | `NULL`  | Two-column start/end window(s), one row per `chromosome`. |
| `loci`             | `NULL`  | Explicit loci: data.frame (`chr`,`start`,`end`) or `"chr:start-end"` strings; overrides `chromosome`. |
| `events`           | `"amp"` | Auto-detect target(s): `"amp"`, `"gain"`, `"loh"`, `"homdel"`. |
| `cluster_gap`      | `5e6`   | Events more than this many bp apart become separate loci. |
| `flank_pct`        | `10`    | Extend each auto-detected window by ±X% of its width. |
| `min_amp_width`    | `1e5`   | Drop auto-detected loci with total event span below this (bp). |

**Optional — event thresholds**

| Argument | Default | Description |
|---|---|---|
| `min_cn_ratio`  | `3`   | Amplification: `copyNumber > min_cn_ratio × ploidy`. |
| `gain_ratio`    | `1.4` | Gain: `> gain_ratio × ploidy`, but not amplified. |
| `loh_thresh`    | `0.5` | LOH: `minorAlleleCopyNumber < loh_thresh`. |
| `homdel_thresh` | `0.5` | Homozygous deletion: `copyNumber < homdel_thresh`. |

**Optional — gene labels**

| Argument | Default | Description |
|---|---|---|
| `genes_to_highlight` | `NULL`  | Gene symbols to label (default oncogene panel if `NULL`). |
| `gene_label_angle`   | `NULL`  | Label rotation in degrees; auto-angles when crowded if `NULL`. |
| `repel_labels`       | `TRUE`  | De-collide labels with ggrepel. |
| `displayExon`        | `FALSE` | Draw exon models (requires `cds_gr`). |
| `cds_gr`             | `NULL`  | GRanges of CDS/exon ranges (`gene_name`), for `displayExon`. |

**Optional — layout & appearance**

| Argument | Default | Description |
|---|---|---|
| `cn_max`               | `NULL`  | Copy-number axis top (auto from data if `NULL`). |
| `gap_frac`             | `0.06`  | Gap between loci as a fraction of total width. |
| `offset_gene`          | `1.15`  | Gene-label height relative to max CN. |
| `ymax_highlight_ratio` | `1.08`  | Height of amp/homdel shading relative to max CN. |
| `karyotype_rel_size`   | `0.048` | Ideogram height relative to the CN axis. |
| `loh_position_ratio`   | `0.5`   | LOH/homdel bar position within the ideogram gap. |
| `highlight_amp`        | `TRUE`  | Shade amplified segments. |
| `highlight_hom_del`    | `TRUE`  | Shade homozygously deleted segments. |
| `plot_width_custom`    | `NULL`  | Override output width (inches). |
| `plot_height_custom`   | `NULL`  | Override output height (inches). |

**Optional — output**

| Argument | Default | Description |
|---|---|---|
| `outdir`  | `NULL`  | Directory to write the PDF; if `NULL`, only the plot object is returned. |
| `save`    | `TRUE`  | Write the PDF when `outdir` is supplied. |
| `verbose` | `FALSE` | Print progress / diagnostics. |

Returns the `ggplot` object (invisibly when a file is written); any written PDF
path is attached as `attr(p, "path")`.

## Status

Early development (v0.0.0.9000). API may change.

## License

MIT © Alhafidz Hamdan
