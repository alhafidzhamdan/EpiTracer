# EpiTracer

<!-- badges: start -->
<!-- badges: end -->

**EpiTracer** detects **episomal (breakage-independent) extrachromosomal DNA
(ecDNA)** from whole-genome sequencing (WGS) data and visualises the structural
rearrangements underlying focal amplifications.

The package currently provides three functions:

- **`call_episomal_ecdna()`** — an episomal ecDNA caller. For each ecDNA
  amplicon it locates the structural-variant breakpoints at the amplicon
  boundaries and flags amplicons whose structure is consistent with the
  *episome* model of formation: a circular amplicon bounded by a duplication
  (DUP) breakpoint, arising from an otherwise non-amplified chromosomal region,
  often leaving a deletion "excision scar" at the origin locus.
- **`plot_amplicon_recon()`** — a multi-locus amplicon plot on a single
  concatenated axis: focuses each amplified locus, lays them side-by-side, and
  draws the structural variants interconnecting separate amplicons (multi-fragment
  / hub ecDNA junctions), with per-locus ideograms.
- **`plot_sv_linear()`** — a linear allele-specific copy-number / structural-
  variant "recon" plotter for one or more loci, drawing CN tracks, SV arcs, a
  karyotype ideogram, LOH / homozygous-deletion bars, and gene labels to PDF.

## Installation

EpiTracer depends on [`gUtils`](https://github.com/mskilab/gUtils) (GitHub) and
several Bioconductor/CRAN packages. Install with:

```r
# install.packages("remotes")
remotes::install_github("alhafidzhamdan/EpiTracer")
```

If the Bioconductor dependencies are not already present:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "regioneR"))
remotes::install_github("mskilab/gUtils")
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

plot_sv_linear(
  sample     = "DO12742T1",
  chromosome = c("chr7", "chr12"),
  karyotype  = "chr_info_hg38.rds",           # ideogram bands
  gene_coord = "gene.coord_strand_name.bed",  # gene coordinates
  wgd_data   = wgd_df,
  cnv_data   = cnv_df,
  sv_data    = sv_df,
  outdir     = "plots"
)
```

## Input format

`call_episomal_ecdna()` expects PURPLE/HMF-style metadata columns:

| Object            | Required metadata |
|-------------------|-------------------|
| `ecdna_gr`        | `ID`, `WGS_ID` |
| `breakpoints_gr`  | `WGS_ID`, `event`, `svclass`, `PURPLE_AF`, `PURPLE_JCN`, `VF`, `PURPLE_CN`, `insLen`, `HOMLEN` |
| `cnv_gr`          | `sample`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber` |
| `cancer_genes_gr` | any (used for overlap annotation) |

## Status

Early development (v0.0.0.9000). API may change.

## License

MIT © Alhafidz Hamdan
