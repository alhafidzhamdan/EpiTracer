# EpiTracer

<!-- badges: start -->
<!-- badges: end -->

**EpiTracer** detects **episomal (breakage-independent) extrachromosomal DNA
(ecDNA)** from whole-genome sequencing (WGS) data and visualises the structural
rearrangements underlying focal amplifications.

The package currently provides two functions:

- **`call_episomal_ecdna()`** — an episomal ecDNA caller. For each ecDNA
  amplicon it locates the structural-variant breakpoints at the amplicon
  boundaries and flags amplicons whose structure is consistent with the
  *episome* model of formation: a circular amplicon bounded by a duplication
  (DUP) breakpoint, arising from an otherwise non-amplified chromosomal region,
  often leaving a deletion "excision scar" at the origin locus.
- **`plot_sv_linear()`** — a linear allele-specific copy-number / structural-
  variant "recon" plotter. It draws one or more loci side-by-side on a single
  concatenated x-axis with CN tracks, SV arcs, karyotype ideograms, LOH /
  homozygous-deletion bars, and gene labels (saved as PDF). Point it at a
  single focused locus (`chromosome` + `chromosome_range`), give explicit
  `loci`, or let it **auto-detect every amplified locus** in the sample — the
  structural variants interconnecting separate amplicons (multi-fragment / hub
  ecDNA junctions) are drawn as arcs spanning the loci.

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

# A single focused locus:
plot_sv_linear(
  sample     = "DO11441T1",
  cnv_data   = cnv_df, sv_data = sv_df, wgd_data = wgd_df,
  karyotype  = "chr_info_hg38.rds",           # ideogram bands
  gene_coord = "gene.coord_strand_name.bed",  # gene coordinates
  chromosome = "chr7",
  chromosome_range = matrix(c(52e6, 56e6), nrow = 1),
  outdir     = "plots"
)

# Every amplified locus, interconnected (multi-fragment ecDNA hub) — auto-detected:
plot_sv_linear(
  sample     = "DUMC12T1",
  cnv_data   = cnv_df, sv_data = sv_df, wgd_data = wgd_df,
  karyotype  = "chr_info_hg38.rds",
  gene_coord = "gene.coord_strand_name.bed",
  outdir     = "plots",                         # loci auto-detected; or loci = c("chr4:...", ...)
  flank_pct  = 10                              # extend each auto-detected window by +/-10%
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
