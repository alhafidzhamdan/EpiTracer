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

EpiTracer is **not on CRAN** — it depends on Bioconductor packages
(`GenomicRanges`, `regioneR`) and on [`gUtils`](https://github.com/mskilab/gUtils)
(GitHub-only). Install it with any tool that resolves Bioconductor + GitHub
dependencies (the `Remotes:` field pulls in `gUtils` automatically):

**BiocManager** (recommended — resolves Bioconductor, GitHub and `Remotes:`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("alhafidzhamdan/EpiTracer")
```

**pak** (also resolves everything in one call):

```r
# install.packages("pak")
pak::pak("alhafidzhamdan/EpiTracer")
```

**remotes / devtools** (enable Bioconductor repos first so the Bioc deps resolve):

```r
install.packages(c("remotes", "BiocManager"))
options(repos = BiocManager::repositories())
remotes::install_github("alhafidzhamdan/EpiTracer")
```

**r-universe** — a CRAN-like experience without CRAN. If the package is built on
a personal [r-universe](https://r-universe.dev), users install with plain
`install.packages()`:

```r
install.packages("EpiTracer",
  repos = c("https://alhafidzhamdan.r-universe.dev",
            "https://bioc.r-universe.dev", getOption("repos")))
```

**Local source tarball** (dependencies must be pre-installed):

```sh
R CMD build EpiTracer
```
```r
install.packages("EpiTracer_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

> The GitHub-based options require the repository to be pushed to
> `github.com/alhafidzhamdan/EpiTracer` first.

### Offline / GitHub-only environments

On a firewalled host where CRAN/Bioconductor are unreachable but GitHub is
whitelisted, how you install depends on whether the dependencies are already
present (check e.g. `requireNamespace("GenomicRanges")`). The bundled hg38
references mean no data download is needed at run time.

**Dependencies already installed** (common on managed / HPC R installs) — pull
just EpiTracer and `gUtils` from GitHub, without touching CRAN/Bioconductor:

```r
remotes::install_github("mskilab/gUtils",           dependencies = FALSE, upgrade = "never")
remotes::install_github("alhafidzhamdan/EpiTracer", dependencies = FALSE, upgrade = "never")
```

**Fully air-gapped** — stage everything on a connected machine, transfer, then
install locally. Either download the whole dependency tree:

```r
# on a machine with internet:
pak::pkg_download("alhafidzhamdan/EpiTracer", dest_dir = "epitracer_pkgs", dependencies = TRUE)
# copy epitracer_pkgs/ across, then on the offline host:
install.packages(dir("epitracer_pkgs", full.names = TRUE), repos = NULL, type = "source")
```

…or install from a clone (dependencies must already be present):

```sh
git clone https://github.com/alhafidzhamdan/EpiTracer.git
R CMD INSTALL EpiTracer
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
