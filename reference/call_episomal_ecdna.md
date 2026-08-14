# Call episomal extrachromosomal DNA from WGS structural variants

Detects ecDNA amplicons whose structure is consistent with the *episome*
(breakage-independent) model of formation: a circular amplicon bounded
by a duplication breakpoint, arising from an otherwise non-amplified
chromosomal region, and often leaving a deletion "excision scar" at the
origin locus.

## Usage

``` r
call_episomal_ecdna(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  mc.cores = 1L,
  verbose = FALSE,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 1e+05
)
```

## Arguments

- ecdna_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of ecDNA amplicon regions with metadata columns `ID` (unique amplicon
  identifier) and `WGS_ID` (sample identifier) — typically the
  AmpliconArchitect amplicon catalogue. If `NULL` (the default),
  focal-amplicon seeds are detected from `cnv_gr` with
  [`detect_amplicon_seeds()`](https://alhafidzhamdan.github.io/EpiTracer/reference/detect_amplicon_seeds.md),
  so EpiTracer can run without AmpliconArchitect.

- breakpoints_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of PURPLE (HMF pipeline) SV breakpoints with metadata columns
  `WGS_ID`, `event`, `svclass` (e.g. "DUP", "DEL"), `PURPLE_AF`,
  `PURPLE_JCN`, `VF`, `PURPLE_CN`, `insLen`, `HOMLEN`.

- cnv_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of PURPLE (HMF pipeline) allele-specific copy-number segments with
  metadata columns `sample`, `copyNumber`, `ploidy`,
  `majorAlleleCopyNumber`, `minorAlleleCopyNumber`.

- cancer_genes_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of cancer gene loci with a `gene` (gene symbol) metadata column, used
  to annotate breakpoints with the overlapping oncogene.

- ext:

  Integer; number of base pairs to extend each amplicon by when
  searching for boundary breakpoints (default `1e7`).

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

- verbose:

  Logical; print per-amplicon progress (default `FALSE`).

- min_cn_ratio, seed_gap, seed_min_width:

  Passed to
  [`detect_amplicon_seeds()`](https://alhafidzhamdan.github.io/EpiTracer/reference/detect_amplicon_seeds.md)
  when `ecdna_gr` is `NULL` (copy-number amplicon threshold
  `copyNumber > min_cn_ratio * ploidy`, gap to merge across, and minimum
  seed width). Ignored when `ecdna_gr` is supplied.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
combining the annotated breakpoints of all amplicons, with
per-breakpoint classification columns including
`duplication_at_boundary`, `duplication_at_boundary_has_highest_VF`,
`episome_region`, `deletion_flanking_boundary`, `episomal`, and
`has_excision_scar` (all character `"TRUE"`/`"FALSE"`).

## Details

Each amplicon called by AmpliconArchitect (supplied in `ecdna_gr`, a
compulsory input) is processed independently (optionally in parallel).
Within each amplicon the function annotates the flanking structural
variant breakpoints with oncogene and allele-specific copy-number
context, then applies the heuristic described in
[`classify_amplicon_episomal()`](https://alhafidzhamdan.github.io/EpiTracer/reference/classify_amplicon_episomal.md).

## See also

[`classify_amplicon_episomal()`](https://alhafidzhamdan.github.io/EpiTracer/reference/classify_amplicon_episomal.md),
[`plot_sv_linear()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_linear.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ecdna_gr        <- readRDS("ecDNA_amplicon_regions.rds")
breakpoints_gr  <- readRDS("SV_catalogue.rds")
cnv_gr          <- readRDS("CN_segments.rds")
cancer_genes_gr <- readRDS("cancer_genes.rds")

episomal <- call_episomal_ecdna(
  ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
  ext = 1e7, mc.cores = 4
)
} # }
```
