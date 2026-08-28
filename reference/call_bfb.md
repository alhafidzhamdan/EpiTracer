# Call breakage-fusion-bridge (BFB) amplicons

Standalone caller for the classical BFB signature (McClintock; genomic
signature after Bignell et al. 2007): (i) fold-back inversions, (ii)
intrachromosomal (no amplified translocation partner), (iii) terminal
amplification with a distal deletion running contiguously to the
(absent) telomere, and (iv) a copy-number staircase (fold-backs spread
across several stepped levels – the feature that distinguishes iterative
BFB from a single BRF event or a focal ecDNA spike).

## Usage

``` r
call_bfb(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 1e+05,
  centromeres = NULL,
  chrom_lengths = NULL,
  bfb_min_del_width = 1e+06,
  bfb_min_del_frac = 0.7,
  bfb_loss_max = 1.5,
  bfb_min_levels = 3,
  bfb_min_spread = 0.3,
  mc.cores = 1
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

- min_cn_ratio, seed_gap, seed_min_width:

  Passed to
  [`detect_amplicon_seeds()`](https://alhafidzhamdan.github.io/EpiTracer/reference/detect_amplicon_seeds.md)
  when `ecdna_gr` is `NULL` (copy-number amplicon threshold
  `copyNumber > min_cn_ratio * ploidy`, gap to merge across, and minimum
  seed width). Ignored when `ecdna_gr` is supplied.

- centromeres, chrom_lengths, bfb_min_del_width, bfb_min_del_frac,
  bfb_loss_max, bfb_min_levels, bfb_min_spread:

  See
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md);
  `centromeres` and `chrom_lengths` are both required (the BFB
  annotation is disabled without them).

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints with `bfb` (`"TRUE"`/`"FALSE"`), `bfb_anchor`
and `n_foldbacks`.

## See also

[`load_centromeres()`](https://alhafidzhamdan.github.io/EpiTracer/reference/load_centromeres.md),
[`load_chrom_lengths()`](https://alhafidzhamdan.github.io/EpiTracer/reference/load_chrom_lengths.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
