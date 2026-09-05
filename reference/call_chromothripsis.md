# Call chromothripsis within (episomal) amplicons

Standalone caller that scores each amplicon's INTERNAL structural
variants for chromothripsis, using the ShatterSeek hallmarks
(Cortes-Ciriano et al., Nat Genet 2020) restricted to the amplified
footprint. It is designed to be run on amplicons already called episomal
(pass them as `ecdna_gr`) to flag the subset that have since shattered –
the ecDNA -\> micronucleus -\> chromothripsis route – but works on any
amplicon catalogue (or `NULL` to auto-detect seeds).

## Usage

``` r
call_chromothripsis(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 5000,
  min_sv = 6L,
  min_oscillations = 3L,
  join_p = 0.05,
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

- min_sv:

  Minimum number of distinct internal SV events for the prevalence
  hallmark (default 6).

- min_oscillations:

  Minimum copy-number direction changes (turning points) across the
  footprint segments for the oscillation hallmark (default 3).

- join_p:

  Significance threshold for the fragment-join randomness test; the
  junction-orientation distribution must NOT differ from uniform at this
  level (default 0.05).

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints with `chromothripsis` (`"TRUE"`/`"FALSE"`),
`chromothripsis_conf` (`"high"`/`"low"`/`"none"`), `n_internal_sv`,
`n_intrachr_sv`, `sv_type_pval`, `cn_oscillations` and
`loh_interspersed`.

## Details

A footprint is called chromothriptic when it carries at least `min_sv`
distinct internal SV events, its four intrachromosomal junction
orientations are close to equally represented (chi-squared
goodness-of-fit p \>= `join_p`, i.e. random fragment joins), and its
rounded copy-number profile changes direction at least
`min_oscillations` times across the footprint (oscillating copy number).
All three give `chromothripsis_conf = "high"`. Random fragment joins are
required for any positive call (they separate chromothripsis from
orientation-biased mechanisms such as BFB); prevalence with random joins
but weak oscillation gives `"low"`. A clean simple episome (one boundary
junction, few internal SVs) fails the prevalence test and is called
`"FALSE"`.

## See also

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md),
[`call_chimeric_amplicon()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chimeric_amplicon.md),
[`call_bfb()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.md)
