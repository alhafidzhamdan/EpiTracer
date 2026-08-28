# Nominate episomes by a VF-early circularisation DUP

Standalone caller that VF-stratifies each amplicon's junctions (as
[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md)
does) and flags a footprint-spanning boundary DUP in the EARLY tier –
the founder or the highest-VF cluster ("cluster 1"). The circularisation
DUP need not be the single highest-VF junction; an inverted- duplication
amplicon's fold-back can carry more copies, so requiring the boundary
DUP merely to be *early* rather than *maximal* is the correct test.

## Usage

``` r
call_founder_boundary(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 1e+05,
  max_k = 4L,
  span_frac = 0.8,
  early_strata = 2L,
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

- max_k:

  Maximum number of VF strata (default 4), as in
  [`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md).

- span_frac:

  A boundary DUP must span at least this fraction of the amplicon
  footprint to count (default 0.8).

- early_strata:

  Highest stratum index still considered "early" (default 2: the founder
  and the top cluster).

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints with `founder_class` (svclass of the max-VF
copy-gaining junction: DUP or an inversion), `founder_vf`,
`early_boundary_dup` (`"TRUE"`/`"FALSE"`), `early_dup_vf`,
`early_dup_stratum`, `n_strata`, `interlocus_tra` (`"TRUE"`/`"FALSE"`; a
translocation to another locus), `max_tra_vf` and `n_ancestral_del`
(deletions excluded from the founder call).

## Details

The founding event of an amplicon must be a copy-GAINING junction: a DUP
(circularisation / tandem) or an inversion (inverted duplication). A
DELETION removes sequence and therefore cannot found an amplification –
a high-VF deletion is an ANCESTRAL event (high-copy only because every
copy carries it), not a founder – so deletions are excluded from the
founding-mechanism call and counted in `n_ancestral_del`.
`founder_class` is thus DUP (simple excision) or an inversion (inverted
duplication, not simple excision).

Each amplicon is assessed on its OWN (intrachromosomal) junctions:
inter-locus translocations are excluded from the founding-mechanism call
and reported separately (`interlocus_tra`, `max_tra_vf`), because a TRA
joins this circle to another amplicon rather than forming it. Two
independently-episomal (DUP-founder) loci linked by a TRA are the
two-ecDNA-recombination (micronuclear chromothripsis) hypothesis.

## See also

[`plot_sv_reconstruction()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_reconstruction.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
