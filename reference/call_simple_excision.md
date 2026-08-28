# Call episomal extrachromosomal DNA from WGS structural variants

Detects ecDNA amplicons whose structure is consistent with the *episome*
(breakage-independent) model of formation: a circular amplicon bounded
by a duplication breakpoint, arising from an otherwise non-amplified
chromosomal region, and often leaving a deletion "excision scar" at the
origin locus.

## Usage

``` r
call_simple_excision(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  mc.cores = 1L,
  verbose = FALSE,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 1e+05,
  flank_baseline = c("chromosome", "ploidy"),
  gain_ratio = 1.4,
  min_flank_width = 2000,
  bridge_gap = seed_gap,
  founder_jcn = 30,
  centromeres = NULL,
  mh_min_homology = 2,
  hr_min_homology = 14
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

- flank_baseline:

  How the "flanks not gained" test is calibrated: `"chromosome"`
  (default) compares the amplicon flanks to the local per-chromosome
  baseline (the width-weighted median copy number of that chromosome's
  non-focally-amplified segments), `"ploidy"` compares them to the
  global sample ploidy. Use `"chromosome"` for focal episomes on a
  polysomic chromosome (e.g. EGFR on the gained chr7 in glioblastoma),
  which `"ploidy"` misses because the flanks sit above tumour ploidy.

- gain_ratio:

  Numeric; a flank is "gained" when its copy number is at or above
  `gain_ratio` times the baseline (default `1.4`).

- min_flank_width:

  Integer; ignore copy-number segments narrower than this (in bp) when
  reading the flanking copy number, so a sub-kb "segmentation shoulder"
  emitted at a sharp amplicon edge is not mistaken for the true flank
  (default `2000`). The nearest segment at least this wide is tested; if
  none qualifies, the abutting segment is used as a fallback.

- bridge_gap:

  Numeric; the maximum non-amplified gap (in bp) an amplicon-connecting
  junction may span before the amplicon is flagged as `bridging` and
  rejected (default: `seed_gap`, i.e. the same distance that separates
  distinct amplicon seeds). A junction whose both breakends sit at
  amplicon-level copy number but that spans a wider gap fuses two
  separate amplicons rather than bounding one episome; a genuine
  internal deletion leaves only a sub-`bridge_gap` gap.

- founder_jcn:

  Numeric; the junction-copy-number above which an SV counts as a
  clonal, high-copy "founder" junction for the micronucleus test
  (default `30`). When a founder-level inversion begins inside a
  founder-level duplication's span (interleaved DUP+INV), the amplicon
  is a micronucleus / chromothripsis product and is called non-episomal.
  Low-JCN (subclonal) interleaving – later internal rearrangement of a
  genuine episome – is ignored.

- centromeres:

  Optional
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of centromere spans (one per chromosome, e.g. from
  [`load_centromeres()`](https://alhafidzhamdan.github.io/EpiTracer/reference/load_centromeres.md)).
  When supplied, an amplicon whose amplified edge terminates *within* a
  centromere (rather than traversing across it) is flagged
  `flag_chromosomal_bridge` and excluded from the episomal call as a
  dicentric / chromosomal-bridge event. `NULL` (default) skips the
  check.

- mh_min_homology, hr_min_homology:

  Numeric; breakpoint-homology (bp) thresholds used to classify the
  circularisation (boundary-DUP) junction into an inferred DSB-repair
  pathway, following Eugen-Olsen et al. (Nucleic Acids Res 2025;
  doi:10.1093/nar/gkaf122). A junction with homology `< mh_min_homology`
  is called `"NHEJ"` (near-blunt), `>= mh_min_homology` and
  `< hr_min_homology` is `"MMEJ"` (short microhomology / alternative
  end-joining), and `>= hr_min_homology` is `"HR"` (long homology /
  homologous recombination). Defaults `2` and `14` are the human
  microhomology and HR minimum-homology lengths reported in that review.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
combining the annotated breakpoints of all amplicons, with
per-breakpoint classification columns including
`duplication_at_boundary`, `duplication_at_boundary_has_highest_VF`,
`episome_region`, `deletion_flanking_boundary`, `episomal`,
`has_excision_scar`, and the diagnostic flags `flag_internal_sv_high_vf`
(an internal SV out-VFs the boundary DUP; informational), the
disqualifying `flag_internal_inversion` (a fold-back inversion does so –
break-fusion- bridge), and `flag_bridging_amplicon` (a junction fuses
two separate amplicons) – all character `"TRUE"`/`"FALSE"`. An amplicon
carrying either disqualifying flag is set `episomal = "FALSE"`. For an
amplicon with a boundary DUP it also reports the
circularisation-junction microhomology as `boundary_homology` (numeric
bp) and `junction_homology_class` (`"NHEJ"`/`"MMEJ"`/`"HR"`; `NA` when
no boundary DUP is found), the inferred DSB-repair pathway that sealed
the circle.

The other amplicon-formation mechanisms are computed by dedicated
standalone callers, not here:
[`call_brf()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_brf.md)
(breakage-replication/fusion / adjacent parallel breakpoints),
[`call_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_micronucleation.md)
(high-VF non-homologous translocation),
[`call_bfb()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.md)
(breakage-fusion-bridge) and
[`call_translocation_bridge_amp()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_translocation_bridge_amp.md)
(translocation-bridge amplification). Join their output to this one by
`WGS_ID` + `ID` to assemble a combined mechanism table.

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

episomal <- call_simple_excision(
  ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
  ext = 1e7, mc.cores = 4
)
} # }
```
