# Call translocation-bridge amplification (TBA) amplicons

Standalone caller for translocation-bridge (TB) amplification (Lee, Kim,
..., Park, *Nature* 2023; doi:10.1038/s41586-023-06057-w): a focal
amplicon whose boundary is defined by an **inter-chromosomal
translocation that is itself amplified** (the "boundary translocation"
carried within the amplicon). In this model an oncogene neighbourhood is
translocated in G1, creating a dicentric chromosome that forms a
chromosome bridge in mitosis; its breakage and repair amplify the
fragment (often as ecDNA), leaving an amplified translocation
demarcating the amplicon edge and co-amplifications on the partner
chromosome. Runs independently of the other mechanism callers.

## Usage

``` r
call_translocation_bridge_amp(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 5000,
  tb_edge_tol = 10000,
  centromeres = NULL,
  chrom_lengths = NULL,
  loh_max = 0.5,
  bridge_loh_min_frac = 0.15,
  nonbridge_loh_max_frac = 0.1,
  bridge_anchor_tol = 3e+06,
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

- tb_edge_tol:

  Integer; how close (bp) a translocation breakend must sit to an
  amplified amplicon edge to count as a boundary translocation (default
  `1e4`).

- centromeres, chrom_lengths:

  Optional
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of centromere spans (e.g.
  [`load_centromeres()`](https://alhafidzhamdan.github.io/EpiTracer/reference/load_centromeres.md))
  and named chromosome-length vector
  ([`load_chrom_lengths()`](https://alhafidzhamdan.github.io/EpiTracer/reference/load_chrom_lengths.md)).
  When both are supplied, each TBA amplicon is additionally scored for
  the paper's confirmatory footprint (see `tb_confident` in the return
  value); without them the confidence columns are `"NA"`.

- loh_max:

  Numeric; a copy-number segment is loss-of-heterozygosity (LOH) when
  its minor-allele copy number is below this (default `0.5`, i.e. ~0).

- bridge_loh_min_frac:

  Numeric; minimum fraction of a "bridge" arm that must be LOH, with an
  LOH block anchored near a landmark (telomere OR centromere), for that
  arm to count as a bridge arm (default `0.15`).

- nonbridge_loh_max_frac:

  Numeric; maximum LOH fraction allowed on the amplicon chromosome's
  opposite ("non-bridge") arm for it to count as spared (default `0.1`).

- bridge_anchor_tol:

  Integer; how close (bp) an LOH block must reach to a telomere or
  centromere to anchor the bridge-arm LOH (default `3e6`). The tolerance
  lets centromere-proximal LOH that starts a megabase or two past the
  centromere (a heterozygous peri-centromeric sliver in between) still
  anchor.

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints with `tba` (`"TRUE"`/`"FALSE"`),
`n_boundary_tra` (number of distinct boundary translocations) and
`tb_partner_chr` (comma-separated partner chromosome(s)). When
`centromeres` and `chrom_lengths` are supplied it also reports the
confirmatory footprint of Lee et al. (*Nature* 2023):
`tb_bridge_arm_loh` (bridge-arm LOH on the amplicon's own arm, to a
telomere or centromere), `tb_partner_arm_loh` (the same on a
translocated partner arm), `tb_nonbridge_spared` (the amplicon
chromosome's opposite arm retains heterozygosity), `tb_confident`
(`"TRUE"` when a bridge arm – amplicon and/or partner – shows LOH AND
the non-bridge arm is spared: the asymmetric pattern distinguishing TB
amplification from symmetric chromothripsis) and `tb_high_confidence`
(`"TRUE"` only for the dual-LOH pattern: LOH on BOTH the amplicon and
partner bridge arms with the non-bridge arm spared – the strongest
signature of a dicentric translocation bridge).

## See also

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md),
[`call_bfb()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.md),
[`call_chimeric_amplicon()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chimeric_amplicon.md)
