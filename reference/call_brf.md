# Call breakage-replication/fusion (BRF) amplicons

Standalone caller for the BRF annotation: an amplicon carrying *adjacent
parallel breakpoints*, the hallmark of breakage-replication/fusion
(Zhang, Mendez-Dorantes, Burns & Pellman, *Nat Genet* **58**, 88-99,
2026; doi:10.1038/s41588-025-02434-5). Runs independently of the
episomal call.

## Usage

``` r
call_brf(
  ecdna_gr = NULL,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  min_cn_ratio = 3,
  seed_gap = 1e+06,
  seed_min_width = 5000,
  max_dist = 20000,
  min_dist = 1,
  max_indep_p = 0.05,
  exclude_insertion_adjacency = TRUE,
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

- max_dist, min_dist, max_indep_p, exclude_insertion_adjacency:

  Passed to
  [`find_parallel_breakpoints()`](https://alhafidzhamdan.github.io/EpiTracer/reference/find_parallel_breakpoints.md);
  the defaults follow the source paper.

- mc.cores:

  Integer; number of cores for
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (default `1`). Values \> 1 are ignored on Windows.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints with `brf` (`"TRUE"`/`"FALSE"`),
`n_parallel_pairs` and `min_indep_p` (the most confident pair's
independence bound; `NA` when there is no pair).

## Details

A pair qualifies only when it clears all three of the source paper's
tests: the breakends are not part of a `+/-` insertion or overlap
adjacency, they share an orientation and sit within `max_dist` of one
another, and their estimated probability of independent origin — the
distance between them divided by the distance out to the nearest
opposite-orientation breakend — is at most `max_indep_p`. See
[`find_parallel_breakpoints()`](https://alhafidzhamdan.github.io/EpiTracer/reference/find_parallel_breakpoints.md).

The independence test is what keeps this specific. Without it the caller
fires on essentially any junction-dense amplicon, because
same-orientation breakends land within 20 kb of each other by chance
once a footprint carries tens of junctions;
`validation/simulate_trajectories.R` shows the failure mode on simulated
chromothriptic amplicons that contain no BRF event by construction.

## See also

[`find_parallel_breakpoints()`](https://alhafidzhamdan.github.io/EpiTracer/reference/find_parallel_breakpoints.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
