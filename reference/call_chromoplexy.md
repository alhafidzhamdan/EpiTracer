# Call chromoplexy (closed balanced rearrangement chains)

Sample-level caller for chromoplexy, the punctuated formation of a
closed chain of balanced, low-copy rearrangements whose breakends meet
in *deletion bridges* (Baca et al., *Cell* 2013,
doi:10.1016/j.cell.2013.03.021; graph criteria after the gGnome
`chromoplexy()` caller, Hadi et al., *Cell* 2020,
doi:10.1016/j.cell.2020.08.006). Each structural variant is
reconstructed as a junction between two breakends; balanced
(ploidy-scaled low-copy) junctions are kept; breakends within `max_dist`
collapse into loci (a locus joining two distinct junctions is a deletion
bridge); and a chromoplexy is called for every connected component of
the locus/junction graph that forms a **clean simple cycle** – every
locus of degree two, with at least `min_pairs` bridges closing the ring.
Runs independently of the amplicon-formation callers.

## Usage

``` r
call_chromoplexy(
  breakpoints_gr,
  cnv_gr,
  min_pairs = 3L,
  min_span = 1e+07,
  max_dist = 10000,
  max_cn = 3,
  footprint_width = 1e+06,
  max_small = 50000,
  min_chromosomes = 3L,
  mc.cores = 1
)
```

## Arguments

- breakpoints_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of SV breakends (two per `event`), as used by
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md):
  metadata `WGS_ID`, `event`, `svclass`, `PURPLE_CN`.

- cnv_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of allele-specific copy-number segments with metadata `sample`,
  `copyNumber`, `ploidy` (the per-chromosome baseline used to
  ploidy-scale the low-copy filter).

- min_pairs:

  Integer; minimum number of deletion bridges (equivalently the cycle
  length) required to call a chromoplexy (default `3L`, per Baca et
  al.).

- min_span:

  Numeric; minimum span (bp) for an intrachromosomal junction to count
  as a long-range chain edge (default `1e7`); interchromosomal junctions
  always qualify.

- max_dist:

  Integer; how close (bp) two breakends must sit to collapse into the
  same locus / deletion bridge (default `1e4`).

- max_cn:

  Numeric; a junction is balanced (kept) when its copy number is at most
  `max_cn * ploidy / 2` (default `3`, i.e. \<=3 copies at diploid
  baseline, scaled up on polysomic chromosomes).

- footprint_width:

  Numeric; padding (bp) applied when reporting locus footprints (default
  `1e6`).

- max_small:

  Numeric; a DUP/DEL junction spanning at most this (bp) is dropped as
  local noise before chain-building (default `5e4`).

- min_chromosomes:

  Integer; minimum number of DISTINCT chromosomes the cycle must span
  (default `3L`). Chromoplexy is characteristically interchromosomal
  (Baca et al.); this rejects closed rings that merely link loci on one
  or two chromosomes. Set to `2L` (or `1L`) to relax.

- mc.cores:

  Integer; samples processed in parallel (default `1`).

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with one row per called chromoplasy cycle: `WGS_ID`, `chromoplexy_id`,
`topology` (always `"cycle"` for called events), `n_junctions`,
`n_bridges`, `n_chromosomes`, `chromosomes`, `max_cn`, `min_cn`,
`min_span`, `n_interchrom`, `num_other` (intruding foreign breakends),
`frac_cp` (cleanliness `n_junctions / (n_junctions + num_other)`),
`footprint` and `events` (comma-separated `event` ids in the cycle, for
joining back onto amplicons). Empty (0 rows) when no chromoplexy is
found.

## See also

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md),
[`call_translocation_bridge_amp()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_translocation_bridge_amp.md)
