# Detect focal-amplicon seeds from copy number alone

Derives candidate focal-amplicon regions ("seeds") directly from
allele-specific copy-number segments, so
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
can run **without** an AmpliconArchitect amplicon catalogue. Per sample,
segments with `copyNumber > min_cn_ratio * ploidy` are merged across
gaps up to `gap` bp; adjacent regions of SIMILAR copy number are then
merged across larger gaps (up to `merge_gap`), so a single amplicon
split by internal deletions is not mistaken for several separate
amplicons (which would let a duplication that is internal to the true
amplicon pass as a boundary DUP). Regions narrower than `min_width` are
dropped; each surviving region is labelled with an `ID` and `WGS_ID`,
matching the `ecdna_gr` contract.

## Usage

``` r
detect_amplicon_seeds(
  cnv_gr,
  min_cn_ratio = 3,
  gap = 2e+06,
  min_width = 5000,
  merge_gap = 3e+06,
  cn_ratio = 0.5,
  breakpoints = NULL,
  link_tol = 10000
)
```

## Arguments

- cnv_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of allele-specific copy-number segments with metadata columns
  `sample`, `copyNumber`, `ploidy`.

- min_cn_ratio:

  Numeric; a segment is amplified where
  `copyNumber > min_cn_ratio * ploidy` (default `3`).

- gap:

  Integer; merge amplified segments separated by at most this many bp
  into one seed (default `1e6`).

- min_width:

  Integer; drop seeds narrower than this (default `5e3`).

- merge_gap:

  Integer; after the initial `gap` merge, also merge two adjacent
  same-chromosome seeds separated by at most this many bp when their
  representative copy numbers are similar (see `cn_ratio`). Models one
  amplicon broken only by internal deletions (default `3e6`). Set to
  `gap` to disable copy-number-aware merging.

- cn_ratio:

  Numeric in `(0, 1]`; two adjacent seeds count as "similar copy number"
  (and are eligible to merge across up to `merge_gap`) when
  `min(cn1, cn2) / max(cn1, cn2) >= cn_ratio` (default `0.5`, i.e.
  within 2x).

- breakpoints:

  Optional per-breakend table (a
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  or data.frame with `seqnames`/`chr`, `start`/`pos`, `WGS_ID`/`sample`,
  `event`, `PURPLE_CN`). When supplied the seeds become **SV-aware**:
  two seeds on the same chromosome joined by a structural variant whose
  *both* breakends sit in amplified copy number
  (`PURPLE_CN > min_cn_ratio * ploidy`, the same bar a segment must
  clear to be a seed) are linked into one amplicon, even across a large
  gap (e.g. a centromere) – so a single rearranged amplicon that
  traverses a centromere is not split. A junction that dips into
  non-amplified sequence on either side does not link. In addition, when
  one breakend of such a junction falls in a seed and the other is
  amplified but sits in a block too narrow to have seeded
  (`< min_width`), the seed is **extended** to reach that far breakend –
  so a boundary DUP that spans a small distal amplified block is
  recognised rather than lost. Each returned seed also carries `cn_only`
  (`TRUE` when no breakend falls inside it, i.e. it is called by copy
  number alone with no supporting junction).

- link_tol:

  Integer; tolerance (bp) when mapping a breakend to a seed for SV-aware
  linking and the `cn_only` test (default `1e4`). A boundary-defining
  breakend typically sits *at* the amplified-segment edge, so the
  reduced seed may begin a base or two after it; strict containment
  would then miss exactly the junction that anchors the boundary.
  Matching within `link_tol` (the same window used to read a breakend's
  `PURPLE_CN`) lets an edge-anchored breakend map to its seed.

## Value

A
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
of amplicon seeds with metadata columns `ID`, `WGS_ID`, and (when
`breakpoints` is supplied) `cn_only` (empty if no amplification is
present).

## See also

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)

## Examples

``` r
seeds <- detect_amplicon_seeds(ex_caller_inputs$cnv_gr)
#> Loading required namespace: GenomeInfoDb
seeds
#> GRanges object with 1 range and 3 metadata columns:
#>       seqnames            ranges strand |      WGS_ID             ID   cn_only
#>          <Rle>         <IRanges>  <Rle> | <character>    <character> <logical>
#>   [1]     chr7 55000000-55500000      * |   EXAMPLE01 EXAMPLE01_amp1      <NA>
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
