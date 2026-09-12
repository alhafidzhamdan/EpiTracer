# Call episomal ecDNA from AmpliconSuite output

Classifies a single AmpliconClassifier **ecDNA feature** as episomal
(simple-excision circle) or not, from its AmpliconArchitect graph and
its ecDNA interval(s). The amplicon extent is taken from the
AmpliconClassifier ecDNA interval BED (AA's own amplicon call) – no
copy-number threshold is applied here. The feature is **episomal** when
a BOUNDARY DUP exists: an intra-chromosomal DUP whose ends bracket the
outer boundary of the ecDNA interval on its chromosome
(`lo <= interval start` and `hi >= interval end`). The boundary DUP
reported is the one with the highest read-pair support among spanning
DUPs. A multi-chromosome ecDNA is **chimeric** episomal when every
chromosome is bracketed by a boundary DUP and joined by \>= 1
translocation.

## Usage

``` r
call_episomal_from_ampliconsuite(
  graph_file,
  intervals,
  tol = 0L,
  bridge_gap = 1e+06,
  min_occupancy = 0.05,
  trim_to_amplified = TRUE,
  min_seg_cn_frac = 0.1,
  structural_qc = TRUE,
  qc_amp_frac = 0.5,
  qc_src_frac = 0.5,
  qc_jcn_frac = 0.1,
  qc_edge_tol = 10000,
  ploidy = NULL,
  cnv_file = NULL,
  min_cn_ratio = 3,
  default_ploidy = 2
)
```

## Arguments

- graph_file:

  Path to the amplicon's AmpliconArchitect `*_graph.txt`.

- intervals:

  The ecDNA feature's intervals: a path to an AmpliconClassifier
  `*_ecDNA_<k>_intervals.bed`, or a
  `data.frame`/`data.table`/GRanges-like object with
  `chr`/`start`/`end`. AmpliconClassifier's exclusive-end convention is
  handled internally.

- tol:

  Integer; coordinate tolerance (bp) for the boundary-DUP bracket test
  (default `0`; graph and interval coordinates come from the same graph
  so no tolerance is normally needed).

- min_occupancy:

  Numeric in `(0, 1]`; minimum fraction of the boundary-DUP bracketed
  span that must be covered by the ecDNA intervals for an episomal call
  (default `0.05`, i.e. reject when \>95% of the bracketed span is
  non-amplified). A simple-excision circle is near-contiguous, so a DUP
  spanning a largely non-amplified region (e.g. a lone DUP across a
  diploid chromosome arm) is rejected as not a self-contained circle.
  The default was set by a PCAWG-vs-PURPLE concordance sweep (optimal at
  2-5%; 5% additionally removes scattered low-occupancy secondary
  features).

- trim_to_amplified:

  Logical; trim the AmpliconClassifier ecDNA interval to its amplified
  extent before the boundary-DUP test (default `TRUE`).
  AmpliconClassifier occasionally pads an interval with flanking diploid
  sequence, pushing the outer extent past the boundary DUP and
  mis-calling a clean episome as complex; trimming uses the AA graph's
  own segment copy numbers to drop leading/trailing near-diploid
  segments (edges only – internal low-CN segments are kept). It can only
  shrink the extent, so it may rescue a false-complex call but never
  breaks a true-episomal one.

- min_seg_cn_frac:

  Numeric in `(0, 1]`; a segment counts as amplified for
  `trim_to_amplified` when its copy number is at least this fraction of
  the amplicon's peak segment copy number (default `0.10`). Relative to
  the amplicon's own amplification, not a global ploidy cutoff.

- structural_qc:

  Logical (default `TRUE`); apply structural quality control to a
  bracketed amplicon before calling it episomal, using the AA graph's
  own segment copy numbers, source (open-end) edges and junction copy
  numbers. Rejects (→ `status = "complex"`) an amplicon that is not
  amplified at both boundary-DUP breakends, is OPEN at the boundary (a
  large source edge at a breakend), is closed by a boundary DUP whose
  copy number is far below the circle CN, or has an INV/TRA at a
  boundary out-copying the DUP. Ported from the PURPLE-native
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md).
  Set `FALSE` for the older bracket-and-occupancy behaviour.

- qc_amp_frac, qc_src_frac, qc_jcn_frac:

  Numeric in `(0, 1]`; the structural-QC thresholds, each relative to
  the width-weighted median copy number across the boundary-DUP span
  (the "circle-body" CN). A boundary breakend must be amplified to
  `>= qc_amp_frac` × circle CN (default `0.5`); a boundary source edge
  `>= qc_src_frac` × circle CN marks an open boundary (default `0.5`);
  the boundary DUP's junction copy number must be `>= qc_jcn_frac` ×
  circle CN (default `0.1`; a low floor that rejects only a boundary DUP
  explaining a tiny fraction of the amplification, while tolerating
  internally re-amplified episomes whose resident locus is boosted above
  the circularising junction's copy number — e.g. DO11854 EGFR).

- qc_edge_tol:

  Integer; window (bp) around a boundary-DUP breakend within which a
  source edge or INV/TRA breakend counts as "at the boundary" (default
  `1e4`).

- ploidy, cnv_file:

  Tumour ploidy for the amplification / tandem-duplication gate. Supply
  `ploidy` directly, or `cnv_file` – a path to a genome-wide copy-number
  BED (e.g. the AA-input `*_CNV_CALLS.bed`) whose width-weighted median
  CN is used. When ploidy is known and structural QC is on, a bracketed
  amplicon whose circle-body CN does not exceed `min_cn_ratio` x ploidy
  is downgraded (`qc_fail = "below_ploidy"`): at low copy number a
  boundary DUP + copy-number step is indistinguishable from an in-place
  tandem duplication by unphased short-read data, so a genuine ecDNA
  must be amplified well above ploidy. Both `NULL` (the default) skips
  the gate.

- min_cn_ratio:

  Numeric; the circle-body CN must exceed this multiple of ploidy to
  pass the amplification gate (default `3`, matching the PURPLE-native
  [`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)).

- default_ploidy:

  Numeric; the ploidy assumed when neither `ploidy` nor `cnv_file`
  resolves a value (default `2`, diploid), so the amplification gate
  still applies an absolute floor of `min_cn_ratio` x `default_ploidy`
  (= 6 at the defaults) rather than being skipped. Set `NULL` to skip
  the gate when ploidy is unknown. When the fallback is used,
  `ploidy_assumed` is `TRUE`.

  A **bridged/ambiguous downgrade** applies even when a bracketing
  boundary DUP exists: if the amplicon's highest read-pair-support
  junction (a breakend inside the ecDNA extent) has its other breakend
  on the **same chromosome**, \> `bridge_gap` beyond the extent, and in
  an amplified segment, the dominant rearrangement co-mobilises a
  co-amplified region, so `status` is set to `"bridged"` and `episomal`
  to `FALSE`. Cross-chromosome TRAs are the chimeric-episome signature
  and are NOT downgraded.

## Value

A one-row
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with `episomal` (logical; `TRUE` only for a clean episome), `status`
(`"episomal"`/`"bridged"`/`"complex"`), `chimeric` (logical), `n_chr`
(chromosomes spanned by the ecDNA), the boundary-DUP record `bd_chr`,
`bd_lo`, `bd_hi`, `bd_pe_support`, `bd_jcn`, `bd_homlen` (`NA` when no
boundary DUP), `bd_occupancy` (ecDNA-interval coverage of the
boundary-DUP span; the occupancy-guard metric), the structural-QC record
`bd_circle_cn` (width-weighted median CN of the boundary-DUP span),
`bd_jcn_frac` (boundary-DUP jcn ÷ circle CN), `bd_src_frac` (boundary
source-edge CN ÷ circle CN), `ploidy` (the ploidy used for the gate,
resolved or the `default_ploidy` fallback), `ploidy_assumed` (logical;
`TRUE` when the fallback was used), `cn_over_ploidy` (circle CN ÷
ploidy), `qc_pass` (logical) and `qc_fail` (`;`-separated failed checks
incl. `below_ploidy`, `NA` if passed or QC off), the bridging-junction
record `bridge_pe_support`, `bridge_partner` (`chr:pos` of the
out-of-amplicon breakend; `NA` unless downgraded), and `n_junctions`.

## See also

[`read_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_ampliconsuite.md),
[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)

## Examples

``` r
if (FALSE) { # \dontrun{
call_episomal_from_ampliconsuite(
  graph_file = "SAMPLE_amplicon1_graph.txt",
  intervals  = "SAMPLE_amplicon1_ecDNA_1_intervals.bed")
} # }
```
