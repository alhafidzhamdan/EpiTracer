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
  min_occupancy = 0.05
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
boundary-DUP span; the occupancy-guard metric), the bridging-junction
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
