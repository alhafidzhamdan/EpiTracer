[← EpiTracer](../README.md)

# `call_episomal_from_ampliconsuite()` — An episomal ecDNA caller for AmpliconSuite output

> Applies EpiTracer's **episome model** (a circle bounded by a duplication
> junction) directly to [AmpliconSuite](https://github.com/AmpliconSuite) output,
> so community-standard ecDNA calls — e.g. from
> [AmpliconRepository](https://ampliconrepository.org) — can be classified as
> **episomal** (simple excision) vs **bridged** / **complex**. It is a
> self-contained module; it does not alter the PURPLE-native
> [`call_simple_excision()`](call_simple_excision.md).

A tool to **classify one AmpliconClassifier ecDNA feature** from its
[AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (AA)
graph and its ecDNA interval(s). The amplicon extent is taken from AmpliconSuite
itself (AA/AC's own call) and the feature is classified by the structure of its boundary junction.

## How it works — the episome model

1. **Extent** = the AmpliconClassifier ecDNA interval BED
   (`*_ecDNA_<k>_intervals.bed`), or the dominant AA cycle from a `*_cycles.txt`
   (via [`aa_cycle_intervals()`](https://alhafidzhamdan.github.io/EpiTracer/reference/aa_cycle_intervals.html)).
   AC's exclusive end is handled internally, so the boundary test needs no fuzzy
   tolerance.
2. **Boundary DUP** = an intra-chromosomal duplication whose ends **bracket** the
   ecDNA's outer boundary on its chromosome (`lo ≤ extent start`,
   `hi ≥ extent end`). SV class is assigned **after coordinate-sorting the two
   breakends** (AA writes edges in arbitrary orientation): a `lo-`/`hi+` pair is
   the head-to-tail circulariser.
3. **Episomal** if ≥ 1 boundary DUP exists. The reported boundary DUP is the one
   with the **highest read-pair support** *among the spanning DUPs* — a selector,
   not a gate: a higher-support DEL/INV elsewhere does not veto the call.
4. **Chimeric episomal** = every chromosome of a multi-chromosome ecDNA is
   bracketed by its own boundary DUP **and** ≥ 1 translocation links them.
5. **Internal structure is ignored** — an episome may be internally chromothriptic
   (within a micronucleus); junction count is never a criterion.

Two guards reject geometry that is not a self-contained circle:

- **Occupancy guard** (`min_occupancy`, default `0.05`): the boundary-DUP
  bracketed span must be > 5 % covered by the ecDNA intervals — a lone DUP across
  a largely diploid region (e.g. a whole chromosome arm) is not a circle. The
  default was set by a PCAWG-vs-PURPLE concordance sweep.
- **Bridged/ambiguous downgrade** (`bridge_gap`, default `1 Mb`): if the
  amplicon's highest read-pair-support junction reaches an **amplified region on
  the same chromosome** beyond the extent, the dominant rearrangement
  co-mobilises outside sequence → `status = "bridged"`. Cross-chromosome
  translocations (the chimeric signature) are not downgraded.

## Input data

Two standard AmpliconSuite files per amplicon:

| file | fields used |
|---|---|
| AA `*_graph.txt` | `sequence` edges → CN segments (`PredictedCopyCount`); `discordant` edges → junctions (junction copy number, `NumberOfReadPairs` = read-pair support, homology) |
| AC `*_ecDNA_<k>_intervals.bed` | the ecDNA extent (one feature per BED) |

Chromosome names may carry a `chr` prefix or not, as long as the graph and BED
agree. Breakpoint homology is frequently absent from AA graphs (written `None`);
for a homology signature use [`call_simple_excision()`](call_simple_excision.md).

## Output

A one-row `data.table`: `episomal` (logical), `status`
(`episomal`/`bridged`/`complex`), `chimeric`, `n_chr`, the boundary-DUP record
(`bd_chr`, `bd_lo`, `bd_hi`, `bd_pe_support`, `bd_jcn`, `bd_homlen`),
`bd_occupancy` (the occupancy-guard metric), the bridging-junction record
(`bridge_pe_support`, `bridge_partner`), and `n_junctions`.

## Companion functions

- [`read_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_ampliconsuite.html)
  — parse a `*_graph.txt` into segments + junctions.
- [`aa_cycle_intervals()`](https://alhafidzhamdan.github.io/EpiTracer/reference/aa_cycle_intervals.html)
  — derive the extent from a bare AA `*_cycles.txt`.

See the [vignette](https://alhafidzhamdan.github.io/EpiTracer/articles/ampliconsuite.html)
for a worked, runnable example.
