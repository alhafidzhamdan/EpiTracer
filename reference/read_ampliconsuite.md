# Read an AmpliconArchitect graph file

Parses an AmpliconArchitect `*_graph.txt` into copy-number segments
(from the `sequence` edges) and structural-variant junctions (from the
`discordant` edges), the latter with a **coordinate-invariant** SV
class. Chromosome names are accepted with or without a `chr` prefix.

## Usage

``` r
read_ampliconsuite(graph_file)
```

## Arguments

- graph_file:

  Path to an AmpliconArchitect `*_graph.txt`.

## Value

A list with three
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)s:

- `segments`:

  `chr`, `start`, `end`, `cn` (predicted copy number).

- `junctions`:

  `cA`,`pA`,`sA` and `cB`,`pB`,`sB` (the two breakends'
  chromosome/position/strand), `jcn` (junction predicted copy number),
  `pe_support` (number of supporting read pairs), `homlen` (breakpoint
  homology length, `NA` if absent), `svclass` (`DUP`/`DEL`/`INV`/`TRA`)
  and `lo`/`hi` (the lower/upper breakend coordinate).

- `source`:

  `chr`, `pos`, `cn` for each AA *source* edge – an open end where
  copy-number flow enters/leaves the amplicon via the rest of the genome
  instead of looping back through a junction (a self-contained circle
  has none). Used by the structural QC of
  [`call_episomal_from_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_from_ampliconsuite.md).

Any table may have zero rows.

## See also

[`call_episomal_from_ampliconsuite()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_from_ampliconsuite.md)
