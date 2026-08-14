# Classify a single ecDNA amplicon as episomal or not

Internal worker called once per amplicon by
[`call_episomal_ecdna()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_ecdna.md).
It annotates the structural variant breakpoints falling within (or just
outside) an ecDNA amplicon, then applies the episome heuristic:

1.  find duplication (DUP) breakpoints that sit at the amplicon
    boundaries and are themselves amplified;

2.  require the boundary DUP to carry the highest variant fraction (VF)
    of any DUP in the amplicon;

3.  require the chromosomal segments immediately flanking the boundaries
    to be non-gained (consistent with a circular episome excised from an
    otherwise diploid region);

4.  flag a shared flanking deletion as a candidate excision scar.

## Usage

``` r
classify_amplicon_episomal(
  this_amplicon_id,
  ecdna_gr,
  breakpoints_gr,
  cnv_gr,
  cancer_genes_gr,
  ext = 1e+07,
  verbose = FALSE
)
```

## Arguments

- this_amplicon_id:

  Character scalar; a single value of `ecdna_gr$ID`.

- ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr:

  See
  [`call_episomal_ecdna()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_episomal_ecdna.md).

- ext:

  Integer; bp to extend the amplicon by when searching for boundary SVs.

- verbose:

  Logical; print progress/diagnostics.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
of annotated breakpoints for the amplicon (may have zero rows), or an
empty data.table if the amplicon has no breakpoints in range.
