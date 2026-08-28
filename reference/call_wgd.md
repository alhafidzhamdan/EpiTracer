# Call whole-genome doubling (WGD) from allele-specific copy number

Applies PURPLE's whole-genome-doubling rule: a sample is WGD when its
**major-allele copy number exceeds 1.5 across at least half of the
bases** on **at least 11 of the 22 autosomes**. Operates directly on the
allele-specific copy-number segments EpiTracer already uses, so WGD
status needs no separate input file.

## Usage

``` r
call_wgd(cnv_gr)
```

## Arguments

- cnv_gr:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html),
  `data.frame` or
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  of copy-number segments carrying `majorAlleleCopyNumber` (and
  `seqnames`/`start`/`end`, taken from the ranges for a GRanges). An
  optional `sample` column splits multi-sample input; when absent the
  whole table is treated as one sample. Chromosome names may be with or
  without the `chr` prefix.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with one row per `sample`: `n_wgd_autosomes` (how many autosomes meet
the \>1.5-over-half-bases test) and `wgd` (logical;
`n_wgd_autosomes >= 11`).

## See also

[`plot_sv_circos()`](https://alhafidzhamdan.github.io/EpiTracer/reference/plot_sv_circos.md)
