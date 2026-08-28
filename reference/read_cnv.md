# Read and reformat copy-number segments for EpiTracer

Turns a copy-number table from (almost) any source into the columns
EpiTracer's plotters and callers expect: `sample`, `seqnames`, `start`,
`end`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`,
`minorAlleleCopyNumber`. Column names are matched case-insensitively, so
a raw PURPLE `*.purple.cnv.somatic.tsv` works directly. PURPLE's
occasional negative copy numbers are clamped to 0; the allele copy
numbers are derived if absent; and, because a PURPLE segment file
carries no ploidy (it lives in `*.purple.purity.tsv`), ploidy is
estimated from the segments unless supplied.

## Usage

``` r
read_cnv(x, sample = NULL, ploidy = NULL, name = NULL)
```

## Arguments

- x:

  A file path (TSV/CSV, `*.purple.cnv.somatic.tsv`, or a `.rds` of a
  GRanges/data.frame) or an in-memory data.frame/data.table.

- sample:

  Optional sample name. When given it **overrides** any `sample` column
  and the file name, labelling every row with this name – use it to give
  CN and SV files a common name when their file-naming differs.

- ploidy:

  Optional numeric ploidy. When given it overrides the estimate.

- name:

  Optional file name used for format detection and to derive the sample
  when `x` is a temporary path (e.g. a Shiny upload) whose own name is
  opaque.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
in EpiTracer's copy-number format.

## Examples

``` r
if (FALSE) { # \dontrun{
cnv <- read_cnv("DO13869T.purple.cnv.somatic.tsv", sample = "DO13869T")
} # }
```
