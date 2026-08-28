# Read and reformat structural variants for EpiTracer

Accepts a BEDPE (with `VF`/`JCN`), a PURPLE/GRIDSS breakend VCF (parsed
by
[`read_purple_sv_vcf()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_purple_sv_vcf.md)),
or a `.rds`, and returns the columns EpiTracer uses: `chrom1`, `start1`,
`chrom2`, `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`,
`name`, `homlen`, `sample`. `svclass` is normalised to
`DEL`/`DUP`/`INS`/`TRA`/`h2hINV`/`t2tINV`. Keeping **VF** (supporting
fragments) and **JCN** (junction copy number) matters – the mechanism
callers use them.

## Usage

``` r
read_sv(x, sample = NULL, name = NULL)
```

## Arguments

- x:

  A file path (BEDPE TSV/CSV, `*.vcf`/`*.vcf.gz`, or `.rds`) or a
  data.frame/data.table.

- sample:

  Optional sample name; overrides any `sample` column and the file name
  (see
  [`read_cnv()`](https://alhafidzhamdan.github.io/EpiTracer/reference/read_cnv.md)).

- name:

  Optional file name for format detection / sample derivation.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
in EpiTracer's SV (BEDPE) format.
