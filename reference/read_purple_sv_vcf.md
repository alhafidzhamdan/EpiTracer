# Convert a PURPLE / GRIDSS breakend VCF to BEDPE, keeping VF and JCN

Parses breakend (`BND`) records from a PURPLE/GRIDSS `*.vcf`/`*.vcf.gz`
into the paired-breakend BEDPE EpiTracer uses, **preserving `VF`**
(variant fragments, from the `VF` INFO field or the tumour sample's `VF`
FORMAT field) and **`JCN`** (from `PURPLE_JCN`). A plain coordinate-only
VCF-to-BEDPE conversion drops these, which degrades the mechanism calls.
Mate breakends are de-duplicated and strand / orientation is read from
the breakend ALT notation (`N[chr:pos[` etc.), so `svclass` follows the
usual convention (`+/-`=DEL, `-/+`=DUP, `+/+`=h2hINV, `-/-`=t2tINV,
inter-chromosomal=TRA).

## Usage

``` r
read_purple_sv_vcf(file)
```

## Arguments

- file:

  Path to a `*.vcf` or `*.vcf.gz`.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with BEDPE columns plus `VF`, `JCN`, `homlen`.
