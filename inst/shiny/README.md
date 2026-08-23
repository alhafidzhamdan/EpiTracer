# EpiTracer web app

Upload copy-number and structural-variant calls, plot them with
`plot_sv_linear()`, and (optionally) call amplicon-formation mechanisms
(`call_simple_excision`, `call_brf`, `call_micronucleation`, `call_bfb`,
`call_translocation_bridge_amp`).

## Run locally

```r
# from the package root, with EpiTracer installed:
install.packages(c("shiny", "bslib", "DT"))
shiny::runApp("inst/shiny")
```

Then upload the two demo files in `inst/shiny/example/` to try it, or your own data.

## Input formats

Two files, each **TSV / CSV / BEDPE** or a PURPLE-style **`.rds`**. Column names are
matched case-insensitively against common aliases, so most callers' output works
as-is.

**Copy-number segments** — one row per segment. Required: `seqnames` (or
`chrom`/`chr`), `start`, `end`, `copyNumber` (or `cn`). Optional but recommended:
`sample`, `ploidy` (default 2), `majorAlleleCopyNumber`, `minorAlleleCopyNumber`
(minor allele is needed for LOH-based calls; defaulted if absent). **Segments must
be non-overlapping** (each `start` = previous `end` + 1).

**Copy number — raw PURPLE** — a `.purple.cnv.somatic.tsv` is auto-detected and
reformatted: `chromosome`→`seqnames`, negative copy numbers floored to 0, the
**sample taken from the filename**, and **ploidy estimated** (length-weighted mean
copy number, since ploidy lives in the separate `.purple.purity.tsv`). Override the
estimate with the Ploidy box if you have the real value.

**Structural variants (BEDPE)** — one row per junction. Required: `chrom1`,
`start1`, `chrom2`, `start2`, `svclass`. Optional: `strand1`, `strand2`, `VF`,
`JCN`, `name`, `homlen`, `sample`. `svclass` is normalised to
`DEL`/`DUP`/`INS`/`TRA`/`h2hINV`/`t2tINV` (inter-chromosomal → `TRA`; ambiguous
inversions resolved from strands). **Keep `VF` (supporting fragments) and `JCN`
(junction copy number)** — the mechanism callers use them; a plain BEDPE without
them degrades the calls.

**Structural variants — PURPLE/GRIDSS VCF** — a `.vcf`/`.vcf.gz` of breakend (BND)
records is parsed into the BEDPE above, **preserving `VF` (INFO or FORMAT) and `JCN`
(`PURPLE_JCN`)** rather than dropping them as a coordinate-only conversion would.
Mate breakends are de-duplicated and strand/orientation is read from the breakend
notation. Single-sample PURPLE CN and SV files whose names differ (e.g.
`DO13869T` vs `DO13869T1.svs`) are matched automatically.

## What it does

- **Sample** picker (intersection of the two files).
- **Region**: auto (all amplified loci), a detected amplicon, a whole chromosome,
  or a custom `chr:start-end`.
- **Plot** with `plot_sv_linear()` (download as PDF).
- **Mechanisms** tab: per-amplicon mechanism table (episomal ecDNA, BRF,
  micronucleation, BFB, translocation-bridge), with junction-homology class and
  translocation-bridge confidence.
- **Detected amplicons** tab: the copy-number seeds.

## Deploy to shinyapps.io

See `deploy.R`. In short: install EpiTracer from GitHub (so the server can rebuild
it), connect your shinyapps.io account, then `source("inst/shiny/deploy.R")`.

## Files

- `app.R` — Shiny UI + server.
- `helpers.R` — upload readers, column normalisation, caller wrappers (unit-testable).
- `example/` — a small **synthetic** demo (`DEMO1`): an EGFR episomal ecDNA on chr7
  and a fold-back array on chr8. `make_demo.R` regenerates it.
