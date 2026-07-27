# Bundled reference data (hg38 / GRCh38)

Small default references used by `plot_sv_linear()` when `karyotype` /
`gene_coord` are not supplied.

- **`chr_info_hg38.rds`** — UCSC hg38 cytoband ideogram (`cytoBand` table:
  chromosome, start, end, `gieStain`). Source: UCSC Genome Browser.
- **`oncogene_coord_hg38.bed`** — coordinates (chr, start, end, strand, gene)
  for the default oncogene panel only, derived from Ensembl GRCh38 release 93
  gene annotation. This is a minimal subset kept small for distribution; supply
  your own `gene_coord` (a data.frame or BED path) to label other genes or to
  use a different genome build.
