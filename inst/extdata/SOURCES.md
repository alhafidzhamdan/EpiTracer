# Bundled reference data (hg38 / hg19 / mm10)

Default references used by `plot_sv_linear()` and `plot_sv_reconstruction()`
when `karyotype` / `gene_coord` are not supplied. The build is selected with the
`genome` argument (`"hg38"` default, `"hg19"`, `"mm10"`).

- **`chr_info_{hg38,hg19,mm10}.rds`** — UCSC `cytoBand` ideogram (chromosome,
  start, end, band name, `gieStain`), main chromosomes only. Source: UCSC Genome
  Browser (hg38/hg19: chr1-22,X,Y; mm10: chr1-19,X,Y).
- **`oncogene_coord_{hg38,hg19,mm10}.bed`** — coordinates (chr, start, end,
  strand, gene) for the default oncogene panel. hg38/hg19 use the human symbols;
  mm10 uses the mouse orthologs (`Trp53` for TP53, otherwise the title-cased
  symbol). A minimal panel kept small for distribution — supply your own
  `gene_coord` (a data.frame or BED path) to label other genes. The hg38 panel
  was derived from Ensembl GRCh38 release 93.

The hg19 and mm10 files are regenerated from the hg38 panel by
`data-raw/make_reference_data.R` (UCSC REST API for cytobands and mm10 refGene;
MyGene.info `genomic_pos_hg19` with a GRCh37 Ensembl fallback for hg19 genes).
