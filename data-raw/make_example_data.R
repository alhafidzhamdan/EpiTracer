## Generate the small, self-contained example datasets shipped with EpiTracer.
##
## Both objects describe the SAME toy amplicon: a textbook *episomal* EGFR ecDNA
## on chr7 in a fictional sample "EXAMPLE01" — a focal amplification bounded by a
## high-VF duplication breakpoint, excised from otherwise diploid flanks, leaving
## a deletion "excision scar". The data are synthetic (no patient data) so the
## package is fully reproducible and the vignette runs anywhere in seconds.
##
## Run from the package root:  Rscript data-raw/make_example_data.R
## (re)creates data/ex_caller_inputs.rda and data/ex_plot_inputs.rda.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

## ---------------------------------------------------------------------------
## 1. Inputs for call_episomal_ecdna(): AmpliconArchitect + PURPLE-style GRanges
## ---------------------------------------------------------------------------

## AmpliconArchitect ecDNA amplicon region (needs $ID, $WGS_ID)
ecdna_gr <- GRanges(
  "chr7", IRanges(55000000, 55500000),
  ID = "EXAMPLE01_EGFR", WGS_ID = "EXAMPLE01"
)

## PURPLE SV breakpoints: a duplication at each amplicon boundary (the
## circularisation junction, high VF) plus a low-VF deletion pair (excision scar)
bp <- data.table(
  seqnames   = "chr7",
  start      = c(54900000, 55000000, 55500000, 55600000),
  end        = c(54900000, 55000000, 55500000, 55600000),
  WGS_ID     = "EXAMPLE01",
  event      = c("DEL1", "DUP1", "DUP1", "DEL1"),
  svclass    = c("DEL",  "DUP",  "DUP",  "DEL"),
  PURPLE_AF  = c(0.4,    0.9,    0.9,    0.4),
  PURPLE_JCN = c(1,      40,     40,     1),
  VF         = c(200,    1000,   1000,   200),
  PURPLE_CN  = c(2,      50,     50,     2),
  insLen     = 0L,
  HOMLEN     = 0L
)
breakpoints_gr <- makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)

## PURPLE allele-specific copy-number segments: amplified core, diploid flanks
cnv <- data.table(
  seqnames = "chr7",
  start    = c(40000000, 55000000, 55500001),
  end      = c(54999999, 55500000, 70000000),
  sample   = "EXAMPLE01",
  copyNumber            = c(2,  50, 2),
  ploidy                = 2,
  majorAlleleCopyNumber = c(1,  49, 1),
  minorAlleleCopyNumber = c(1,  1,  1)
)
cnv_gr <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

## Cancer-gene loci used to annotate the boundary breakpoints
cancer_genes_gr <- GRanges(
  "chr7", IRanges(55019017, 55211628), gene = "EGFR"
)

ex_caller_inputs <- list(
  ecdna_gr        = ecdna_gr,
  breakpoints_gr  = breakpoints_gr,
  cnv_gr          = cnv_gr,
  cancer_genes_gr = cancer_genes_gr
)

## ---------------------------------------------------------------------------
## 2. Inputs for plot_sv_linear(): the same amplicon in the plotter's format
##    (karyotype + gene_coord fall back to the bundled hg38 references)
## ---------------------------------------------------------------------------

cnv_data <- data.frame(
  sample   = "EXAMPLE01",
  seqnames = "chr7",
  start    = c(40e6, 55e6, 55.5e6 + 1),
  end      = c(54.99e6, 55.5e6, 70e6),
  copyNumber            = c(2, 30, 2),
  ploidy                = 2,
  majorAlleleCopyNumber = c(1, 29, 1),
  minorAlleleCopyNumber = c(1, 1,  1)
)

sv_data <- data.frame(
  chrom1 = "7", start1 = 55e6,
  chrom2 = "7", start2 = 55.5e6,
  strand1 = "-", strand2 = "+",
  svclass = "DUP", VF = 1000, JCN = 40, sample = "EXAMPLE01"
)

wgd_data <- data.frame(sample = "EXAMPLE01", Polyploidy = "No")

ex_plot_inputs <- list(
  cnv_data = cnv_data,
  sv_data  = sv_data,
  wgd_data = wgd_data
)

## ---------------------------------------------------------------------------
## 3. A kataegis example for the optional SNV panel of plot_sv_linear()
##    A tight cluster of C>N substitutions (short intermutation distances -> a
##    rainfall-plot dip) over the EGFR amplicon, plus a few scattered background
##    SNVs for contrast. SNV column convention: sampleID / allelic_freq / type.
## ---------------------------------------------------------------------------

set.seed(2)
## scattered background SNVs (outside the amplicon) -> long intermutation distances
bg_pos  <- sort(round(c(runif(4, 54.10e6, 54.90e6), runif(4, 55.60e6, 55.95e6))))
## a kataegis shower over the amplicon: ~25 SNVs with 0.2-4 kb gaps, C>N (APOBEC-like)
gaps    <- round(runif(24, 200, 4000))
kat_pos <- 55240000L + as.integer(c(0, cumsum(gaps)))
nb <- length(bg_pos); nk <- length(kat_pos)

ex_snv <- data.frame(
  sampleID     = "EXAMPLE01",
  seqnames     = "chr7",
  start        = c(bg_pos, kat_pos),
  allelic_freq = c(runif(nb, 0.20, 0.35), runif(nk, 0.35, 0.50)),
  variant_cn   = c(rep(2L, nb), rep(24L, nk)),
  type         = "SNV",
  ref          = c(rep("T", nb), rep("C", nk)),           # kataegis: C>N (APOBEC-like)
  mut          = c(rep("C", nb), rep(c("T", "G"), length.out = nk)),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------------
usethis::use_data(ex_caller_inputs, ex_plot_inputs, ex_snv, overwrite = TRUE)
