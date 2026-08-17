## ---------------------------------------------------------------------------
## EpiTracer validation — montage of every episomal ecDNA called across CCLE
##
## Draws the rearrangements at each amplicon that `cell_line_benchmark.R` called
## episomal, from the public CCLE AmpliconArchitect reconstructions. The arc
## spanning each amplified locus is the circularisation junction; simple
## episomes show one clean arc over diploid flanks, busier panels carry more
## internal structure.
##
## Requires: the CCLE benchmark output (validation/output/cell_line_episomal_calls.tsv)
## and the extracted CCLE archive (AmpliconRepository project 6a5a432efaab0afcb5586742).
##
## USAGE (from the package root):
##   Rscript validation/plot_ccle_episomal.R /path/to/results/samples
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(GenomicRanges); library(data.table)
  library(ggplot2); library(patchwork)
})
source("validation/aa_to_epitracer.R")

args <- commandArgs(trailingOnly = TRUE)
samples_root <- if (length(args)) args[1] else stop("pass the CCLE results/samples folder")
epi <- fread("validation/output/cell_line_episomal_calls.tsv")[episomal == TRUE]

## bundled hg38 oncogene panel, to centre each panel on its oncogene
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
onc_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)

plot_one <- function(wgs, graphfile) {
  sd <- file.path(samples_root, wgs)
  g  <- file.path(sd, paste0(wgs, "_reconstruction_results"), graphfile)
  cp <- sub("_graph.txt$", "_cycles.txt", g)
  cnv <- file.path(sd, paste0(wgs, "_CNV_CALLS.bed"))
  inp <- aa_to_epitracer_inputs(g, cp, wgs, cnv_bed = cnv)
  amp <- inp$ecdna_gr
  if (!length(amp)) return(NULL)
  ## prefer the amplicon interval overlapping a panel oncogene; else the widest
  hits <- amp[countOverlaps(amp, onc_gr) > 0]
  pick <- if (length(hits)) hits else amp
  w  <- tapply(width(pick), as.character(seqnames(pick)), sum)
  mc <- names(w)[which.max(w)]
  onchr <- pick[as.character(seqnames(pick)) == mc]
  lo <- min(start(onchr)); up <- max(end(onchr)); span <- up - lo
  pad <- max(span * 0.2, 3e5)
  rng <- matrix(c(lo - pad, up + pad), nrow = 1)
  pin <- aa_to_plot_inputs(g, wgs, cnv_bed = cnv)
  plot_sv_linear(sample = wgs, cnv_data = pin$cnv_data, sv_data = pin$sv_data,
                 genome = "hg38", chromosome = mc, chromosome_range = rng) +
    labs(title = sprintf("%s  (%s)", sub("_.*", "", wgs), mc)) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 8), axis.text = element_text(size = 7))
}

plots <- Filter(Negate(is.null), lapply(seq_len(nrow(epi)), function(i) {
  message("plotting ", epi$WGS_ID[i])
  tryCatch(plot_one(epi$WGS_ID[i], basename(epi$amplicon[i])),
           error = function(e) { message("  ! ", conditionMessage(e)); NULL })
}))

fig <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = sprintf("All %d episomal ecDNA called by EpiTracer across 329 CCLE cell lines", length(plots)),
    subtitle = "each panel: copy number + SV arcs at the amplicon; the arc spanning the locus is the circularisation junction",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/ccle_episomal_montage.png", fig,
       width = 15, height = 15, dpi = 120, bg = "white", limitsize = FALSE)
message("Wrote validation/output/ccle_episomal_montage.png (", length(plots), " panels)")
