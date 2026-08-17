## ---------------------------------------------------------------------------
## EpiTracer validation — montage of every episomal ecDNA called across CCLE
##
## Draws the rearrangements at each amplicon that `cell_line_benchmark.R` called
## episomal, from the public CCLE AmpliconArchitect reconstructions. Each panel
## spans EVERY chromosome the amplicon touches: an amplicon with inter-
## chromosomal (TRA) junctions is drawn as a multi-locus panel, so the
## translocation arcs are visible and the reader can judge whether the call is a
## genuine single-locus episome or a chimeric/complex amplicon.
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

## one locus per chromosome the amplicon touches (footprint + junction breakends,
## padded); ordered widest-first so the main locus leads.
loci_for <- function(inp) {
  amp <- inp$ecdna_gr; b <- inp$breakpoints_gr
  chrs <- unique(c(as.character(seqnames(amp)), as.character(seqnames(b))))
  rows <- lapply(chrs, function(ch) {
    p <- numeric()
    ai <- amp[as.character(seqnames(amp)) == ch]; if (length(ai)) p <- c(p, start(ai), end(ai))
    bi <- b[as.character(seqnames(b)) == ch];   if (length(bi)) p <- c(p, start(bi))
    if (!length(p)) return(NULL)
    lo <- min(p); hi <- max(p); pad <- max((hi - lo) * 0.2, 5e5)
    data.frame(chr = ch, start = max(1, lo - pad), end = hi + pad,
               span = (hi - lo) + 2 * pad, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, Filter(Negate(is.null), rows))
  rows[order(-rows$span), , drop = FALSE]
}

plot_one <- function(wgs, graphfile) {
  sd <- file.path(samples_root, wgs)
  g  <- file.path(sd, paste0(wgs, "_reconstruction_results"), graphfile)
  cp <- sub("_graph.txt$", "_cycles.txt", g)
  cnv <- file.path(sd, paste0(wgs, "_CNV_CALLS.bed"))
  inp <- aa_to_epitracer_inputs(g, cp, wgs, cnv_bed = cnv)
  if (!length(inp$ecdna_gr)) return(NULL)
  L <- loci_for(inp)
  has_tra <- any(inp$breakpoints_gr$svclass == "TRA")
  pin <- aa_to_plot_inputs(g, wgs, cnv_bed = cnv)
  tag <- sprintf("%s  (%s%s)", sub("_.*", "", wgs), paste(L$chr, collapse = "+"),
                 if (has_tra) ", TRA" else "")
  plot_sv_linear(sample = wgs, cnv_data = pin$cnv_data, sv_data = pin$sv_data,
                 genome = "hg38", loci = L[, c("chr", "start", "end")]) +
    labs(title = tag) +
    theme(plot.title = element_text(size = 9.5, face = "bold"),
          axis.title = element_text(size = 8), axis.text = element_text(size = 6.5))
}

plots <- Filter(Negate(is.null), lapply(seq_len(nrow(epi)), function(i) {
  message("plotting ", epi$WGS_ID[i])
  tryCatch(plot_one(epi$WGS_ID[i], basename(epi$amplicon[i])),
           error = function(e) { message("  ! ", conditionMessage(e)); NULL })
}))

fig <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = sprintf("All %d episomal ecDNA called by EpiTracer across 329 CCLE cell lines", length(plots)),
    subtitle = "each panel spans every chromosome the amplicon touches; amplicons with a TRA are drawn multi-locus so translocation arcs show",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/ccle_episomal_montage.png", fig,
       width = 16, height = 22, dpi = 120, bg = "white", limitsize = FALSE)
message("Wrote validation/output/ccle_episomal_montage.png (", length(plots), " panels)")
