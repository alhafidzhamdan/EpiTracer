## ---------------------------------------------------------------------------
## EpiTracer validation — montage of every per-locus EPISOMAL ecDNA in CCLE
##
## Draws each amplified locus that classify_loci.R tags `episomal`, one panel per
## amplicon (fused episomes of the same amplicon share a multi-locus panel). Reads
## the per-locus benchmark output (validation/output/cell_line_loci.tsv) and the
## extracted CCLE archive.
##
## USAGE (from the package root):
##   Rscript validation/plot_episomal_loci.R /path/to/results/samples
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(ggplot2); library(patchwork)
})
source("validation/aa_to_epitracer.R")

args <- commandArgs(trailingOnly = TRUE)
samples_root <- if (length(args)) args[1] else stop("pass the CCLE results/samples folder")
loci <- fread("validation/output/cell_line_loci.tsv")[mechanism == "episomal"]
loci[, key := paste(sample, amplicon)]

plots <- lapply(split(loci, by = "key"), function(g) {
  wgs <- g$sample[1]
  sd  <- file.path(samples_root, wgs)
  gph <- file.path(sd, paste0(wgs, "_reconstruction_results"), paste0(wgs, "_", g$amplicon[1], "_graph.txt"))
  cnv <- file.path(sd, paste0(wgs, "_CNV_CALLS.bed"))
  if (!file.exists(gph)) return(NULL)
  pin <- aa_to_plot_inputs(gph, wgs, cnv_bed = cnv)
  L <- data.frame(chr = g$chr,
                  start = pmax(1, g$start - pmax(0.2 * (g$end - g$start), 3e5)),
                  end = g$end + pmax(0.2 * (g$end - g$start), 3e5))
  onc <- paste(unique(g$oncogenes[nzchar(g$oncogenes)]), collapse = ",")
  tag <- sprintf("%s  (%s%s%s)", sub("_.*", "", wgs), paste(g$chr, collapse = "+"),
                 if (nzchar(onc)) paste0(", ", onc) else "",
                 if (any(g$fused)) ", fused" else "")
  tryCatch(
    plot_sv_linear(sample = wgs, cnv_data = pin$cnv_data, sv_data = pin$sv_data,
                   genome = "hg38", loci = L) +
      labs(title = tag) +
      theme(plot.title = element_text(size = 9.5, face = "bold"),
            axis.title = element_text(size = 8), axis.text = element_text(size = 6.5)),
    error = function(e) { message("  ! ", wgs, ": ", conditionMessage(e)); NULL })
})
plots <- Filter(Negate(is.null), plots)

fig <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = sprintf("Episomal ecDNA across 329 CCLE cell lines (%d loci; per-locus classifier)", nrow(loci)),
    subtitle = "each panel: an amplicon's episomal locus/loci (copy number + SV arcs); fused = two episomes joined by a subclonal translocation",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/ccle_episomal_loci_montage.png", fig,
       width = 16, height = 20, dpi = 120, bg = "white", limitsize = FALSE)
message("Wrote validation/output/ccle_episomal_loci_montage.png (", length(plots), " panels)")
