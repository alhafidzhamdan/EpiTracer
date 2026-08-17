## ---------------------------------------------------------------------------
## EpiTracer validation — GBM39 episomal EGFR ecDNA (positive-control figure)
##
## Draws the structural rearrangements at the EGFR amplicon of GBM39, the
## canonical simple EGFR ecDNA line, straight from its public AmpliconArchitect
## reconstruction. The episome signature is visible at a glance: a single
## circularisation junction (DUP) spanning the highly amplified EGFR locus, with
## diploid flanks either side.
##
## DATA: AmpliconRepository "AmpliconSuite benchmarking set hg38"
##   (project 6a5a4a970664f9111b586742), sample GBM39 — ships the AA
##   *_graph.txt / *_cycles.txt and the genome-wide CNVkit *_CNV_CALLS.bed.
##   Download and extract, then point this script at the sample folder.
##
## USAGE (from the package root):
##   Rscript validation/plot_gbm39_egfr.R /path/to/results/samples/GBM39
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(EpiTracer); library(ggplot2) })
source("validation/aa_to_epitracer.R")

args <- commandArgs(trailingOnly = TRUE)
sd <- if (length(args)) args[1] else stop("pass the GBM39 sample folder")

graphs <- list.files(sd, "_amplicon[0-9]+_graph\\.txt$", full.names = TRUE, recursive = TRUE)
cnv    <- list.files(sd, "CNV_CALLS\\.bed$", full.names = TRUE, recursive = TRUE)
stopifnot(length(graphs) > 0, length(cnv) > 0)

pin <- aa_to_plot_inputs(graphs, "GBM39", cnv_bed = cnv[1])
message("SVs at EGFR amplicon: ", nrow(pin$sv_data),
        " (", paste(names(table(pin$sv_data$svclass)), table(pin$sv_data$svclass),
                    collapse = ", "), ")")

p <- plot_sv_linear(
  sample = "GBM39", cnv_data = pin$cnv_data, sv_data = pin$sv_data,
  chromosome = "chr7", chromosome_range = matrix(c(54.6e6, 56.2e6), nrow = 1)) +
  labs(title = "GBM39: episomal EGFR ecDNA",
       subtitle = "one circularisation junction (arc) spans the amplified EGFR locus; diploid flanks")

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/gbm39_egfr.png", p, width = 9, height = 5.2, dpi = 135, bg = "white")
message("Wrote validation/output/gbm39_egfr.png")
