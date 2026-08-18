## Montage of every EGFR episomal ecDNA in the cohort (one panel per sample,
## drawn at the sample's own EGFR-amplicon extent). Reads the cohort SV BEDPE +
## allele-specific CN and the episomal-catalogue table produced by
## call_episomal_cohort.R. USAGE:
##   Rscript validation/plot_gbm_egfr_episomes.R sv_bedpe.rds cn_segments.rds
suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(ggplot2); library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
sv <- as.data.table(readRDS(args[1])); cn <- as.data.table(readRDS(args[2]))
cat_tab <- fread("validation/output/gbm_episomal_catalogue.tsv")
egfr <- sort(unique(cat_tab[gene == "EGFR"]$WGS_ID))
message(length(egfr), " EGFR-episome samples")

panel <- function(s) {
  cs <- cn[sample == s & seqnames == "chr7"]
  pl <- stats::median(cs$ploidy)
  amp <- cs[copyNumber > 3 * pl & end > 54.5e6 & start < 55.6e6]
  if (!nrow(amp)) return(NULL)
  lo <- min(amp$start); up <- max(amp$end); pad <- max(0.2 * (up - lo), 2e5)
  mx <- round(max(amp$copyNumber))
  tryCatch(
    plot_sv_linear(sample = s, cnv_data = cn[sample == s], sv_data = sv[sample == s],
                   genome = "hg38", chromosome = "chr7",
                   chromosome_range = matrix(c(lo - pad, up + pad), nrow = 1)) +
      labs(title = sprintf("%s  (CN %d)", s, mx)) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            axis.title = element_text(size = 7), axis.text = element_text(size = 6)),
    error = function(e) { message("  ! ", s, ": ", conditionMessage(e)); NULL })
}
plots <- Filter(Negate(is.null), lapply(egfr, panel))
fig <- wrap_plots(plots, ncol = 5) +
  plot_annotation(
    title = sprintf("Episomal EGFR ecDNA in the glioblastoma cohort — all %d samples", length(plots)),
    subtitle = "each panel drawn at its own EGFR-amplicon extent; the arc spanning the locus is the circularisation boundary DUP",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/gbm_egfr_episomes_all.png", fig,
       width = 22, height = 3.1 * ceiling(length(plots) / 5), dpi = 105, bg = "white", limitsize = FALSE)
message("Wrote validation/output/gbm_egfr_episomes_all.png (", length(plots), " panels)")
