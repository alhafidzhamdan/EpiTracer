## Montage of every EGFR amplicon in the cohort that call_simple_excision() did
## NOT classify as episomal, drawn at each sample's full EGFR-amplicon extent (so
## missed episomes vs genuinely non-episomal EGFR amplifications can be told
## apart). Reads the cohort SV BEDPE + allele-specific CN and the caller output
## (cohort_original_caller_per_chr.tsv). USAGE:
##   Rscript validation/plot_gbm_egfr_nonepisomal.R sv_bedpe.rds cn_segments.rds
suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(ggplot2); library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
sv <- as.data.table(readRDS(args[1])); cn <- as.data.table(readRDS(args[2]))
calls <- fread("validation/output/cohort_original_caller_per_chr.tsv")
egfr <- calls[grepl("(^|,)EGFR(,|$)", genes)]
epi_samples <- unique(egfr[episomal == TRUE]$WGS_ID)
samples <- sort(setdiff(unique(egfr$WGS_ID), epi_samples))
message(length(samples), " EGFR-amplicon samples not classified episomal")

EGFR <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628))
panel <- function(s) {
  cs <- cn[sample == s & seqnames == "chr7"]
  pl <- stats::median(cs$ploidy)
  amp <- cs[copyNumber > 3 * pl]
  if (!nrow(amp)) return(NULL)
  red <- GenomicRanges::reduce(GenomicRanges::GRanges("chr7", IRanges::IRanges(amp$start, amp$end)),
                               min.gapwidth = 1e6)
  hit <- red[GenomicRanges::countOverlaps(red, EGFR) > 0]
  if (!length(hit)) hit <- red[which.min(abs((start(red) + end(red)) / 2 - 55.1e6))]
  if (!length(hit)) return(NULL)
  lo <- min(GenomicRanges::start(hit)); up <- max(GenomicRanges::end(hit))
  dup <- sv[sample == s & chrom1 == "chr7" & chrom2 == "chr7" & svclass == "DUP" &
            pmin(start1, start2) <= up + 2e5 & pmax(start1, start2) >= lo - 2e5]
  if (nrow(dup)) { lo <- min(lo, dup$start1, dup$start2); up <- max(up, dup$start1, dup$start2) }
  pad <- max(0.12 * (up - lo), 2e5)
  mx <- round(max(cs[start <= up & end >= lo]$copyNumber))
  tryCatch(
    plot_sv_linear(sample = s, cnv_data = cn[sample == s], sv_data = sv[sample == s],
                   genome = "hg38", chromosome = "chr7",
                   chromosome_range = matrix(c(lo - pad, up + pad), nrow = 1)) +
      labs(title = sprintf("%s  (CN %d)", s, mx)) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            axis.title = element_text(size = 7), axis.text = element_text(size = 6)),
    error = function(e) { message("  ! ", s, ": ", conditionMessage(e)); NULL })
}
plots <- Filter(Negate(is.null), lapply(samples, panel))
fig <- wrap_plots(plots, ncol = 6) +
  plot_annotation(
    title = sprintf("EGFR amplicons NOT classified episomal — glioblastoma cohort (%d samples)", length(plots)),
    subtitle = "detected EGFR amplicons that failed the episome test; each drawn at its full amplicon extent",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/gbm_egfr_nonepisomal.png", fig,
       width = 25, height = 2.9 * ceiling(length(plots) / 6), dpi = 100, bg = "white", limitsize = FALSE)
message("Wrote validation/output/gbm_egfr_nonepisomal.png (", length(plots), " panels)")
