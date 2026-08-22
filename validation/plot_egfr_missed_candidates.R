## Eyeball-confirmation montage for the 8 EGFR amplicons the flank test drops as
## candidate missed episomes (flank_sliver + flank_pass_other from
## egfr_nonepisomal_reasons.tsv). Each panel is drawn at its full amplicon extent
## and titled with the failure reason, so a genuine boundary-spanning DUP arc can
## be confirmed by eye before any caller change.
## USAGE: Rscript validation/plot_egfr_missed_candidates.R sv.rds cn.rds
suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(ggplot2); library(patchwork)
})
args <- commandArgs(trailingOnly = TRUE)
sv <- as.data.table(readRDS(args[1])); cn <- as.data.table(readRDS(args[2]))
rs <- fread("validation/output/egfr_nonepisomal_reasons.tsv")
cand <- rs[reason %in% c("flank_sliver", "flank_pass_other")][order(reason, -maxCN)]
message(nrow(cand), " candidate missed episomes to plot")

EGFR <- GenomicRanges::GRanges("chr7", IRanges::IRanges(55019017, 55211628))
panel <- function(s, reason, maxcn, fail_flank, fail_cn, fail_width, thresh) {
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
  ## extend to the boundary DUP breakends (SV BEDPE stores chrom prefix-LESS -> "7")
  dup <- sv[sample == s & chrom1 == "7" & chrom2 == "7" & svclass == "DUP" &
            pmin(start1, start2) <= up + 2e5 & pmax(start1, start2) >= lo - 2e5]
  if (nrow(dup)) { lo <- min(lo, dup$start1, dup$start2); up <- max(up, dup$start1, dup$start2) }
  pad <- max(0.12 * (up - lo), 2e5)
  sub <- if (reason == "flank_sliver")
           sprintf("flank sliver: %s CN %.1f / %s bp (thr %.1f)", fail_flank, fail_cn, fail_width, thresh)
         else "flanks pass on sorted recompute (ordering)"
  tryCatch(
    plot_sv_linear(sample = s, cnv_data = cn[sample == s], sv_data = sv[sample == s],
                   genome = "hg38", chromosome = "chr7",
                   chromosome_range = matrix(c(lo - pad, up + pad), nrow = 1)) +
      labs(title = sprintf("%s  (CN %d)", s, maxcn), subtitle = sub) +
      theme(plot.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 7, colour = "#b2182b"),
            axis.title = element_text(size = 7), axis.text = element_text(size = 6)),
    error = function(e) { message("  ! ", s, ": ", conditionMessage(e)); NULL })
}
plots <- Filter(Negate(is.null), lapply(seq_len(nrow(cand)), function(i) {
  r <- cand[i]; panel(r$WGS_ID, r$reason, r$maxCN, r$fail_flank, r$fail_cn, r$fail_width, r$thresh)
}))
fig <- wrap_plots(plots, ncol = 4) +
  plot_annotation(
    title = sprintf("Candidate MISSED EGFR episomes — dropped by the flank test (%d)", length(plots)),
    subtitle = "each has a dominant boundary DUP; flank test tripped by a sub-kb shoulder or boundary-ordering. Confirm the locus-spanning DUP arc.",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11, colour = "grey30")))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
ggsave("validation/output/gbm_egfr_missed_candidates.png", fig,
       width = 20, height = 3.4 * ceiling(length(plots) / 4), dpi = 110, bg = "white", limitsize = FALSE)
message("Wrote validation/output/gbm_egfr_missed_candidates.png (", length(plots), " panels)")
