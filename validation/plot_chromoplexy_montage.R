## ---------------------------------------------------------------------------
## Montage of every called CHROMOPLEXY cycle (closed balanced rearrangement
## chains; Baca et al., Cell 2013), one panel per cycle.
##
## For each cycle we draw the whole chromosomes it spans, with ONLY the cycle's
## own junctions overlaid (subset by the caller's `events` list), so the closed
## ring chr_a -> chr_b -> chr_c -> chr_a is legible rather than buried in the
## sample's full rearrangement set. amp + LOH events are shaded for context;
## oncogene labels are suppressed (chromoplexy is not amplicon-centric).
##
## USAGE (from the package root):
##   Rscript validation/plot_chromoplexy_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"

sv <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))
pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]
cn[, seqnames := pfx(as.character(seqnames))]
cl <- load_chrom_lengths("hg38"); chr_order <- names(cl)

cpx <- fread("validation/output/cohort_chromoplexy_per_chr.tsv")
setorder(cpx, -n_junctions, -n_chromosomes, WGS_ID)
message("Chromoplexy cycles: ", nrow(cpx), " in ", uniqueN(cpx$WGS_ID), " samples")

## suppress oncogene labels (empty 5-col gene_coord); chromoplexy isn't focal
empty_genes <- data.frame(chr = character(), start = numeric(), end = numeric(),
                          strand = character(), gene = character())

per_chr_ploidy <- function(cs) {
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  p <- base$b[match(cs$seqnames, base$seqnames)]; p[is.na(p)] <- 2; p
}
## reconstruct the per-sample event ids EXACTLY as call_episomal_cohort.R did, so
## a cycle's `events` ("SV3,SV11,...") map back to rows of the sample's bedpe
sample_events <- function(vs) if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))

plots <- list()
for (i in seq_len(nrow(cpx))) {
  s <- cpx$WGS_ID[i]
  chrs <- strsplit(cpx$chromosomes[i], ",")[[1]]
  chrs <- chrs[order(match(chrs, chr_order))]
  cs <- cn[sample == s]; vs <- sv[sample == s]
  if (!nrow(cs) || !nrow(vs)) next
  ev <- sample_events(vs)
  cyc_ev <- strsplit(cpx$events[i], ",")[[1]]
  vs_cyc <- as.data.frame(vs[ev %in% cyc_ev])            # only this cycle's junctions
  if (!nrow(vs_cyc)) next
  cnv_data <- as.data.frame(cs); cnv_data$ploidy <- per_chr_ploidy(cs)
  cr <- do.call(rbind, lapply(chrs, function(c) c(1, as.numeric(cl[[c]]))))
  ttl <- sprintf("%s | %s | %d junctions, %d bridges | frac_cp %.2f",
                 s, cpx$chromosomes[i], cpx$n_junctions[i], cpx$n_bridges[i], cpx$frac_cp[i])
  message(ttl)
  p <- tryCatch(
    plot_sv_linear(sample = s, cnv_data = cnv_data, sv_data = vs_cyc, genome = "hg38",
                   chromosome = chrs, chromosome_range = cr, events = c("amp", "loh"),
                   genes_to_highlight = NULL, gene_coord = empty_genes,
                   save = FALSE, verbose = FALSE),
    error = function(e) { message("  plot failed: ", conditionMessage(e)); NULL })
  if (!is.null(p)) plots[[length(plots) + 1]] <- p +
    ggplot2::ggtitle(ttl) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"))
}
message("panels: ", length(plots))

## paginated montage: 2 columns x 4 rows per page, one PDF
if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 2; nrow <- 4; per <- ncol * nrow
  outpdf <- "validation/output/chromoplexy_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 15 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Chromoplexy cycles %d-%d of %d (closed balanced rearrangement chains; Baca et al., Cell 2013)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("patchwork not available; montage skipped")
