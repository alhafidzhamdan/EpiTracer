## ---------------------------------------------------------------------------
## Montage of the NON-clean episomal calls flagged by the focality audit
## (validation/scan_episomal_focality.R): arm-scale and/or fold-back-dominated.
## Each is drawn zoomed to its own amplicon footprint (amplified CN segments,
## padded), so the local structure -- boundary DUP vs internal fold-backs, CN
## oscillation -- is visible.
##
## USAGE:
##   Rscript validation/plot_suspect_episomes_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

scan <- fread("validation/output/episomal_focality_scan.tsv")[foldback_dominated == TRUE]
scan[, cat := "fold-back-dominated"]
setorder(scan, -max_fb_vf)
message("fold-back-dominated episomal calls to plot: ", nrow(scan))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
gene_coord <- data.frame(chr = onc$chr, start = onc$start, end = onc$end,
                         strand = onc$strand, gene = onc$gene)
cl <- load_chrom_lengths("hg38")

build_inputs <- function(s) {
  cs <- cn[sample == s]; vs <- sv_raw[sample == s]; if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  pcol <- base$b[match(cs$seqnames, base$seqnames)]; pcol[is.na(pcol)] <- 2
  cs$ploidy <- pcol
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s, copyNumber = cs$copyNumber,
                    ploidy = pcol, majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
                    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- as.character(vs$name)
  mk <- function(ch, pos, st) data.frame(chr = pfx(ch), pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = st, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN, VF = vs$VF,
    insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen))
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4)); ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, bgr = bgr, cs = cs, vs = vs)
}

plots <- list()
for (s in unique(scan$WGS_ID)) {
  inp <- build_inputs(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  for (i in which(scan$WGS_ID == s)) {
    id <- scan$ID[i]; ch <- scan$chr[i]
    seed <- ecdna_gr[ecdna_gr$ID == id]; if (!length(seed)) next
    amp <- cs[seqnames == ch & copyNumber > 3 * ploidy &
                end >= min(start(seed)[as.character(seqnames(seed)) == ch]) &
                start <= max(end(seed)[as.character(seqnames(seed)) == ch])]
    if (!nrow(amp)) next
    lo <- min(amp$start); hi <- max(amp$end); pad <- max(3e5, 0.25 * (hi - lo))
    cr <- matrix(c(max(1, lo - pad), min(as.numeric(cl[[ch]]), hi + pad)), nrow = 1)
    gg <- if (nchar(scan$gene[i])) scan$gene[i] else NA
    ttl <- sprintf("%s | %s | %s | %s | fp %.1f Mb, %d FB, bDUP %s vs FB %s",
                   s, ch, ifelse(is.na(gg) | gg == "", "(no oncogene)", gg), scan$cat[i],
                   scan$footprint_mb[i], scan$n_foldbacks[i], scan$bdup_vf[i], scan$max_fb_vf[i])
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = as.data.frame(cs[seqnames == ch]),
                     sv_data = as.data.frame(vs), genome = "hg38",
                     chromosome = ch, chromosome_range = cr, events = c("amp", "loh"),
                     genes_to_highlight = if (is.na(gg)) NULL else unlist(strsplit(gg, ",")),
                     gene_coord = gene_coord, save = FALSE, verbose = FALSE),
      error = function(e) { message("  panel failed ", s, ": ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      plots[[length(plots) + 1]] <- p + ggtitle(ttl) +
        theme(plot.title = element_text(size = 8, face = "bold"))
      message(ttl)
    }
  }
}
message("panels: ", length(plots))

if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 4; nrow <- 4; per <- ncol * nrow
  outpdf <- "validation/output/foldback_dominated_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 8 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Fold-back-dominated episomal calls %d-%d of %d (internal fold-back VF > boundary-DUP VF; amplicon zoom)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = theme(plot.title = element_text(face = "bold", size = 14)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("patchwork not available or no panels; montage skipped")
