## ---------------------------------------------------------------------------
## Plot every HIGH-CONFIDENCE TBA amplicon (dual-LOH: bridge arm on BOTH the
## amplicon and the partner chromosome, with the non-bridge arm spared).
##
## For each high-confidence amplicon we draw its amplicon chromosome together
## with its translocation partner chromosome, whole-chromosome, with amp + LOH
## events shaded so the two LOH bridge arms and the spared arm are visible.
##
## USAGE (from the package root):
##   Rscript validation/plot_tba_highconf.R \
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

onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cgr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
ce <- load_centromeres("hg38"); cl <- load_chrom_lengths("hg38")
chr_order <- names(cl)

hc <- fread("validation/output/cohort_original_caller_per_chr.tsv")[tb_high_confidence == TRUE]
message("High-confidence TBA amplicons: ", nrow(hc), " in ", uniqueN(hc$WGS_ID), " samples")

per_chr_ploidy <- function(cs) {
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  p <- base$b[match(cs$seqnames, base$seqnames)]; p[is.na(p)] <- 2; p
}

## Run the TBA caller for a sample and return the per-breakpoint annotation, so we
## can read the TRUE set of chromosomes each amplicon spans (its own breakend
## chromosomes AND its bridge partners) rather than guessing from a single seed.
run_tba <- function(s) {
  cs <- cn[sample == s]; vs <- sv[sample == s]
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s, copyNumber = cs$copyNumber,
    ploidy = per_chr_ploidy(cs), majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))
  mk <- function(ch, pos, st) data.frame(chr = ch, pos = pos, WGS_ID = s, event = ev, svclass = vs$svclass,
    bp_strand = st, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN, VF = vs$VF, insLen = 0,
    HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen))
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4)); ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  tba <- tryCatch(call_translocation_bridge_amp(ecdna_gr = NULL, breakpoints_gr = bgr, cnv_gr = cnv_gr,
    cancer_genes_gr = cgr, centromeres = ce, chrom_lengths = cl), error = function(e) NULL)
  if (is.null(tba)) return(NULL)
  list(t = as.data.table(as.data.frame(tba)), cnv_data = { d <- as.data.frame(cs); d$ploidy <- per_chr_ploidy(cs); d },
       vs = as.data.frame(vs))
}

plots <- list(); titles <- character()
for (s in unique(hc$WGS_ID)) {
  r <- run_tba(s); if (is.null(r)) { message(s, ": TBA run failed"); next }
  for (amp_id in hc[WGS_ID == s]$ID) {
    ta <- r$t[ID == amp_id & tba == "TRUE"]
    ## true chromosome set = amplicon's own breakend chromosomes + bridge partners
    chrs <- unique(c(as.character(ta$seqnames),
                     unlist(strsplit(paste(unique(ta$tb_partner_chr), collapse = ";"), ";"))))
    chrs <- chrs[!is.na(chrs) & nzchar(chrs)]
    chrs <- chrs[order(match(chrs, chr_order))]
    gene <- hc[WGS_ID == s & ID == amp_id]$genes[1]; gh <- if (nzchar(gene)) strsplit(gene, ",")[[1]] else NULL
    msg <- sprintf("%s  %s  chrs=%s  %s", s, amp_id, paste(chrs, collapse = "+"),
                   if (!is.null(gh)) paste(gh, collapse = ",") else "(no oncogene)")
    message(msg)
    ## whole-chromosome windows (1 .. chromosome length) for every linked chromosome
    cr <- do.call(rbind, lapply(chrs, function(c) c(1, as.numeric(cl[[c]]))))
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = r$cnv_data, sv_data = r$vs,
                     genome = "hg38", chromosome = chrs, chromosome_range = cr,
                     events = c("amp", "loh"), genes_to_highlight = gh, save = FALSE, verbose = FALSE),
      error = function(e) { message("  plot failed: ", conditionMessage(e)); NULL })
    if (!is.null(p)) { plots[[length(plots) + 1]] <- p + ggplot2::ggtitle(msg); titles <- c(titles, msg) }
  }
}

## single-canvas montage as one PDF page (2 columns)
if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 2; nrow <- ceiling(length(plots) / ncol)
  montage <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = sprintf("High-confidence TBA amplicons (dual bridge-arm LOH): %d in %d samples",
                      length(plots), uniqueN(hc$WGS_ID)),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16)))
  outpdf <- "validation/output/tba_highconf_montage.pdf"
  ggplot2::ggsave(outpdf, montage, width = 15 * ncol, height = 4.2 * nrow,
                  device = grDevices::cairo_pdf, limitsize = FALSE, bg = "white")
  message("Wrote ", outpdf)
} else message("patchwork not available; montage skipped")
