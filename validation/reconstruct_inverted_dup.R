## ---------------------------------------------------------------------------
## VF-stratified reconstruction of the INVERTED-DUPLICATION calls (founder = a
## fold-back inversion, from validation/nominate_founder_boundary.R). One page per
## amplicon: the founder ("Max VF") panel shows the dominant inversion; the
## boundary DUP, if present, appears demoted in a lower (cluster 1) panel. This is
## the view for deciding inverted-duplication ecDNA vs reject.
##
## USAGE:
##   Rscript validation/reconstruct_inverted_dup.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

inv <- fread("validation/output/founder_nomination.tsv")[founder_class %in% c("h2hINV", "t2tINV")]
setorder(inv, -founder_vf)
message("inverted-duplication calls to reconstruct: ", nrow(inv))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
gene_coord <- data.frame(chr = onc$chr, start = onc$start, end = onc$end, strand = onc$strand, gene = onc$gene)
cl <- load_chrom_lengths("hg38")

build <- function(s) {
  cs <- cn[sample == s]; vs <- sv_raw[sample == s]; if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  pcol <- base$b[match(cs$seqnames, base$seqnames)]; pcol[is.na(pcol)] <- 2; cs$ploidy <- pcol
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
for (s in unique(inv$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  for (i in which(inv$WGS_ID == s)) {
    id <- inv$ID[i]
    seed <- ecdna_gr[ecdna_gr$ID == id]; if (!length(seed)) next
    best <- NULL
    for (ch in unique(as.character(seqnames(seed)))) {
      amp <- cs[seqnames == ch & copyNumber > 3 * ploidy &
                  end >= min(start(seed)[as.character(seqnames(seed)) == ch]) &
                  start <= max(end(seed)[as.character(seqnames(seed)) == ch])]
      if (!nrow(amp)) next
      w <- max(amp$end) - min(amp$start)
      if (is.null(best) || w > best$w) best <- list(ch = ch, lo = min(amp$start), hi = max(amp$end), w = w)
    }
    if (is.null(best)) next
    ch <- best$ch; lo <- best$lo; hi <- best$hi; pad <- max(2e6, 0.15 * (hi - lo))
    cr <- matrix(c(max(1, lo - pad), min(as.numeric(cl[[ch]]), hi + pad)), nrow = 1)
    gg <- if (nchar(inv$gene[i])) inv$gene[i] else NA
    ttl <- sprintf("%s | %s | %s | founder = %s (VF %.0f)   [inverted-duplication]",
                   s, ch, ifelse(is.na(gg) | gg == "", "(no oncogene)", gg),
                   inv$founder_class[i], inv$founder_vf[i])
    p <- tryCatch(
      plot_sv_reconstruction(sample = s, cnv_data = as.data.frame(cs[seqnames == ch]),
        sv_data = as.data.frame(vs), genome = "hg38", chromosome = ch, chromosome_range = cr,
        events = c("amp", "loh"), genes_to_highlight = if (is.na(gg)) NULL else unlist(strsplit(gg, ",")),
        gene_coord = gene_coord, vf_col = "VF", k = "auto", isolate_founder = TRUE,
        founder_offscale = TRUE, save = FALSE, verbose = FALSE),
      error = function(e) { message("  failed ", s, "/", id, ": ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      plots[[length(plots) + 1]] <- p + patchwork::plot_annotation(
        title = ttl, theme = theme(plot.title = element_text(face = "bold", size = 12)))
      message(ttl)
    }
  }
}
message("pages: ", length(plots))

if (length(plots)) {
  outpdf <- "validation/output/inverted_dup_reconstructions.pdf"
  grDevices::cairo_pdf(outpdf, width = 13, height = 11, onefile = TRUE)
  for (p in plots) print(p)
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("no pages")
