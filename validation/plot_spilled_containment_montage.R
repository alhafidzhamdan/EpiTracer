## ---------------------------------------------------------------------------
## Montage of the episomal calls that FAILED amplified-junction containment
## (containment > 1.5; from validation/quantify_containment.R), zoomed to each
## amplicon. The window is bounded by the amplicon's own HIGH-VF junctions (not
## the whole chromosome), and low-VF translocations are dropped so passenger arcs
## to non-homologous chromosomes do not clutter the view. Discrete amplified
## islands with diploid gaps read as multi-fragment ecDNA; a continuous elevated
## backbone spanning the region reads as retained / arm-scale chromothripsis.
##
## USAGE:
##   Rscript validation/plot_spilled_containment_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

sel <- fread("validation/output/episomal_containment.tsv")[containment > 1.5]
setorder(sel, -containment)
message("span-spilled episomal calls to plot: ", nrow(sel))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
gene_coord <- data.frame(chr = onc$chr, start = onc$start, end = onc$end,
                         strand = onc$strand, gene = onc$gene)
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
for (s in unique(sel$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  for (i in which(sel$WGS_ID == s)) {
    id <- sel$ID[i]; ch <- sel$chr[i]; pc <- nochr(ch)
    seed <- ecdna_gr[ecdna_gr$ID == id]; if (!length(seed)) next
    amp <- cs[seqnames == ch & copyNumber > 3 * ploidy &
                end >= min(start(seed)[as.character(seqnames(seed)) == ch]) &
                start <= max(end(seed)[as.character(seqnames(seed)) == ch])]
    if (!nrow(amp)) next
    lo <- min(amp$start); hi <- max(amp$end)

    ## drop low-VF translocations (passengers to non-homologous chromosomes)
    touch <- vs[(chrom1 == pc & start1 >= lo & start1 <= hi) | (chrom2 == pc & start2 >= lo & start2 <= hi)]
    maxvf <- if (nrow(touch)) max(touch$VF, na.rm = TRUE) else 0
    thr <- 0.1 * maxvf
    vs_keep <- vs[!(svclass == "TRA" & VF < thr)]

    ## window bounded by the amplicon's kept junctions on the primary chr
    kt <- vs_keep[(chrom1 == pc & start1 >= lo & start1 <= hi) | (chrom2 == pc & start2 >= lo & start2 <= hi)]
    pos <- c(lo, hi, kt$start1[kt$chrom1 == pc], kt$start2[kt$chrom2 == pc])
    wlo <- min(pos); whi <- max(pos); pad <- max(2e6, 0.12 * (whi - wlo))
    cr <- matrix(c(max(1, wlo - pad), min(as.numeric(cl[[ch]]), whi + pad)), nrow = 1)

    gg <- if (nchar(sel$gene[i])) sel$gene[i] else NA
    ttl <- sprintf("%s | %s | %s | fp %.1f Mb, amp-span %.0f Mb, containment %.1f",
                   s, ch, ifelse(is.na(gg) | gg == "", "(no oncogene)", gg),
                   sel$footprint_mb[i], sel$amp_span_mb[i], sel$containment[i])
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = as.data.frame(cs[seqnames == ch]),
                     sv_data = as.data.frame(vs_keep), genome = "hg38",
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
  outpdf <- "validation/output/spilled_containment_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 8 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Containment-failed episomal calls %d-%d of %d (amplicon zoom, low-VF TRAs dropped; discrete islands = ecDNA, continuous backbone = chromothripsis)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = theme(plot.title = element_text(face = "bold", size = 12)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("patchwork not available or no panels; montage skipped")
