## ---------------------------------------------------------------------------
## Focality / boundary-dominance audit of every EPISOMAL call.
##
## Motivated by HMF001167T2 (chr3q), where the simple-excision caller flagged the
## highest-CN sub-window of an ARM-SCALE chromothriptic event as a focal episome.
## For each episomal amplicon we measure, on its primary (highest-footprint)
## chromosome:
##   footprint_mb    -- width of the amplified footprint (CN > 3x ploidy)
##   jspan_mb        -- genomic span of the breakends of every SV touching the
##                      footprint (captures junctions that spill beyond it)
##   arm_frac        -- jspan / length of the chromosome arm it sits on
##   bdup_vf         -- VF of the boundary DUP spanning the footprint (circularisation)
##   max_fb_vf       -- max VF of an internal fold-back inversion in the footprint
## Flags: arm_scale (jspan_mb > jspan_mb_thr OR arm_frac > arm_frac_thr) and
##   foldback_dominated (max_fb_vf > bdup_vf). Either => suspect (likely arm-scale
##   chromothripsis / complex amplification miscalled as episomal).
##
## USAGE:
##   Rscript validation/scan_episomal_focality.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
jspan_mb_thr <- 15      # focal episomes span < a few Mb; 15 Mb is generous
arm_frac_thr <- 0.30    # or > 30% of the chromosome arm

sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))
epi_cat <- fread("validation/output/episomal_chromothripsis.tsv")[episomal == TRUE]
setkey(epi_cat, WGS_ID, ID)

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
cen <- load_centromeres("hg38"); cl <- load_chrom_lengths("hg38")
cen_dt <- as.data.table(as.data.frame(cen))[, .(chr = as.character(seqnames),
             cen_lo = min(start), cen_hi = max(end)), by = seqnames][, seqnames := NULL]

build_cnv_bp <- function(s) {
  cs <- cn[sample == s]; vs <- sv_raw[sample == s]; if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  pcol <- base$b[match(cs$seqnames, base$seqnames)]; pcol[is.na(pcol)] <- 2
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

rows <- list()
for (s in unique(epi_cat$WGS_ID)) {
  inp <- build_cnv_bp(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  cs$ploidy <- base$b[match(cs$seqnames, base$seqnames)]
  for (id in epi_cat[WGS_ID == s]$ID) {
    seed <- ecdna_gr[ecdna_gr$ID == id]; if (!length(seed)) next
    ## primary chromosome = seed chr with the widest amplified footprint
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
    pc <- nochr(best$ch); lo <- best$lo; hi <- best$hi; fw <- hi - lo

    ## SVs with a breakend inside the footprint on the primary chr
    b1 <- vs$chrom1 == pc & vs$start1 >= lo & vs$start1 <= hi
    b2 <- vs$chrom2 == pc & vs$start2 >= lo & vs$start2 <= hi
    touch <- vs[b1 | b2]
    if (!nrow(touch)) next
    ## all breakend positions ON the primary chr among touching SVs (captures spill)
    pos <- c(touch$start1[touch$chrom1 == pc], touch$start2[touch$chrom2 == pc])
    jspan <- max(pos) - min(pos)

    ## boundary DUP spanning the footprint (intrachromosomal, spans >= 80%)
    dup <- touch[svclass == "DUP" & chrom1 == pc & chrom2 == pc]
    span_dup <- dup[pmin(start1, start2) <= lo + 0.2 * fw & pmax(start1, start2) >= hi - 0.2 * fw]
    bdup_vf <- if (nrow(span_dup)) max(span_dup$VF) else NA_real_
    ## internal fold-back inversions in the footprint
    fb <- touch[svclass %in% c("h2hINV", "t2tINV") & chrom1 == pc & chrom2 == pc]
    max_fb_vf <- if (nrow(fb)) max(fb$VF) else NA_real_

    ## arm fraction
    cc <- cen_dt[chr == best$ch]
    armlen <- if (nrow(cc)) { if ((lo + hi) / 2 < cc$cen_lo) cc$cen_lo else as.numeric(cl[[best$ch]]) - cc$cen_hi } else as.numeric(cl[[best$ch]])
    arm_frac <- jspan / max(1, armlen)

    arm_scale <- (jspan / 1e6 > jspan_mb_thr) || (arm_frac > arm_frac_thr)
    fb_dom <- is.finite(max_fb_vf) && (!is.finite(bdup_vf) || max_fb_vf > bdup_vf)
    rows[[length(rows) + 1]] <- data.table(
      WGS_ID = s, ID = id, gene = epi_cat[.(s, id)]$genes, chr = best$ch,
      footprint_mb = round(fw / 1e6, 2), jspan_mb = round(jspan / 1e6, 2),
      arm_frac = round(arm_frac, 2), n_junctions = nrow(touch), n_foldbacks = nrow(fb),
      bdup_vf = bdup_vf, max_fb_vf = max_fb_vf, has_spanning_dup = nrow(span_dup) > 0,
      arm_scale = arm_scale, foldback_dominated = fb_dom, suspect = arm_scale || fb_dom,
      chromothripsis = epi_cat[.(s, id)]$chromothripsis)
  }
}
res <- rbindlist(rows, fill = TRUE)
fwrite(res, "validation/output/episomal_focality_scan.tsv", sep = "\t")

## ---- report ---------------------------------------------------------------
cat("\n================ FOCALITY / BOUNDARY-DOMINANCE AUDIT ================\n")
cat("episomal amplicons scanned:", nrow(res), "\n")
cat("  arm_scale (jspan >", jspan_mb_thr, "Mb OR >", arm_frac_thr, "of arm):",
    sum(res$arm_scale), sprintf("(%.0f%%)\n", 100 * mean(res$arm_scale)))
cat("  foldback_dominated (internal FB VF > boundary-DUP VF):",
    sum(res$foldback_dominated), sprintf("(%.0f%%)\n", 100 * mean(res$foldback_dominated)))
cat("  no spanning boundary DUP at all:", sum(!res$has_spanning_dup), "\n")
cat("  SUSPECT (either flag):", sum(res$suspect), sprintf("(%.0f%%)\n", 100 * mean(res$suspect)))
cat("\nsuspect breakdown vs chromothripsis annotation:\n")
print(table(suspect = res$suspect, chromothripsis = res$chromothripsis))
cat("\nclean focal episomes (not suspect) by oncogene:\n")
print(sort(table(unlist(strsplit(res[suspect == FALSE & gene != ""]$gene, ","))), decreasing = TRUE))
cat("\nSUSPECT calls by oncogene:\n")
print(sort(table(unlist(strsplit(res[suspect == TRUE & gene != ""]$gene, ","))), decreasing = TRUE))
cat("\nworst offenders (largest junction span):\n")
print(head(res[order(-jspan_mb), .(WGS_ID, gene, chr, footprint_mb, jspan_mb, arm_frac,
      n_foldbacks, bdup_vf, max_fb_vf, arm_scale, foldback_dominated)], 15))

## ---- plot -----------------------------------------------------------------
res[, cat := fifelse(arm_scale & foldback_dominated, "both",
             fifelse(arm_scale, "arm-scale",
             fifelse(foldback_dominated, "fold-back dominated", "focal (clean)")))]
res[, fb_ratio := max_fb_vf / bdup_vf]
cols <- c("focal (clean)" = "#2c7fb8", "arm-scale" = "#f16913",
          "fold-back dominated" = "#6a51a3", "both" = "#c2181b")
p <- ggplot(res, aes(pmax(jspan_mb, 0.1), pmax(fb_ratio, 0.05, na.rm = TRUE), colour = cat)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_vline(xintercept = jspan_mb_thr, linetype = 2, colour = "grey55") +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = cols, name = NULL) +
  labs(title = "Focality audit of episomal calls",
       subtitle = sprintf("%d of %d suspect (arm-scale and/or fold-back-dominated)", sum(res$suspect), nrow(res)),
       x = "junction span on primary chromosome (Mb, log)",
       y = "max internal fold-back VF / boundary-DUP VF (log; >1 = FB dominates)") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("validation/output/episomal_focality_scan.png", p, width = 11, height = 7, dpi = 150, bg = "white")
cat("\nWrote validation/output/episomal_focality_scan.tsv + .png\n")
