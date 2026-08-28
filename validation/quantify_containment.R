## ---------------------------------------------------------------------------
## Amplified-junction containment across the episomal catalogue.
##
## Mechanism-sound test of "self-contained circle" vs "chromosomally-retained /
## arm-spread". For each episomal amplicon, on its primary chromosome footprint
## [lo,hi], take the STRUCTURAL junctions -- SVs whose BOTH breakends land in
## ELEVATED copy number (CN > per-chr baseline + 1, so rearranged arm counts but
## true-diploid passengers do not) -- and measure how far they reach:
##   amp_span_mb    -- span of structural-junction breakends on the primary chr
##   containment    -- amp_span_mb / footprint_mb  (~1 = contained circle; >>1 = spill)
##   spill_mb       -- amp_span_mb - footprint_mb  (elevated junctions beyond footprint)
##   n_offchr       -- structural junctions whose partner is ELEVATED on another chr
## A contained amplicon (low containment, no off-chr elevated partners) reads as a
## self-contained circle; a spilled one reads as retained/arm-scale chromothripsis.
##
## USAGE:
##   Rscript validation/quantify_containment.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))
epi_cat <- fread("validation/output/episomal_chromothripsis.tsv")[episomal == TRUE]

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)

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
  basev <- base$b; names(basev) <- as.character(base$seqnames)
  list(cnv_gr = cnv_gr, bgr = bgr, cs = cs, vs = vs, base = basev)
}

## CN at a set of (chr,pos) points = max CN of the overlapping segment (+/-10kb)
cn_at <- function(chrs, pos, cnv_gr) {
  q <- GRanges(pfx(chrs), IRanges(pmax(1, pos - 1e4), pos + 1e4))
  ov <- findOverlaps(q, cnv_gr)
  out <- rep(NA_real_, length(q))
  m <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  out[as.integer(names(m))] <- m; out
}

rows <- list()
for (s in unique(epi_cat$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs; base <- inp$base; cnv_gr <- inp$cnv_gr
  base_at <- function(chrs) { b <- base[pfx(chrs)]; b[is.na(b)] <- 2; as.numeric(b) }
  for (id in epi_cat[WGS_ID == s]$ID) {
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
    pc <- nochr(best$ch); lo <- best$lo; hi <- best$hi; fw <- hi - lo

    touch <- vs[(chrom1 == pc & start1 >= lo & start1 <= hi) | (chrom2 == pc & start2 >= lo & start2 <= hi)]
    if (!nrow(touch)) next
    ## elevation at each breakend
    cn1 <- cn_at(touch$chrom1, touch$start1, cnv_gr); cn2 <- cn_at(touch$chrom2, touch$start2, cnv_gr)
    el1 <- cn1 > base_at(touch$chrom1) + 1; el2 <- cn2 > base_at(touch$chrom2) + 1
    el1[is.na(el1)] <- FALSE; el2[is.na(el2)] <- FALSE
    structural <- touch[el1 & el2]
    ## structural-junction breakend positions on the primary chr
    posp <- c(structural$start1[structural$chrom1 == pc], structural$start2[structural$chrom2 == pc])
    amp_span <- if (length(posp)) max(posp) - min(posp) else fw
    n_offchr <- nrow(structural[chrom1 != chrom2 & (chrom1 != pc | chrom2 != pc)])
    rows[[length(rows) + 1]] <- data.table(
      WGS_ID = s, ID = id, gene = epi_cat[WGS_ID == s & ID == id]$genes[1], chr = best$ch,
      footprint_mb = round(fw / 1e6, 2), amp_span_mb = round(amp_span / 1e6, 2),
      containment = round(amp_span / max(1, fw), 2), spill_mb = round((amp_span - fw) / 1e6, 2),
      n_structural = nrow(structural), n_offchr_amp = n_offchr)
  }
}
res <- rbindlist(rows, fill = TRUE)
res[, spilled := containment > 1.5 | n_offchr_amp > 0]
fwrite(res, "validation/output/episomal_containment.tsv", sep = "\t")

## ---- report ---------------------------------------------------------------
cat("\n============ AMPLIFIED-JUNCTION CONTAINMENT ACROSS EPISOMAL CATALOGUE ============\n")
cat("episomal amplicons:", nrow(res), "\n\n")
cat("containment ratio (amplified-junction span / footprint) distribution:\n")
print(table(cut(res$containment, c(0, 1.1, 1.5, 2, 5, 1e6),
    labels = c("<=1.1 (tight)", "1.1-1.5", "1.5-2", "2-5", ">5"))))
cat("\nSPILLED (containment>1.5 OR off-chr amplified partner):",
    sum(res$spilled), sprintf("(%.0f%%)\n", 100 * mean(res$spilled)))
cat("  with an off-chromosome amplified partner:", sum(res$n_offchr_amp > 0), "\n\n")
cat("Ground-truth checks:\n")
for (w in c("C3N-02788T1", "HMF001167T2", "DO13072T1", "DUMC05T1")) {
  r <- res[WGS_ID == w]
  if (nrow(r)) for (k in seq_len(nrow(r)))
    cat(sprintf("  %-12s %-6s fp %.1f Mb, amp-span %.1f Mb, containment %.2f -> %s\n",
        r$WGS_ID[k], r$gene[k], r$footprint_mb[k], r$amp_span_mb[k], r$containment[k],
        ifelse(r$spilled[k], "SPILLED", "contained")))
}
cat("\nSPILLED calls by oncogene:\n")
print(sort(table(unlist(strsplit(res[spilled == TRUE & gene != ""]$gene, ","))), decreasing = TRUE))
cat("\nworst offenders (largest containment ratio):\n")
print(head(res[order(-containment), .(WGS_ID, gene, chr, footprint_mb, amp_span_mb, containment,
      n_structural, n_offchr_amp)], 12))

## ---- plot -----------------------------------------------------------------
p <- ggplot(res, aes(pmax(footprint_mb, 0.1), pmax(containment, 0.9), colour = spilled)) +
  geom_hline(yintercept = 1.5, linetype = 2, colour = "grey55") +
  geom_jitter(width = 0.02, height = 0, size = 2.3, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#c2181b"),
                      name = NULL, labels = c("contained (circle)", "spilled (retained/arm-scale)")) +
  labs(title = "Amplified-junction containment of episomal calls",
       subtitle = sprintf("%d of %d spilled (amplified junctions reach beyond the footprint)", sum(res$spilled), nrow(res)),
       x = "footprint width (Mb, log)", y = "containment = amplified-junction span / footprint (log)") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("validation/output/episomal_containment.png", p, width = 10, height = 7, dpi = 150, bg = "white")
cat("\nWrote validation/output/episomal_containment.tsv + .png\n")
