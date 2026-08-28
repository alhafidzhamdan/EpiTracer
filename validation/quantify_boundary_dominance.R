## ---------------------------------------------------------------------------
## Boundary-dominance quantification across the episomal catalogue.
##
## For every EPISOMAL amplicon, on its primary (highest-footprint) chromosome:
##   bdup_vf  -- VF of the boundary DUP spanning the footprint (circularisation)
## then count the amplicon's INTERNAL SVs (breakend in the footprint, excluding
## the boundary DUP itself) whose VF EXCEEDS bdup_vf -- i.e. junctions more
## amplified than the circularisation boundary, which the simple-excision model
## says should not exist. Broken down by SV class.
##
## USAGE:
##   Rscript validation/quantify_boundary_dominance.R \
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
  list(cnv_gr = cnv_gr, bgr = bgr, cs = cs, vs = vs)
}

rows <- list()
for (s in unique(epi_cat$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
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
    dup <- touch[svclass == "DUP" & chrom1 == pc & chrom2 == pc]
    span_dup <- dup[pmin(start1, start2) <= lo + 0.2 * fw & pmax(start1, start2) >= hi - 0.2 * fw]
    if (!nrow(span_dup)) { has_boundary <- FALSE; bdup_vf <- NA_real_; bname <- "" }
    else { has_boundary <- TRUE; bi <- which.max(span_dup$VF); bdup_vf <- span_dup$VF[bi]; bname <- span_dup$name[bi] }

    internal <- touch[name != bname]                # exclude the boundary DUP itself
    internal <- internal[!duplicated(name)]
    gt <- if (has_boundary) internal[VF > bdup_vf] else internal[0]
    rows[[length(rows) + 1]] <- data.table(
      WGS_ID = s, ID = id, gene = epi_cat[WGS_ID == s & ID == id]$genes[1], chr = best$ch,
      footprint_mb = round(fw / 1e6, 2), has_boundary = has_boundary, bdup_vf = bdup_vf,
      n_internal = nrow(internal), n_internal_gt = nrow(gt),
      gt_DUP = nrow(gt[svclass == "DUP"]), gt_DEL = nrow(gt[svclass == "DEL"]),
      gt_INV = nrow(gt[svclass %in% c("h2hINV", "t2tINV")]), gt_TRA = nrow(gt[svclass == "TRA"]),
      max_internal_vf = if (nrow(internal)) max(internal$VF) else NA_real_)
  }
}
res <- rbindlist(rows, fill = TRUE)
res[, ratio := round(max_internal_vf / bdup_vf, 2)]
fwrite(res, "validation/output/boundary_dominance.tsv", sep = "\t")

## ---- report ---------------------------------------------------------------
wb <- res[has_boundary == TRUE]
cat("\n============ BOUNDARY-DOMINANCE ACROSS EPISOMAL CATALOGUE ============\n")
cat("episomal amplicons:", nrow(res), " | with an identifiable boundary DUP:", nrow(wb),
    " | no spanning boundary DUP:", sum(!res$has_boundary), "\n\n")
cat("TOTAL internal SVs (breakend in footprint, excl. boundary):", sum(wb$n_internal), "\n")
cat("TOTAL internal SVs with VF > boundary-DUP VF:", sum(wb$n_internal_gt),
    sprintf("(%.1f%% of internal SVs)\n", 100 * sum(wb$n_internal_gt) / max(1, sum(wb$n_internal))))
cat("  of those, by class:  DUP", sum(wb$gt_DUP), " DEL", sum(wb$gt_DEL),
    " INV(fold-back)", sum(wb$gt_INV), " TRA", sum(wb$gt_TRA), "\n\n")
cat("amplicons with >=1 internal SV exceeding the boundary (dominance VIOLATED):",
    sum(wb$n_internal_gt >= 1), sprintf("(%.0f%% of %d)\n", 100 * mean(wb$n_internal_gt >= 1), nrow(wb)))
cat("distribution of internal-SVs-exceeding-boundary per amplicon:\n")
print(table(cut(wb$n_internal_gt, c(-1, 0, 1, 2, 5, 10, 1000),
    labels = c("0", "1", "2", "3-5", "6-10", ">10"))))
cat("\nviolated amplicons by oncogene:\n")
print(sort(table(unlist(strsplit(wb[n_internal_gt >= 1 & gene != ""]$gene, ","))), decreasing = TRUE))
cat("\nworst offenders (most internal SVs above the boundary):\n")
print(head(wb[order(-n_internal_gt), .(WGS_ID, gene, chr, footprint_mb, bdup_vf,
      n_internal, n_internal_gt, gt_INV, ratio)], 12))

## ---- plot -----------------------------------------------------------------
wb[, viol := n_internal_gt >= 1]
p <- ggplot(wb, aes(pmax(n_internal, 0.5), pmax(n_internal_gt, 0.5), colour = viol)) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey70") +
  geom_jitter(width = 0.06, height = 0.06, size = 2.3, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#c2181b"),
                      name = "boundary dominance", labels = c("intact", "violated")) +
  labs(title = "Internal SVs exceeding the boundary-DUP VF, per episomal amplicon",
       subtitle = sprintf("%d of %d amplicons have >=1 internal SV above the boundary; %d internal SVs total exceed it",
                          sum(wb$viol), nrow(wb), sum(wb$n_internal_gt)),
       x = "internal SVs in footprint (log)", y = "internal SVs with VF > boundary DUP (log)") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("validation/output/boundary_dominance.png", p, width = 10, height = 7, dpi = 150, bg = "white")
cat("\nWrote validation/output/boundary_dominance.tsv + .png\n")
