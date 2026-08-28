## ---------------------------------------------------------------------------
## Two-locus review of the TRA-founder episomal calls that DO have a boundary DUP
## (from validation/nominate_founder_boundary.R). For each, the founding
## translocation's two ends are drawn side by side on a concatenated axis, so you
## can see whether the partner locus is AMPLIFIED with its own boundary DUP
## (multi-fragment ecDNA -> episomal) or a non-amplified bridge target
## (translocation-bridge amplification -> not simple excision). Low-VF passenger
## TRAs are dropped.
##
## USAGE:
##   Rscript validation/plot_tra_coamp_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

nom <- fread("validation/output/founder_nomination.tsv")[founder_class == "TRA" & early_boundary_dup == TRUE]
setorder(nom, -founder_vf)
message("TRA-with-boundary-DUP calls to review: ", nrow(nom))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
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
cn_at <- function(cs, ch, pos) { seg <- cs[seqnames == ch & start <= pos + 1e4 & end >= pos - 1e4]; if (nrow(seg)) max(seg$copyNumber) else NA_real_ }
base_at <- function(cs, ch) { b <- cs[seqnames == ch]$ploidy[1]; if (is.na(b)) 2 else b }

plots <- list()
for (s in unique(nom$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  for (i in which(nom$WGS_ID == s)) {
    id <- nom$ID[i]
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
    ch1 <- best$ch; pc1 <- nochr(ch1); lo <- best$lo; hi <- best$hi

    ## founding TRA = max-VF translocation with a breakend in the footprint
    tra <- vs[svclass == "TRA" & ((chrom1 == pc1 & start1 >= lo & start1 <= hi) |
                                  (chrom2 == pc1 & start2 >= lo & start2 <= hi))]
    if (!nrow(tra)) next
    ft <- tra[which.max(VF)]
    if (ft$chrom1 == pc1 && ft$start1 >= lo && ft$start1 <= hi) { pchr <- ft$chrom2; ppos <- ft$start2 }
    else { pchr <- ft$chrom1; ppos <- ft$start1 }
    pch2 <- pfx(pchr)
    partner_cn <- cn_at(cs, pch2, ppos); partner_amp <- isTRUE(partner_cn > base_at(cs, pch2) + 1)

    ## windows: primary footprint, and partner region (amplified extent near ppos)
    pad1 <- max(2e6, 0.2 * (hi - lo))
    pamp <- cs[seqnames == pch2 & copyNumber > 3 * ploidy & start <= ppos + 5e6 & end >= ppos - 5e6]
    if (nrow(pamp)) { plo <- min(pamp$start); phi <- max(pamp$end) } else { plo <- ppos - 3e6; phi <- ppos + 3e6 }
    pad2 <- max(2e6, 0.2 * (phi - plo))
    chrs <- c(ch1, pch2)
    cr <- rbind(c(max(1, lo - pad1), min(as.numeric(cl[[ch1]]), hi + pad1)),
                c(max(1, plo - pad2), min(as.numeric(cl[[pch2]]), phi + pad2)))
    ## drop low-VF passenger TRAs
    maxvf <- max(vs[(chrom1 == pc1 & start1 >= lo & start1 <= hi) | (chrom2 == pc1 & start2 >= lo & start2 <= hi)]$VF)
    vs_keep <- vs[!(svclass == "TRA" & VF < 0.1 * maxvf)]

    gg <- if (nchar(nom$gene[i])) nom$gene[i] else NA
    ttl <- sprintf("%s | %s<->%s | %s | TRA founder VF %.0f | partner CN %.1f = %s",
                   s, ch1, pch2, ifelse(is.na(gg) | gg == "", "(no oncogene)", gg), ft$VF,
                   ifelse(is.na(partner_cn), NA, partner_cn),
                   ifelse(partner_amp, "AMPLIFIED (multi-fragment ecDNA?)", "not amplified (bridge?)"))
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = as.data.frame(cs[seqnames %in% chrs]),
                     sv_data = as.data.frame(vs_keep), genome = "hg38",
                     chromosome = chrs, chromosome_range = cr, events = c("amp", "loh"),
                     genes_to_highlight = if (is.na(gg)) NULL else unlist(strsplit(gg, ",")),
                     gene_coord = gene_coord, save = FALSE, verbose = FALSE),
      error = function(e) { message("  panel failed ", s, "/", id, ": ", conditionMessage(e)); NULL })
    if (!is.null(p)) { plots[[length(plots) + 1]] <- p + ggtitle(ttl) +
        theme(plot.title = element_text(size = 8, face = "bold")); message(ttl) }
  }
}
message("panels: ", length(plots))

if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 2; nrow <- 4; per <- ncol * nrow
  outpdf <- "validation/output/tra_coamp_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 13 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("TRA-founder episomes with a boundary DUP %d-%d of %d (partner amplified = multi-fragment ecDNA; not amplified = bridge)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = theme(plot.title = element_text(face = "bold", size = 13)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("no panels")
