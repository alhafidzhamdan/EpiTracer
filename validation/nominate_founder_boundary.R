## ---------------------------------------------------------------------------
## Re-nominate the episomal catalogue with the founder-boundary rule.
##
## For every EPISOMAL amplicon, run call_founder_boundary() (VF-stratified test
## for a circularisation DUP) and classify by the founder (max-VF) junction class:
##   DUP        -> founder IS the boundary DUP           => episome CONFIRMED
##   h2h/t2tINV -> founder is a fold-back inversion      => inverted-duplication (REJECT)
##   DEL/TRA    -> founder is neither                    => other (review)
## Reports how many calls flip, and whether the early boundary DUP (when present)
## sits in the founder slot (stratum 1) or is demoted to cluster 1 (stratum 2).
##
## USAGE:
##   Rscript validation/nominate_founder_boundary.R \
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
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)

build <- function(s) {
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
  list(cnv_gr = cnv_gr, bgr = bgr)
}

rows <- list()
for (s in unique(epi_cat$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  r <- tryCatch(call_founder_boundary(NULL, inp$bgr, inp$cnv_gr, cancer_genes_gr), error = function(e) NULL)
  if (is.null(r) || !nrow(as.data.frame(r))) next
  d <- as.data.table(as.data.frame(r))
  d[, early_dup_stratum := as.integer(early_dup_stratum)]
  a <- d[!is.na(founder_class), {
    v <- early_dup_stratum[!is.na(early_dup_stratum)]
    .(founder_class = founder_class[1], founder_vf = as.numeric(founder_vf[1]),
      early_boundary_dup = any(early_boundary_dup == "TRUE"),
      early_dup_stratum = if (length(v)) min(v) else NA_integer_,
      n_strata = max(n_strata), interlocus_tra = any(interlocus_tra == "TRUE"),
      max_tra_vf = { t <- as.numeric(max_tra_vf); if (any(!is.na(t))) max(t, na.rm = TRUE) else NA_real_ })
  }, by = ID]
  a <- a[ID %in% epi_cat[WGS_ID == s]$ID]
  if (nrow(a)) { a[, WGS_ID := s]; rows[[length(rows) + 1]] <- a }
}
res <- rbindlist(rows, fill = TRUE)
res <- merge(res, epi_cat[, .(WGS_ID, ID, gene = genes)], by = c("WGS_ID", "ID"), all.x = TRUE)
## Per-locus verdict on OWN (intrachromosomal) junctions -- TRA excluded from the
## founder call. Simple-excision episome = founder is a boundary DUP, or a DEL
## variant of the circle. INV founder = inverted duplication (not simple excision).
## founder = the copy-gaining junction that formed the amplicon (DUP or inversion);
## deletions are ancestral and excluded upstream. Simple-excision episome = DUP
## founder; inverted duplication = inversion founder.
res[, verdict := fifelse(founder_class == "DUP", "episome (simple excision)",
             fifelse(founder_class %in% c("h2hINV", "t2tINV"), "inverted-dup (not episome)",
             "unresolved (no copy-gaining founder)"))]
res[, episome := verdict == "episome (simple excision)"]
fwrite(res, "validation/output/founder_nomination.tsv", sep = "\t")

cat("\n================ FOUNDER-BOUNDARY RE-NOMINATION (per-locus, TRA excluded) ================\n")
cat("amplicons scored:", nrow(res), "\n\n")
cat("founder junction class (own intrachromosomal junctions):\n"); print(table(res$founder_class, useNA = "ifany"))
cat("\nverdict:\n"); print(table(res$verdict))
cat("\nEPISOMES (simple excision: DUP or DEL founder):", sum(res$episome),
    sprintf("(%.0f%%)\n", 100 * mean(res$episome)))
cat("NOT episomal (inverted-duplication, INV founder):", sum(res$founder_class %in% c("h2hINV","t2tINV")), "\n\n")
cat("EPISOMES by oncogene:\n")
print(sort(table(unlist(strsplit(res[episome == TRUE & gene != ""]$gene, ","))), decreasing = TRUE))
cat("\ninverted-duplication (rejected) by oncogene:\n")
print(sort(table(unlist(strsplit(res[founder_class %in% c("h2hINV","t2tINV") & gene != ""]$gene, ","))), decreasing = TRUE))

## Two-ecDNA-recombination candidates: episomal loci joined to another locus by a TRA
cat("\n--- TWO-ecDNA RECOMBINATION (micronuclear chromothripsis) HYPOTHESIS ---\n")
cat("episomal loci carrying an inter-locus TRA:", sum(res$episome & res$interlocus_tra), "\n")
samp_epi_tra <- res[episome == TRUE & interlocus_tra == TRUE]
## pairs: samples where >=2 episomal loci are TRA-linked (both ends independently episomal)
cand <- samp_epi_tra[, .N, by = WGS_ID][N >= 2]
cat("samples with >=2 TRA-linked EPISOMAL loci (both-episomal recombination):", nrow(cand), "\n")
print(res[episome == TRUE & interlocus_tra == TRUE & WGS_ID %in% cand$WGS_ID,
    .(WGS_ID, ID, gene, founder_class, founder_vf, max_tra_vf)][order(WGS_ID, -founder_vf)])

## plot: founder class composition, and confirmed vs flipped by oncogene
res[, oncogene := vapply(strsplit(gene, ","), function(g) if (length(g) && g[1] != "") g[1] else "(none)", "")]
top <- res[, .N, by = oncogene][order(-N)][seq_len(min(.N, 10))]$oncogene
pl <- res[oncogene %in% top]
pl[, oncogene := factor(oncogene, levels = rev(top))]
p <- ggplot(pl, aes(oncogene, fill = verdict)) + geom_bar() + coord_flip() +
  scale_fill_manual(values = c("episome (simple excision)" = "#2c7fb8",
      "inverted-dup (not episome)" = "#c2181b", "unresolved" = "#f0a848"), name = NULL) +
  labs(title = "Founder-boundary re-nomination (per-locus, TRA excluded)",
       subtitle = sprintf("%d/%d episomes (DUP/DEL founder); %d inverted-duplications; %d episomal loci TRA-linked",
           sum(res$episome), nrow(res), sum(res$founder_class %in% c("h2hINV","t2tINV")),
           sum(res$episome & res$interlocus_tra)),
       x = NULL, y = "amplicons") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))
ggsave("validation/output/founder_nomination.png", p, width = 10, height = 6, dpi = 150, bg = "white")
cat("\nWrote validation/output/founder_nomination.tsv + .png\n")
