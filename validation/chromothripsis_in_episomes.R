## ---------------------------------------------------------------------------
## Chromothripsis WITHIN episomal ecDNA across the GBM cohort.
##
## For every amplicon the published caller flags episomal (call_simple_excision,
## per-chromosome flank baseline), score its INTERNAL structural variants for
## chromothripsis (call_chromothripsis; ShatterSeek-style hallmarks over the
## amplified footprint). Reports how many episomes have since shattered
## (ecDNA -> micronucleus -> chromothripsis) and at which oncogenes.
##
## USAGE (from the package root):
##   Rscript validation/chromothripsis_in_episomes.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
centromeres <- load_centromeres("hg38")

## identical input conversion to validation/call_episomal_cohort.R
build_inputs <- function(s) {
  cs <- cn[sample == s]; vs <- sv[sample == s]
  if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  ploidy_col <- base$b[match(cs$seqnames, base$seqnames)]; ploidy_col[is.na(ploidy_col)] <- 2
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s,
                    copyNumber = cs$copyNumber, ploidy = ploidy_col,
                    majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
                    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))
  mk <- function(ch, pos, strand) data.frame(chr = ch, pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = strand, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN,
    VF = vs$VF, insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4))
  ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, breakpoints_gr = bgr)
}

samples <- intersect(unique(sv$sample), unique(cn$sample))
message("Scoring chromothripsis in episomes across ", length(samples), " samples ...")

res_list <- lapply(samples, function(s) {
  inp <- build_inputs(s); if (is.null(inp)) return(NULL)
  a <- list(ecdna_gr = NULL, breakpoints_gr = inp$breakpoints_gr,
            cnv_gr = inp$cnv_gr, cancer_genes_gr = cancer_genes_gr)
  epi <- tryCatch(do.call(call_simple_excision, c(a, list(centromeres = centromeres))), error = function(e) NULL)
  if (is.null(epi) || !nrow(as.data.frame(epi))) return(NULL)
  ct  <- tryCatch(do.call(call_chromothripsis, a), error = function(e) NULL)
  de <- as.data.table(as.data.frame(epi))
  ea <- de[, .(episomal = any(episomal == "TRUE"),
               genes = paste(unique(na.omit(gene)), collapse = ",")), by = .(WGS_ID, ID)]
  if (!is.null(ct) && nrow(ct)) {
    ca <- as.data.table(as.data.frame(ct))[, .(
      chromothripsis = any(chromothripsis == "TRUE"),
      conf = fifelse(any(chromothripsis_conf == "high"), "high",
             fifelse(any(chromothripsis_conf == "low"), "low", "none")),
      n_internal_sv = max(n_internal_sv), n_intrachr_sv = max(n_intrachr_sv),
      sv_type_pval = suppressWarnings(min(sv_type_pval, na.rm = TRUE)),
      cn_oscillations = max(cn_oscillations),
      loh = any(loh_interspersed == "TRUE")), by = .(WGS_ID, ID)]
    ca[!is.finite(sv_type_pval), sv_type_pval := NA_real_]
    ea <- merge(ea, ca, by = c("WGS_ID", "ID"), all.x = TRUE)
  }
  ea
})
res <- rbindlist(Filter(Negate(is.null), res_list), fill = TRUE)

epi <- res[episomal == TRUE]
epi[is.na(chromothripsis), chromothripsis := FALSE]
epi[is.na(conf), conf := "none"]
for (cc in c("n_internal_sv", "n_intrachr_sv", "cn_oscillations")) epi[is.na(get(cc)), (cc) := 0L]

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(epi, "validation/output/episomal_chromothripsis.tsv", sep = "\t")

## ---- report ---------------------------------------------------------------
cat("\n================ CHROMOTHRIPSIS WITHIN EPISOMAL ecDNA ================\n")
cat("EPISOMAL amplicons:", nrow(epi), "in", uniqueN(epi$WGS_ID), "samples\n")
cat("  chromothriptic:", sum(epi$chromothripsis),
    sprintf("(%.1f%%)\n", 100 * mean(epi$chromothripsis)))
cat("  confidence tier:\n"); print(table(factor(epi$conf, c("high", "low", "none"))))
cat("\nchromothriptic (shattered) episomes by oncogene:\n")
print(sort(table(unlist(strsplit(epi[chromothripsis == TRUE & genes != ""]$genes, ","))), decreasing = TRUE))
cat("\nsimple (non-chromothriptic) episomes by oncogene:\n")
print(sort(table(unlist(strsplit(epi[chromothripsis == FALSE & genes != ""]$genes, ","))), decreasing = TRUE))
cat("\ninternal-SV / oscillation summary among episomal amplicons:\n")
print(summary(epi[, .(n_internal_sv, cn_oscillations)]))

## ---- diagnostic plot ------------------------------------------------------
epi[, conf := factor(conf, levels = c("high", "low", "none"))]
epi[, gene1 := vapply(strsplit(genes, ","), function(g) if (length(g) && g[1] != "") g[1] else "(none)", "")]
cols <- c(high = "#c2181b", low = "#f16913", none = "#9e9e9e")

pA <- ggplot(epi, aes(n_internal_sv, cn_oscillations, colour = conf)) +
  geom_hline(yintercept = 3, linetype = 2, colour = "grey55") +
  geom_vline(xintercept = 6, linetype = 2, colour = "grey55") +
  geom_jitter(width = 0.25, height = 0.18, size = 2, alpha = 0.85) +
  scale_colour_manual(values = cols, name = "chromothripsis", drop = FALSE) +
  labs(title = "Each episomal amplicon by its internal-SV burden and CN oscillation",
       subtitle = "dashed lines = calling thresholds (>=6 internal SVs, >=3 CN direction changes)",
       x = "distinct internal SV events", y = "copy-number oscillations (turning points)") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))

topg <- epi[chromothripsis == TRUE & genes != "",
            .(n = uniqueN(paste(WGS_ID, ID))), by = gene1][order(-n)][seq_len(min(.N, 12))]
pB <- ggplot(topg, aes(stats::reorder(gene1, n), n)) +
  geom_col(fill = "#c2181b") + coord_flip() +
  labs(title = "Chromothriptic episomes by oncogene", x = NULL, y = "amplicons") +
  theme_bw(base_size = 12) + theme(plot.background = element_rect(fill = "white", colour = NA))

fig <- if (requireNamespace("patchwork", quietly = TRUE)) patchwork::wrap_plots(pA, pB, widths = c(2, 1)) else pA
ggsave("validation/output/episomal_chromothripsis.png", fig, width = 15, height = 6, dpi = 150, bg = "white")
cat("\nWrote validation/output/episomal_chromothripsis.tsv + .png\n")
