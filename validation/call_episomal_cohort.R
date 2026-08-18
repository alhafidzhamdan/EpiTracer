## ---------------------------------------------------------------------------
## Run the ORIGINAL call_episomal_ecdna() across a WGS cohort
##
## Converts a cohort's SV BEDPE + allele-specific copy number into the GRanges
## the published caller expects (breakpoints_gr / cnv_gr / cancer_genes_gr) and
## runs call_episomal_ecdna(ecdna_gr = NULL, ...) per sample, letting the caller
## auto-detect amplicon seeds from copy number.
##
## `ploidy_mode`:
##   "global"  - use the CN table's own (PURPLE) ploidy in every flank test.
##               This is the caller AS-PUBLISHED. NB: in GBM chr7 is polysomic,
##               so EGFR-episome flanks read as "gained" and the episome is
##               missed (the boundary DUP is found, but the diploid-flank test
##               fails). See GBM39.
##   "per_chr" - set the per-segment ploidy to the local per-chromosome baseline,
##               so the flank test is calibrated to a gained chromosome (recovers
##               the EGFR episomes). Same fix as the per-locus classifier.
##
## USAGE (from the package root):
##   Rscript validation/call_episomal_cohort.R sv_bedpe.rds cn_segments.rds [global|per_chr]
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
sv <- as.data.table(readRDS(args[1])); cn <- as.data.table(readRDS(args[2]))
ploidy_mode <- if (length(args) >= 3) args[3] else "global"

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]

onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)

build_inputs <- function(s) {
  cs <- cn[sample == s]; vs <- sv[sample == s]
  if (!nrow(vs) || !nrow(cs)) return(NULL)
  ploidy_col <- cs$ploidy
  if (ploidy_mode == "per_chr") {                      # local per-chromosome baseline
    base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
    ploidy_col <- base$b[match(cs$seqnames, base$seqnames)]
    ploidy_col[is.na(ploidy_col)] <- 2
  }
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s,
                    copyNumber = cs$copyNumber, ploidy = ploidy_col,
                    majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
                    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))
  mk <- function(ch, pos) data.frame(chr = ch, pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN, VF = vs$VF,
    insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1), mk(vs$chrom2, vs$start2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  ## PURPLE_CN = max CN in a +/-10 kb window (a breakend at an amplicon boundary
  ## must take the amplified side, not the diploid flank)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4))
  ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, breakpoints_gr = bgr)
}

samples <- intersect(unique(sv$sample), unique(cn$sample))
message("Running call_episomal_ecdna (ploidy_mode=", ploidy_mode, ") on ", length(samples), " samples ...")
res <- rbindlist(lapply(samples, function(s) {
  inp <- build_inputs(s); if (is.null(inp)) return(NULL)
  r <- tryCatch(call_episomal_ecdna(ecdna_gr = NULL, breakpoints_gr = inp$breakpoints_gr,
                                    cnv_gr = inp$cnv_gr, cancer_genes_gr = cancer_genes_gr),
                error = function(e) NULL)
  if (is.null(r) || !nrow(as.data.frame(r))) return(NULL)
  d <- as.data.table(as.data.frame(r))
  d[, .(episomal = any(episomal == "TRUE"), excision_scar = any(has_excision_scar == "TRUE"),
        genes = paste(unique(na.omit(gene)), collapse = ",")), by = .(WGS_ID, ID)]
}), fill = TRUE)

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(res, sprintf("validation/output/cohort_original_caller_%s.tsv", ploidy_mode), sep = "\t")
epi <- res[episomal == TRUE]
cat("\namplicons detected:", nrow(res), "in", uniqueN(res$WGS_ID), "samples\n")
cat("EPISOMAL amplicons:", nrow(epi), "in", uniqueN(epi$WGS_ID), "samples",
    "(", sum(epi$excision_scar), "with excision scar )\n")
cat("episomal by oncogene:\n")
print(sort(table(unlist(strsplit(epi[genes != ""]$genes, ","))), decreasing = TRUE))
