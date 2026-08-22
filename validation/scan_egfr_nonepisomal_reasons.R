## For every EGFR amplicon the caller detected but did NOT classify episomal,
## report WHY it failed, by reading call_simple_excision()'s own internal
## annotation columns plus a recomputation of the flank test. Buckets:
##   no_boundary_dup   - no single DUP spans both amplicon boundaries
##                       (no self-ligation junction -> BFB / complex / other)
##   dup_not_dominant  - a boundary DUP exists but is NOT the highest-VF DUP
##                       (flank test never reached)
##   flank_sliver      - passed boundary+VF; failed flank test ONLY on a
##                       sub-2kb shoulder segment  => candidate missed episome
##   flank_wide        - passed boundary+VF; failed flank test on a >=2kb,
##                       genuinely gained flank    => plausibly non-episomal
## USAGE: Rscript validation/scan_egfr_nonepisomal_reasons.R sv.rds cn.rds
suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges)
})
args <- commandArgs(trailingOnly = TRUE)
sv <- as.data.table(readRDS(args[1])); cn <- as.data.table(readRDS(args[2]))
pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]
onc <- read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                  sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)

calls <- fread("validation/output/cohort_original_caller_global.tsv")
egfr <- calls[grepl("(^|,)EGFR(,|$)", genes)]
epi_samples <- unique(egfr[episomal == TRUE]$WGS_ID)
samples <- sort(setdiff(unique(egfr$WGS_ID), epi_samples))
message(length(samples), " non-episomal EGFR-amplicon samples to scan")
SLIVER <- 2000L  # width (bp) below which a failing flank is deemed a segmentation shoulder

build_inputs <- function(s) {
  cs <- cn[sample == s]; vs <- sv[sample == s]
  if (!nrow(vs) || !nrow(cs)) return(NULL)
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s, copyNumber = cs$copyNumber,
    ploidy = cs$ploidy, majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))
  mk <- function(ch, pos) data.frame(chr = ch, pos = pos, WGS_ID = s, event = ev, svclass = vs$svclass,
    PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN, VF = vs$VF, insLen = 0,
    HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1), mk(vs$chrom2, vs$start2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN, VF = bp$VF,
    insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4))
  ov <- findOverlaps(win, cnv_gr); cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, breakpoints_gr = bgr)
}

chr_baseline_of <- function(cs7, pl) {
  bg <- cs7[copyNumber < 3 * pl]; if (!nrow(bg)) bg <- cs7
  w <- pmax(1, bg$end - bg$start); o <- order(bg$copyNumber)
  max(1, bg$copyNumber[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]])
}

scan_one <- function(s) {
  inp <- build_inputs(s); if (is.null(inp)) return(NULL)
  r <- tryCatch(as.data.table(as.data.frame(
        call_simple_excision(ecdna_gr = NULL, breakpoints_gr = inp$breakpoints_gr,
          cnv_gr = inp$cnv_gr, cancer_genes_gr = cancer_genes_gr, flank_baseline = "chromosome"))),
        error = function(e) NULL)
  if (is.null(r) || !nrow(r)) return(data.table(WGS_ID = s, reason = "caller_error"))
  amp_ids <- unique(r[gene == "EGFR"]$ID)          # amplicon(s) carrying EGFR
  if (!length(amp_ids)) return(data.table(WGS_ID = s, reason = "no_egfr_amp"))
  out <- rbindlist(lapply(amp_ids, function(aid) {
    a <- r[ID == aid & seqnames == "chr7"]
    has_bnd <- any(a$duplication_at_boundary == "TRUE")
    has_vf  <- any(a$duplication_at_boundary_has_highest_VF == "TRUE")
    maxcn   <- round(max(a$PURPLE_CN, na.rm = TRUE))
    row <- data.table(WGS_ID = s, ID = aid, maxCN = maxcn,
                      episomal = any(a$episomal == "TRUE"),
                      excision_scar = any(a$has_excision_scar == "TRUE"),
                      boundary_dup = has_bnd, dup_highest_vf = has_vf,
                      reason = NA_character_, fail_flank = NA_character_,
                      fail_cn = NA_real_, fail_width = NA_integer_, thresh = NA_real_)
    if (!has_bnd) { row$reason <- "no_boundary_dup"; return(row) }
    if (!has_vf)  { row$reason <- "dup_not_dominant"; return(row) }
    ## reached flank test -> recompute which flank failed
    bidx <- a[duplication_at_boundary_has_highest_VF == "TRUE"]$start
    prox <- min(bidx); dist <- max(bidx)
    cs7 <- cn[sample == s & seqnames == "chr7"][order(start)]; pl <- median(cs7$ploidy)
    thr <- 1.4 * chr_baseline_of(cs7, pl); row$thresh <- round(thr, 2)
    fb <- cs7[end <= prox][end == max(end)]; fa <- cs7[start >= dist][start == min(start)]
    fails <- rbindlist(list(
      if (nrow(fb) && fb$copyNumber >= thr) data.table(side="prox", cn=fb$copyNumber, w=fb$end-fb$start),
      if (nrow(fa) && fa$copyNumber >= thr) data.table(side="dist", cn=fa$copyNumber, w=fa$end-fa$start)
    ))
    if (!nrow(fails)) { row$reason <- "flank_pass_other"; return(row) }  # shouldn't happen
    narrow <- fails[w == min(w)][1]
    row$fail_flank <- paste(fails$side, collapse = "+")
    row$fail_cn <- round(narrow$cn, 2); row$fail_width <- as.integer(narrow$w)
    row$reason <- if (all(fails$w < SLIVER)) "flank_sliver" else "flank_wide"
    row
  }), fill = TRUE)
  ## keep the most "episome-like" amplicon per sample (sliver > wide > dup_not_dominant > no_dup)
  rank <- c(flank_sliver = 1, flank_wide = 2, dup_not_dominant = 3, no_boundary_dup = 4,
            flank_pass_other = 5)
  out[order(rank[reason])][1]
}

res <- rbindlist(lapply(samples, scan_one), fill = TRUE)
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "validation/output/egfr_nonepisomal_reasons.tsv", sep = "\t")

cat("\n=== Why the", nrow(res), "EGFR amplicons are NOT episomal ===\n")
print(res[, .N, by = reason][order(-N)])
cat("\n--- flank_sliver: candidate MISSED episomes (boundary DUP + dominant VF; failed on <2kb shoulder) ---\n")
print(res[reason == "flank_sliver",
          .(WGS_ID, maxCN, excision_scar, fail_flank, fail_cn, fail_width, thresh)][order(-maxCN)])
cat("\n--- flank_wide: boundary DUP present but a genuinely gained (>=2kb) flank ---\n")
print(res[reason == "flank_wide",
          .(WGS_ID, maxCN, excision_scar, fail_flank, fail_cn, fail_width, thresh)][order(-maxCN)])
cat("\nWrote validation/output/egfr_nonepisomal_reasons.tsv\n")
