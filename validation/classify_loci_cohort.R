## ---------------------------------------------------------------------------
## Per-locus ecDNA/amplicon mechanism classifier for a WGS cohort
##
## Same mechanism logic as classify_loci.R (episomal / translocation-bridge /
## LTA / BFB / complex) but driven by a cohort's own SV BEDPE + allele-specific
## copy-number segments (PURPLE-style) rather than AmpliconArchitect graphs.
## Because allele-specific CN is available, it also reports LOH EVIDENCE:
##   * tp53_17p_loh  - >=75% LOH of 17p up to TP53 (the LTA initiating loss)
##   * subtelo_loh   - subtelomeric LOH on the amplicon's arm (the TB bridge-arm LOH)
##   * tra_to_17p    - a translocation at the locus with its partner in 17p
## so LTA / translocation-bridge calls can be checked against the authoritative
## (Cell 2024 / Nature 2023) LOH signatures.
##
## USAGE (from the package root):
##   Rscript validation/classify_loci_cohort.R example_data/all_348_SV_bedpe.rds \
##            example_data/all_353_CN_segments.rds   # -> output/cohort_loci.tsv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges)
})

TP53_CHR <- "chr17"; TP53_END <- 7.7e6            # 17p terminal .. TP53 (hg38)

classify_cohort <- function(sv, cn,
                            amp_ratio = 3, min_cn = 6, gap = 1e6, min_locus_width = 2e5,
                            boundary_tol_frac = 0.05, foldback_dist = 5e4, foldback_min = 3L,
                            gain_mult = 1.4, min_vf = 5, min_bdup_cn = 2,
                            complexity_max = 3L, hi_junc_frac = 0.25, tra_frac = 0.25,
                            loss_frac = 0.6, loh_thresh = 0.5, centromeres = NULL, centromere_pad = 5e5) {
  sv <- as.data.table(sv); cn <- as.data.table(cn)
  ## harmonise seqnames to "chr" style (CN already has it; SV uses "1")
  pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]
  sv <- sv[VF >= min_vf]

  out <- list()
  for (s in unique(cn$sample)) {
    cs <- cn[sample == s]; vs <- sv[sample == s]
    ploidy <- stats::median(cs$ploidy, na.rm = TRUE); if (!is.finite(ploidy)) ploidy <- 2
    ## sample-level 17p/TP53 LOH (weighted fraction of 17p<->TP53 that is LOH)
    p17 <- cs[seqnames == TP53_CHR & start < TP53_END]
    tp53_loh <- if (nrow(p17)) {
      w <- pmax(1, pmin(p17$end, TP53_END) - p17$start)
      sum(w[p17$minorAlleleCopyNumber < loh_thresh]) / sum(w) >= 0.75
    } else FALSE

    ## amplified loci: focal CN>max(min_cn, amp_ratio*ploidy), merged within gap
    a <- cs[copyNumber > pmax(min_cn, amp_ratio * ploidy)]
    if (!nrow(a)) next
    red <- GenomicRanges::reduce(GRanges(a$seqnames, IRanges(a$start, a$end)),
                                 min.gapwidth = gap)
    red <- red[width(red) >= min_locus_width]
    if (!length(red)) next

    for (i in seq_along(red)) {
      ch <- as.character(seqnames(red[i])); lo <- start(red[i]); up <- end(red[i]); span <- up - lo
      maxcn <- max(cs[seqnames == ch & start <= up & end >= lo]$copyNumber)
      base  <- max(2, round(ploidy))
      intra <- vs[chrom1 == ch & chrom2 == ch &
                  pmin(start1, start2) <= up & pmax(start1, start2) >= lo]
      ## boundary DUP reaching both edges
      btol <- max(2e5, boundary_tol_frac * span)
      d <- intra[svclass == "DUP"]
      bd_cn <- NA_real_
      if (nrow(d)) {
        aa <- pmin(d$start1, d$start2); bb <- pmax(d$start1, d$start2)
        ok <- abs(aa - lo) <= btol & abs(bb - up) <= btol & d$JCN >= min_bdup_cn
        if (any(ok)) bd_cn <- max(d$JCN[ok])
      }
      inv <- intra[svclass %in% c("h2hINV", "t2tINV")]
      fb_near <- if (nrow(inv)) sum(abs(inv$start1 - inv$start2) < foldback_dist) else 0L
      n_hi_junc <- sum(intra$JCN >= hi_junc_frac * maxcn)
      ## diploid flanks
      lf <- cs[seqnames == ch & end <= lo][which.max(end)]$copyNumber
      rf <- cs[seqnames == ch & start >= up][which.min(start)]$copyNumber
      lf <- if (length(lf)) lf else base; rf <- if (length(rf)) rf else base
      diploid_flanks <- lf < gain_mult * base & rf < gain_mult * base
      ## translocations at the locus
      tr <- vs[svclass == "TRA" &
               ((chrom1 == ch & start1 >= lo & start1 <= up) |
                (chrom2 == ch & start2 >= lo & start2 <= up))]
      max_tra_cn <- if (nrow(tr)) max(tr$JCN) else 0
      substantial_tra <- max_tra_cn >= tra_frac * maxcn
      bfb <- fb_near >= foldback_min
      ## TRA partner in 17p (LTA connection)
      tra_to_17p <- nrow(tr[(chrom1 == TP53_CHR & start1 < 24e6) |
                            (chrom2 == TP53_CHR & start2 < 24e6)]) > 0
      ## subtelomeric LOH on this arm (TB bridge-arm signature): LOH in the
      ## telomere-proximal 10 Mb on the amplicon side of the centromere
      cc <- if (!is.null(centromeres)) centromeres[centromeres$chr == ch, ] else NULL
      cen_mid <- if (!is.null(cc) && nrow(cc)) (cc$start + cc$end) / 2 else NA
      arm_end <- if (!is.na(cen_mid) && lo < cen_mid) 0 else max(cs[seqnames == ch]$end, na.rm = TRUE)
      subt <- cs[seqnames == ch & abs((start + end) / 2 - arm_end) < 1e7]
      subtelo_loh <- if (nrow(subt)) mean(subt$minorAlleleCopyNumber < loh_thresh) > 0.5 else FALSE
      cen_adj <- if (!is.null(cc) && nrow(cc)) any(up >= cc$start - centromere_pad & lo <= cc$end + centromere_pad) else FALSE
      arm_loss <- any(cs[seqnames == ch]$copyNumber < loss_frac * base &
                      (cs[seqnames == ch]$end - cs[seqnames == ch]$start) > 1e6)

      mech <- if (bfb && substantial_tra) "LTA"
      else if (is.na(bd_cn)) {
        if (substantial_tra) "translocation-bridge" else if (bfb) "BFB" else "complex"
      } else if (bfb) "BFB"
      else if (n_hi_junc > complexity_max) "complex"
      else if (!diploid_flanks) "complex"
      else if (cen_adj) "complex"
      else if (max_tra_cn > bd_cn) "translocation-bridge"
      else "episomal"

      out[[length(out) + 1L]] <- data.table(
        sample = s, chr = ch, start = lo, end = up, width_mb = round(span / 1e6, 2),
        max_cn = round(maxcn), ploidy = round(ploidy, 2), boundary_dup_cn = round(bd_cn, 1),
        foldbacks = fb_near, n_hi_junc = n_hi_junc, diploid_flanks = diploid_flanks,
        max_tra_cn = round(max_tra_cn, 1), mechanism = mech,
        tp53_17p_loh = tp53_loh, tra_to_17p = tra_to_17p, subtelo_loh = subtelo_loh,
        arm_loss = arm_loss)
    }
  }
  rbindlist(out, fill = TRUE)
}

## ---- run -------------------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  sv <- readRDS(args[1]); cn <- readRDS(args[2])
  ## hg38 centromeres from the bundled karyotype
  ki <- readRDS(system.file("extdata", "chr_info_hg38.rds", package = "EpiTracer"))
  acen <- ki[ki$gieStain == "acen", ]
  cen <- data.frame(chr = tapply(as.character(acen$seqnames), acen$seqnames, `[`, 1),
                    start = tapply(acen$start, acen$seqnames, min),
                    end = tapply(acen$end, acen$seqnames, max))
  cen <- cen[!is.na(cen$chr), ]
  res <- classify_cohort(sv, cn, centromeres = cen)
  dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
  fwrite(res, "validation/output/cohort_loci.tsv", sep = "\t")
  cat("major amplified loci:", nrow(res), "in", length(unique(res$sample)), "samples\n")
  print(res[, .N, by = mechanism][order(-N)])
  cat("\nepisomal loci:", res[mechanism == "episomal", .N], "\n")
  cat("LTA loci with TP53 17p LOH support:", res[mechanism == "LTA" & tp53_17p_loh == TRUE, .N],
      "/", res[mechanism == "LTA", .N], "\n")
  cat("translocation-bridge loci with subtelomeric-arm LOH:",
      res[mechanism == "translocation-bridge" & subtelo_loh == TRUE, .N],
      "/", res[mechanism == "translocation-bridge", .N], "\n")
}
