## ---------------------------------------------------------------------------
## Run the ORIGINAL call_simple_excision() across a WGS cohort
##
## Converts a cohort's SV BEDPE + allele-specific copy number into the GRanges
## the published caller expects (breakpoints_gr / cnv_gr / cancer_genes_gr) and
## runs call_simple_excision(ecdna_gr = NULL, ...) per sample, letting the caller
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
centromeres <- load_centromeres("hg38")
chrom_lengths <- load_chrom_lengths("hg38")

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
  mk <- function(ch, pos, strand) data.frame(chr = ch, pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = strand, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN,
    VF = vs$VF, insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
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
message("Running call_simple_excision (ploidy_mode=", ploidy_mode, ") on ", length(samples), " samples ...")
## Each amplicon-formation mechanism now has its OWN standalone caller; run all
## four and join their per-amplicon summaries on WGS_ID + ID.
res <- rbindlist(lapply(samples, function(s) {
  inp <- build_inputs(s); if (is.null(inp)) return(NULL)
  args <- list(ecdna_gr = NULL, breakpoints_gr = inp$breakpoints_gr,
               cnv_gr = inp$cnv_gr, cancer_genes_gr = cancer_genes_gr)
  epi <- tryCatch(do.call(call_simple_excision, c(args, list(centromeres = centromeres))), error = function(e) NULL)
  if (is.null(epi) || !nrow(as.data.frame(epi))) return(NULL)
  brf <- tryCatch(do.call(call_brf, args), error = function(e) NULL)
  mnc <- tryCatch(do.call(call_micronucleation, args), error = function(e) NULL)
  bfb <- tryCatch(do.call(call_bfb, c(args, list(centromeres = centromeres, chrom_lengths = chrom_lengths))), error = function(e) NULL)
  tba <- tryCatch(do.call(call_translocation_bridge_amp, c(args, list(centromeres = centromeres, chrom_lengths = chrom_lengths))), error = function(e) NULL)
  de <- as.data.table(as.data.frame(epi))
  agg <- de[, .(episomal = any(episomal == "TRUE"), episomal_type2 = any(episomal_type2 == "TRUE"),
        excision_scar = any(has_excision_scar == "TRUE"),
        flag_micronucleus = any(flag_micronucleus == "TRUE"),
        flag_chromosomal_bridge = any(flag_chromosomal_bridge == "TRUE"),
        flag_no_boundary_sv = any(flag_no_boundary_sv == "TRUE"),
        flag_inv_at_boundary = any(flag_inv_at_boundary == "TRUE"),
        flag_tra_at_boundary = any(flag_tra_at_boundary == "TRUE"),
        flag_cn_only = any(flag_cn_only == "TRUE"),
        flag_internal_inversion = any(flag_internal_inversion == "TRUE"),
        flag_bridging = any(flag_bridging_amplicon == "TRUE"),
        flag_internal_sv_high_vf = any(flag_internal_sv_high_vf == "TRUE"),
        boundary_homology = suppressWarnings(max(boundary_homology, na.rm = TRUE)),
        junction_homology_class = na.omit(junction_homology_class)[1],
        genes = paste(unique(na.omit(gene)), collapse = ",")), by = .(WGS_ID, ID)]
  agg[!is.finite(boundary_homology), boundary_homology := NA_real_]
  if (!is.null(brf) && nrow(brf)) agg <- merge(agg,
    as.data.table(as.data.frame(brf))[, .(brf = any(brf == "TRUE"), n_parallel_pairs = max(n_parallel_pairs)), by = .(WGS_ID, ID)],
    by = c("WGS_ID", "ID"), all.x = TRUE)
  if (!is.null(mnc) && nrow(mnc)) agg <- merge(agg,
    as.data.table(as.data.frame(mnc))[, .(micronucleation = any(micronucleation == "TRUE")), by = .(WGS_ID, ID)],
    by = c("WGS_ID", "ID"), all.x = TRUE)
  if (!is.null(bfb) && nrow(bfb)) agg <- merge(agg,
    as.data.table(as.data.frame(bfb))[, .(bfb = any(bfb == "TRUE"), n_foldbacks = max(n_foldbacks),
      bfb_anchor = paste(unique(bfb_anchor[bfb_anchor != "none"]), collapse = ",")), by = .(WGS_ID, ID)],
    by = c("WGS_ID", "ID"), all.x = TRUE)
  if (!is.null(tba) && nrow(tba)) agg <- merge(agg,
    as.data.table(as.data.frame(tba))[, .(tba = any(tba == "TRUE"), n_boundary_tra = max(n_boundary_tra),
      tb_partner_chr = paste(unique(tb_partner_chr[tb_partner_chr != ""]), collapse = ";"),
      tb_bridge_arm_loh = any(tb_bridge_arm_loh == "TRUE"),
      tb_partner_arm_loh = any(tb_partner_arm_loh == "TRUE"),
      tb_nonbridge_spared = any(tb_nonbridge_spared == "TRUE"),
      tb_confident = any(tb_confident == "TRUE"),
      tb_high_confidence = any(tb_high_confidence == "TRUE")), by = .(WGS_ID, ID)],
    by = c("WGS_ID", "ID"), all.x = TRUE)
  ## combined initiation-mechanism label (episomal / micronucleation, mutually exclusive)
  agg[, mechanism := fifelse(episomal & micronucleation, "candidate_episomal_with_micronucleation",
                     fifelse(micronucleation, "micronucleation_chromothripsis",
                     fifelse(episomal, "episomal_ecDNA", "unclassified")))]
  agg
}), fill = TRUE)

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(res, sprintf("validation/output/cohort_original_caller_%s.tsv", ploidy_mode), sep = "\t")
epi <- res[episomal == TRUE]
cat("\namplicons detected:", nrow(res), "in", uniqueN(res$WGS_ID), "samples\n")
cat("EPISOMAL amplicons:", nrow(epi), "in", uniqueN(epi$WGS_ID), "samples",
    "(", sum(epi$excision_scar), "with excision scar )\n")
cat("episomal by oncogene:\n")
print(sort(table(unlist(strsplit(epi[genes != ""]$genes, ","))), decreasing = TRUE))

cat("\nepisomal circularisation-junction repair pathway",
    "(boundary-DUP breakpoint homology; Eugen-Olsen et al., NAR 2025):\n")
print(table(epi$junction_homology_class, useNA = "ifany"))
cat("  boundary homology (bp) summary among episomal:\n")
print(summary(epi$boundary_homology))

brf <- res[brf == TRUE]
cat("\nBRF-annotated amplicons (adjacent parallel breakpoints):", nrow(brf),
    "in", uniqueN(brf$WGS_ID), "samples\n")
cat("  of which also episomal:", nrow(res[brf == TRUE & episomal == TRUE]),
    "| non-episomal:", nrow(res[brf == TRUE & episomal == FALSE]), "\n")
cat("BRF by oncogene:\n")
print(sort(table(unlist(strsplit(brf[genes != ""]$genes, ","))), decreasing = TRUE))

cat("\ninitiation mechanism (episomal vs micronucleation mutually exclusive):\n")
print(table(res$mechanism))
mn <- res[micronucleation == TRUE]
cat("\nmicronucleation+chromothripsis amplicons:", nrow(mn), "in", uniqueN(mn$WGS_ID), "samples\n")
cat("  candidate episomal + micronucleation:", nrow(res[mechanism == "candidate_episomal_with_micronucleation"]), "\n")
cat("micronucleation by oncogene:\n")
print(sort(table(unlist(strsplit(mn[genes != ""]$genes, ","))), decreasing = TRUE))

bfb <- res[bfb == TRUE]
cat("\nBFB amplicons (intrachromosomal terminal fold-back staircase with distal deletion to the telomere):",
    nrow(bfb), "in", uniqueN(bfb$WGS_ID), "samples\n")
cat("  anchor:"); print(table(bfb$bfb_anchor))
cat("BFB by oncogene:\n")
print(sort(table(unlist(strsplit(bfb[genes != ""]$genes, ","))), decreasing = TRUE))

tba <- res[tba == TRUE]
cat("\nTBA amplicons (amplified boundary translocation; Lee et al., Nature 2023):",
    nrow(tba), "in", uniqueN(tba$WGS_ID), "samples\n")
cat("  confident TBA (bridge-arm LOH on amplicon and/or partner arm + non-bridge arm spared):",
    nrow(tba[tb_confident == TRUE]), "in", uniqueN(tba[tb_confident == TRUE]$WGS_ID), "samples\n")
cat("  HIGH-confidence TBA (dual-LOH: BOTH amplicon and partner bridge arms + non-bridge spared):",
    nrow(tba[tb_high_confidence == TRUE]), "in", uniqueN(tba[tb_high_confidence == TRUE]$WGS_ID), "samples\n")
cat("    bridge_arm_loh:", nrow(tba[tb_bridge_arm_loh == TRUE]),
    "| partner_arm_loh:", nrow(tba[tb_partner_arm_loh == TRUE]),
    "| nonbridge_spared:", nrow(tba[tb_nonbridge_spared == TRUE]), "\n")
cat("TBA by oncogene:\n")
print(sort(table(unlist(strsplit(tba[genes != ""]$genes, ","))), decreasing = TRUE))
cat("confident-TBA by oncogene:\n")
print(sort(table(unlist(strsplit(tba[tb_confident == TRUE & genes != ""]$genes, ","))), decreasing = TRUE))
