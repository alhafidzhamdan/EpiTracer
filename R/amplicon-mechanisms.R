## ---------------------------------------------------------------------------
## Standalone amplicon-mechanism callers: breakage-replication/fusion (BRF),
## micronucleation + chromothripsis, and breakage-fusion-bridge (BFB).
##
## These sit ALONGSIDE call_simple_excision() rather than inside it: each is a
## top-level caller that takes the same inputs (an amplicon catalogue or NULL to
## auto-detect seeds, plus breakpoints_gr / cnv_gr / cancer_genes_gr) and returns
## the annotated breakpoint table with that mechanism's columns. All four share
## the amplicon set-up (annotate_amplicon) and per-chromosome amplified-region
## coordinates (.amplicon_chr_coords).
## ---------------------------------------------------------------------------

#' Per-amplicon breakpoint annotation (shared set-up)
#'
#' Internal helper shared by [call_simple_excision()] and the standalone mechanism
#' callers ([call_brf()], [call_micronucleation()], [call_bfb()]). For one
#' amplicon it builds the structural-variant breakpoint table falling within (or
#' just outside, by `ext`) the amplicon, annotated with oncogene overlap and
#' allele-specific copy number, and returns it together with the sample's full
#' breakpoint set, copy-number segments, ploidy and the amplicon ranges.
#'
#' @param this_amplicon_id Character scalar; one value of `ecdna_gr$ID`.
#' @param ecdna_gr,breakpoints_gr,cnv_gr,cancer_genes_gr,ext See
#'   [call_simple_excision()].
#' @return A list with `bp` (annotated per-breakend data.table for the amplicon),
#'   `sample_bp` (all sample breakends, a GRanges), `cnv` (sample CN GRanges),
#'   `ploidy` (numeric) and `amplicon_gr` (the amplicon ranges); or `NULL` if the
#'   amplicon has no breakpoints in range.
#' @keywords internal
annotate_amplicon <- function(this_amplicon_id, ecdna_gr, breakpoints_gr,
                              cnv_gr, cancer_genes_gr, ext = 1e7) {
  this_amplicon_gr <- ecdna_gr[ecdna_gr$ID %in% this_amplicon_id]
  this_amplicon_gr$gene <- NULL; this_amplicon_gr$group <- NULL; this_amplicon_gr$arm <- NULL

  this_sample <- unique(this_amplicon_gr$WGS_ID)
  this_sample_cnv_gr <- cnv_gr[cnv_gr$sample == this_sample]

  this_sample_breakpoints <- breakpoints_gr[breakpoints_gr$WGS_ID == this_sample] %>% gr2dt()
  if (!"bp_strand" %in% names(this_sample_breakpoints))
    this_sample_breakpoints$bp_strand <- NA_character_
  this_sample_breakpoints <- this_sample_breakpoints %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::select(seqnames, start, end, event, svclass, bp_strand,
                  AF = PURPLE_AF, JCN = PURPLE_JCN, VF, PURPLE_CN, insLen, homLen = HOMLEN) %>%
    to_granges()

  bp <- (this_sample_breakpoints %$% GenomicRanges::trim((this_amplicon_gr + ext))) %>%
    gr2dt() %>%
    dplyr::filter(seqnames %in% as.character(this_amplicon_gr@seqnames)) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(ID != "") %>%
    dplyr::select(-c(strand, width))
  bp <- ((bp %>% to_granges()) %$% cancer_genes_gr) %>% gr2dt() %>% dplyr::select(-c(strand, width))
  bp <- ((bp %>% to_granges()) %$% this_sample_cnv_gr) %>% gr2dt() %>% dplyr::select(-c(strand, width))
  if (nrow(bp) == 0) return(NULL)
  bp$PURPLE_CN <- as.numeric(bp$PURPLE_CN)
  list(bp = data.table::as.data.table(bp), sample_bp = this_sample_breakpoints,
       cnv = this_sample_cnv_gr, ploidy = as.numeric(bp$ploidy[1]),
       amplicon_gr = this_amplicon_gr)
}

## Amplified-region coordinates for one chromosome of an amplicon: the first/last
## CN segment with copyNumber > min_cn_ratio * (segment) ploidy within the
## amplicon's span on that chromosome.
.amplicon_chr_coords <- function(cnv_dt, amplicon_gr, chr, min_cn_ratio = 3) {
  cc <- cnv_dt[seqnames == chr & copyNumber > min_cn_ratio * ploidy]
  if (!nrow(cc)) return(list(has_amp_region = FALSE))
  on <- as.character(GenomeInfoDb::seqnames(amplicon_gr)) == chr
  amp_lo <- min(GenomicRanges::start(amplicon_gr)[on])
  amp_hi <- max(GenomicRanges::start(amplicon_gr)[on] + GenomicRanges::width(amplicon_gr)[on])
  lo <- cc[start >= amp_lo]; hi <- cc[end <= amp_hi]
  if (!nrow(lo) || !nrow(hi)) return(list(has_amp_region = FALSE))
  list(has_amp_region = TRUE, min_amp_coord = min(lo$start), max_amp_coord = max(hi$end))
}

## Engine: seed-detect (if needed), then map a per-chromosome `detector` over
## every amplicon, attaching the detector's columns to the breakpoint table.
.run_amplicon_detector <- function(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                    ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, detector) {
  if (is.null(ecdna_gr)) {
    ecdna_gr <- detect_amplicon_seeds(cnv_gr, min_cn_ratio = min_cn_ratio,
                                      gap = seed_gap, min_width = seed_min_width,
                                      breakpoints = breakpoints_gr)
    if (length(ecdna_gr) == 0L)
      stop("No focal amplicons detected in cnv_gr (need copyNumber > ",
           min_cn_ratio, " x ploidy).", call. = FALSE)
  }
  ids <- unique(ecdna_gr$ID)
  res <- parallel::mclapply(ids, function(id) {
    ctx <- annotate_amplicon(id, ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr, ext)
    if (is.null(ctx) || !nrow(ctx$bp)) return(NULL)
    cnv_dt <- gr2dt(ctx$cnv)
    data.table::rbindlist(lapply(unique(ctx$bp$seqnames), function(ch) {
      this_chr <- ctx$bp[seqnames == ch]
      coords <- .amplicon_chr_coords(cnv_dt, ctx$amplicon_gr, ch, min_cn_ratio)
      cols <- detector(this_chr, coords, ctx, ch, min_cn_ratio)
      for (nm in names(cols)) this_chr[[nm]] <- cols[[nm]]
      this_chr
    }), fill = TRUE)
  }, mc.cores = mc.cores)
  data.table::rbindlist(Filter(Negate(is.null), res), fill = TRUE)
}

## ---- detectors -------------------------------------------------------------

.detect_brf <- function(this_chr, coords, ctx, chr, min_cn_ratio) {
  if (!isTRUE(coords$has_amp_region)) return(list(n_parallel_pairs = 0L, brf = "FALSE"))
  pl <- ctx$ploidy
  pb <- this_chr[start >= coords$min_amp_coord - 2e4 & start <= coords$max_amp_coord + 2e4 &
                   PURPLE_CN > min_cn_ratio * pl]
  ori <- if ("bp_strand" %in% names(pb) && any(!is.na(pb$bp_strand))) pb$bp_strand
         else ifelse(pb$svclass == "h2hINV", "+", ifelse(pb$svclass == "t2tINV", "-", NA_character_))
  n <- nrow(find_parallel_breakpoints(pb$start, ori, pb$event))
  list(n_parallel_pairs = n, brf = if (n >= 1) "TRUE" else "FALSE")
}

.detect_micronucleation <- function(this_chr, coords, ctx, chr, min_cn_ratio) {
  if (!isTRUE(coords$has_amp_region)) return(list(micronucleation = "FALSE"))
  pl <- ctx$ploidy
  tra_amp <- this_chr[svclass == "TRA" & PURPLE_CN > min_cn_ratio * pl]
  mn <- "FALSE"
  if (nrow(tra_amp)) {
    amp_vf <- this_chr[, .(vf = VF[1]), by = event]$vf
    vf75 <- if (length(amp_vf)) stats::quantile(amp_vf, 0.75, names = FALSE, na.rm = TRUE) else Inf
    allbp <- gr2dt(ctx$sample_bp)
    for (ev in unique(tra_amp[VF >= vf75]$event)) {
      partner <- allbp[event == ev & !(seqnames %in% chr)]
      if (nrow(partner) && any(partner$PURPLE_CN > min_cn_ratio * pl)) { mn <- "TRUE"; break }
    }
  }
  list(micronucleation = mn)
}

.detect_bfb <- function(this_chr, coords, ctx, chr, min_cn_ratio, centromeres, chrom_lengths,
                        bfb_min_del_width, bfb_min_del_frac, bfb_loss_max,
                        bfb_min_levels, bfb_min_spread) {
  out <- list(n_foldbacks = 0L, bfb = "FALSE", bfb_anchor = "none")
  if (!isTRUE(coords$has_amp_region)) return(out)
  pl <- ctx$ploidy; min_amp_coord <- coords$min_amp_coord; max_amp_coord <- coords$max_amp_coord
  fb <- this_chr[svclass %in% c("h2hINV", "t2tINV") & PURPLE_CN > min_cn_ratio * pl &
                   start >= min_amp_coord - 2e4 & start <= max_amp_coord + 2e4]
  n_foldbacks <- length(unique(fb$event)); out$n_foldbacks <- n_foldbacks
  if (n_foldbacks < 1 || is.null(centromeres) || !length(centromeres) ||
      is.null(chrom_lengths) || !(chr %in% names(chrom_lengths))) return(out)
  cen <- centromeres[as.character(GenomeInfoDb::seqnames(centromeres)) %in% chr]
  if (!length(cen)) return(out)
  L <- chrom_lengths[[chr]]
  cen_lo <- min(GenomicRanges::start(cen)); cen_hi <- max(GenomicRanges::end(cen))
  ## intrachromosomal: no amplified cross-chromosome TRA partner
  allbp <- gr2dt(ctx$sample_bp)
  tra_ev <- unique(this_chr[svclass == "TRA" &
              start >= min_amp_coord - 2e4 & start <= max_amp_coord + 2e4]$event)
  for (ev in tra_ev) {
    partner <- allbp[event == ev & !(seqnames %in% chr)]
    if (nrow(partner) && any(partner$PURPLE_CN > min_cn_ratio * pl)) return(out)
  }
  arm <- if (max_amp_coord <= cen_lo) "p" else if (min_amp_coord >= cen_hi) "q" else NA
  if (is.na(arm)) return(out)
  chr_cn_bfb <- gr2dt(ctx$cnv); chr_cn_bfb <- chr_cn_bfb[seqnames == chr]
  if (arm == "p") { t_lo <- 0;             t_hi <- min_amp_coord; anchor <- "p-telomere" }
  else            { t_lo <- max_amp_coord; t_hi <- L;             anchor <- "q-telomere" }
  term <- chr_cn_bfb[end >= t_lo & start <= t_hi]
  lost <- term[copyNumber < bfb_loss_max]
  lost_w <- if (nrow(lost)) sum(pmax(0, pmin(lost$end, t_hi) - pmax(lost$start, t_lo))) else 0
  term_w <- max(1, t_hi - t_lo); lost_frac <- lost_w / term_w
  tel_seg <- if (arm == "p") term[which.min(term$start)] else term[which.max(term$end)]
  tel_lost <- nrow(tel_seg) > 0 && all(tel_seg$copyNumber < bfb_loss_max)
  distal_del <- lost_w >= bfb_min_del_width && lost_frac >= bfb_min_del_frac && tel_lost
  ## copy-number staircase (spread fold-backs at >= bfb_min_levels stepped levels)
  staircase <- FALSE
  amp_seg <- chr_cn_bfb[end >= min_amp_coord & start <= max_amp_coord &
                          copyNumber > min_cn_ratio * pl]
  amp_seg <- amp_seg[(pmin(amp_seg$end, max_amp_coord) - pmax(amp_seg$start, min_amp_coord)) >= 5e4]
  if (n_foldbacks >= 2 && nrow(amp_seg) >= 3) {
    fbpos <- fb$start
    spread <- (max(fbpos) - min(fbpos)) / max(1, max_amp_coord - min_amp_coord)
    lv <- sort(amp_seg$copyNumber, decreasing = TRUE); kept <- lv[1]
    for (v in lv[-1]) if (v < tail(kept, 1) * 0.7) kept <- c(kept, v)
    rng <- max(amp_seg$copyNumber) / max(1, min(amp_seg$copyNumber))
    if (spread >= bfb_min_spread && length(kept) >= bfb_min_levels && rng >= 2) staircase <- TRUE
  }
  if (distal_del && staircase) { out$bfb <- "TRUE"; out$bfb_anchor <- anchor }
  out
}

## ---- exported callers ------------------------------------------------------

#' Call breakage-replication/fusion (BRF) amplicons
#'
#' Standalone caller for the BRF annotation: an amplicon carrying *adjacent
#' parallel breakpoints* (two breakends of the same orientation within 20 kb, from
#' distinct junctions), the hallmark of breakage-replication/fusion
#' (Mendez-Dorantes, Zhang, Burns & Pellman, *Nat Genet* 2026;
#' doi:10.1038/s41588-025-02434-5). Runs independently of the episomal call.
#'
#' @inheritParams call_simple_excision
#' @return A [data.table::data.table] of annotated breakpoints with `brf`
#'   (`"TRUE"`/`"FALSE"`) and `n_parallel_pairs`.
#' @seealso [find_parallel_breakpoints()], [call_simple_excision()]
#' @export
call_brf <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                     ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                     mc.cores = 1) {
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, .detect_brf)
}

#' Call micronucleation + chromothripsis amplicons
#'
#' Standalone caller for the micronucleation signature: an amplicon joined to a
#' NON-HOMOLOGOUS chromosome by a high-VF interchromosomal translocation with both
#' breakends amplified ("high VF" = top quartile of the amplicon's own junctions)
#' -- the amplicon incorporates a fragment of another chromosome, the hallmark of
#' assembly inside a micronucleus after chromosome shattering.
#'
#' @inheritParams call_simple_excision
#' @return A [data.table::data.table] of annotated breakpoints with
#'   `micronucleation` (`"TRUE"`/`"FALSE"`).
#' @seealso [call_simple_excision()]
#' @export
call_micronucleation <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                 ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                                 mc.cores = 1) {
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, .detect_micronucleation)
}

#' Call breakage-fusion-bridge (BFB) amplicons
#'
#' Standalone caller for the classical BFB signature (McClintock; genomic
#' signature after Bignell et al. 2007): (i) fold-back inversions, (ii)
#' intrachromosomal (no amplified translocation partner), (iii) terminal
#' amplification with a distal deletion running contiguously to the (absent)
#' telomere, and (iv) a copy-number staircase (fold-backs spread across several
#' stepped levels -- the feature that distinguishes iterative BFB from a single
#' BRF event or a focal ecDNA spike).
#'
#' @inheritParams call_simple_excision
#' @param centromeres,chrom_lengths,bfb_min_del_width,bfb_min_del_frac,bfb_loss_max,bfb_min_levels,bfb_min_spread
#'   See [call_simple_excision()]; `centromeres` and `chrom_lengths` are both
#'   required (the BFB annotation is disabled without them).
#' @return A [data.table::data.table] of annotated breakpoints with `bfb`
#'   (`"TRUE"`/`"FALSE"`), `bfb_anchor` and `n_foldbacks`.
#' @seealso [load_centromeres()], [load_chrom_lengths()], [call_simple_excision()]
#' @export
call_bfb <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                     ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                     centromeres = NULL, chrom_lengths = NULL,
                     bfb_min_del_width = 1e6, bfb_min_del_frac = 0.7, bfb_loss_max = 1.5,
                     bfb_min_levels = 3, bfb_min_spread = 0.3, mc.cores = 1) {
  det <- function(this_chr, coords, ctx, ch, min_cn_ratio)
    .detect_bfb(this_chr, coords, ctx, ch, min_cn_ratio, centromeres, chrom_lengths,
                bfb_min_del_width, bfb_min_del_frac, bfb_loss_max, bfb_min_levels, bfb_min_spread)
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, det)
}
