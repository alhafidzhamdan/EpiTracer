## ---------------------------------------------------------------------------
## Standalone amplicon-mechanism callers: breakage-replication/fusion (BRF),
## micronucleation + chromothripsis, breakage-fusion-bridge (BFB), and
## translocation-bridge amplification (TBA).
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

.detect_brf <- function(this_chr, coords, ctx, chr, min_cn_ratio,
                        max_dist = 2e4, min_dist = 1, max_indep_p = 0.05,
                        exclude_insertion_adjacency = TRUE) {
  none <- list(n_parallel_pairs = 0L, min_indep_p = NA_real_, brf = "FALSE")
  if (!isTRUE(coords$has_amp_region)) return(none)
  pl <- ctx$ploidy
  pb <- this_chr[start >= coords$min_amp_coord - max_dist &
                   start <= coords$max_amp_coord + max_dist &
                   PURPLE_CN > min_cn_ratio * pl]
  ori <- if ("bp_strand" %in% names(pb) && any(!is.na(pb$bp_strand))) pb$bp_strand
         else ifelse(pb$svclass == "h2hINV", "+", ifelse(pb$svclass == "t2tINV", "-", NA_character_))
  pp <- find_parallel_breakpoints(pb$start, ori, pb$event, max_dist = max_dist,
                                  min_dist = min_dist, max_indep_p = max_indep_p,
                                  exclude_insertion_adjacency = exclude_insertion_adjacency)
  n <- nrow(pp)
  if (!n) return(none)
  list(n_parallel_pairs = n, min_indep_p = min(pp$indep_p),
       brf = "TRUE")
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

## Chromothripsis test over ONE amplicon's internal structural variants
## (ShatterSeek-style hallmarks; Cortes-Ciriano et al., Nat Genet 2020),
## restricted to the amplified footprint [min_amp_coord, max_amp_coord] on this
## chromosome. Intended to be run on amplicons ALREADY called episomal, to flag
## the subset that have subsequently shattered (ecDNA -> micronucleus ->
## chromothripsis). Three hallmarks are scored:
##   (1) prevalence  -- >= `min_sv` distinct SV events with a breakend inside the
##       footprint (a clean simple episome carries only its boundary junction and
##       a handful of internal SVs, so it fails here and is called FALSE);
##   (2) random joins -- the four intrachromosomal junction orientations
##       (DEL +/-, DUP -/+, h2hINV +/+, t2tINV -/-) are ~equally represented, so a
##       chi-squared goodness-of-fit against a uniform expectation FAILS to reject
##       (p >= `join_p`); a mechanism with a dominant orientation (e.g. BFB's
##       fold-backs) is rejected;
##   (3) CN oscillation -- the rounded copy-number profile changes direction at
##       least `min_oscillations` times across the footprint segments (turning
##       points), the ShatterSeek oscillating-copy-number signature.
## `loh_interspersed` (allele-specific CN willing) is reported as supporting
## evidence but is not required. Random fragment joins (2) are required for any
## positive call; confidence is "high" when copy-number oscillation (3) also
## holds, "low" with prevalence + joins alone, else "none".
.detect_chromothripsis <- function(this_chr, coords, ctx, chr, min_cn_ratio,
                                   min_sv = 6L, min_oscillations = 3L, join_p = 0.05) {
  out <- list(n_internal_sv = 0L, n_intrachr_sv = 0L, sv_type_pval = NA_real_,
              cn_oscillations = 0L, loh_interspersed = "FALSE",
              chromothripsis = "FALSE", chromothripsis_conf = "none")
  if (!isTRUE(coords$has_amp_region)) return(out)
  lo <- coords$min_amp_coord; hi <- coords$max_amp_coord
  intern <- this_chr[start >= lo & start <= hi]
  if (!nrow(intern)) return(out)

  ## (1) prevalence of clustered internal SVs
  out$n_internal_sv <- length(unique(intern$event))

  ## (2) randomness of fragment joins: the 4 intrachromosomal orientation classes
  intra <- intern[svclass %in% c("DEL", "DUP", "h2hINV", "t2tINV")]
  intra <- intra[!duplicated(event)]
  out$n_intrachr_sv <- nrow(intra)
  if (nrow(intra) >= 4) {
    tab <- table(factor(intra$svclass, levels = c("DEL", "DUP", "h2hINV", "t2tINV")))
    out$sv_type_pval <- tryCatch(
      suppressWarnings(stats::chisq.test(tab)$p.value),
      error = function(e) NA_real_)
  }

  ## (3) copy-number oscillation across the footprint (turning points in the
  ## rounded CN profile of the amplicon's own segments)
  seg <- gr2dt(ctx$cnv)[seqnames == chr & end >= lo & start <= hi]
  seg <- seg[order(start)]
  if (nrow(seg) >= 3) {
    v <- rle(round(seg$copyNumber))$values
    if (length(v) >= 3) {
      s <- sign(diff(v)); s <- s[s != 0]
      if (length(s) >= 2) out$cn_oscillations <- sum(s[-length(s)] != s[-1])
    }
    ## supporting: interspersed LOH (het / LOH segments alternate)
    if ("minorAlleleCopyNumber" %in% names(seg)) {
      loh_runs <- rle(seg$minorAlleleCopyNumber < 0.5)
      if (length(loh_runs$lengths) >= 3 && any(loh_runs$values) && any(!loh_runs$values))
        out$loh_interspersed <- "TRUE"
    }
  }

  ## decision -- random fragment joins are required for any positive call (the
  ## hallmark separating chromothripsis from ordered, orientation-biased
  ## mechanisms such as BFB); copy-number oscillation upgrades it to high.
  prevalence_ok <- out$n_internal_sv >= min_sv
  joins_ok      <- !is.na(out$sv_type_pval) && out$sv_type_pval >= join_p
  osc_ok        <- out$cn_oscillations >= min_oscillations
  if (prevalence_ok && joins_ok && osc_ok) {
    out$chromothripsis <- "TRUE"; out$chromothripsis_conf <- "high"
  } else if (prevalence_ok && joins_ok) {
    out$chromothripsis <- "TRUE"; out$chromothripsis_conf <- "low"
  }
  out
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

## Translocation-bridge amplification (Lee et al., Nature 2023): an amplicon
## whose boundary is defined by an inter-chromosomal translocation that is itself
## amplified (the "boundary translocation" carried within the amplicon), the
## footprint of a dicentric chromosome bridge formed by translocation.
##
## Loss-of-heterozygosity along one arm of a chromosome, the bridge footprint of
## a translocation-bridge. Returns TRUE when the arm carries LOH covering at least
## `frac` of its length AND a SUBSTANTIAL LOH segment (>= `loh_min_seg` bp) reaches
## within `anchor_tol` bp of a chromosome landmark on that arm -- the telomere OR
## the centromere (the bridge can break on the telomeric or the centromeric side of
## the amplicon). The `anchor_tol` window tolerates a small heterozygous segment
## sitting right at the landmark (e.g. a peri-centromeric sliver) so centromere-
## proximal LOH that begins a megabase or two past the centromere still anchors;
## `loh_min_seg` ignores tiny/point LOH slivers so an isolated sub-megabase LOH
## near an otherwise-heterozygous terminus does NOT anchor. `minorAlleleCopyNumber
## < loh_max` marks an LOH segment.
.arm_bridge_loh <- function(cn_all, chr, arm, cen_lo, cen_hi, L, loh_max, frac,
                            anchor_tol = 3e6, loh_min_seg = 2e5) {
  d <- cn_all[seqnames == chr]
  if (!nrow(d)) return(FALSE)
  if (arm == "p") { lo <- 0; hi <- cen_lo; tel <- 0; cen <- cen_lo }
  else            { lo <- cen_hi; hi <- L;  tel <- L; cen <- cen_hi }
  seg <- d[end >= lo & start <= hi]
  if (!nrow(seg)) return(FALSE)
  arm_w <- max(1, hi - lo)
  loh <- seg[minorAlleleCopyNumber < loh_max]
  if (!nrow(loh)) return(FALSE)
  loh_w <- sum(pmax(0, pmin(loh$end, hi) - pmax(loh$start, lo)))
  ## anchor only on SUBSTANTIAL LOH segments (ignore tiny slivers/points): a
  ## substantial LOH block must reach within anchor_tol of the telomere or centromere
  sub <- loh[(pmin(end, hi) - pmax(start, lo)) >= loh_min_seg]
  anchored <- FALSE
  if (nrow(sub)) {
    gap_tel <- min(pmin(abs(sub$start - tel), abs(sub$end - tel)))
    gap_cen <- min(pmin(abs(sub$start - cen), abs(sub$end - cen)))
    anchored <- gap_tel <= anchor_tol || gap_cen <= anchor_tol
  }
  anchored && (loh_w / arm_w) >= frac
}

.detect_tba <- function(this_chr, coords, ctx, chr, min_cn_ratio, tb_edge_tol,
                        centromeres = NULL, chrom_lengths = NULL, loh_max = 0.5,
                        bridge_loh_min_frac = 0.15, nonbridge_loh_max_frac = 0.1,
                        bridge_anchor_tol = 3e6) {
  out <- list(n_boundary_tra = 0L, tba = "FALSE", tb_partner_chr = "",
              tb_bridge_arm_loh = "NA", tb_partner_arm_loh = "NA",
              tb_nonbridge_spared = "NA", tb_confident = "FALSE",
              tb_high_confidence = "FALSE")
  if (!isTRUE(coords$has_amp_region)) return(out)
  pl <- ctx$ploidy
  min_c <- coords$min_amp_coord; max_c <- coords$max_amp_coord
  ## breakends amplified as part of the amplicon and anchored at an amplified edge
  at_edge <- this_chr[PURPLE_CN > min_cn_ratio * pl &
                        (abs(start - min_c) <= tb_edge_tol | abs(start - max_c) <= tb_edge_tol)]
  if (!nrow(at_edge)) return(out)
  tra_ev <- unique(at_edge[svclass == "TRA"]$event)          # boundary translocation(s)
  if (!length(tra_ev)) return(out)
  allbp <- gr2dt(ctx$sample_bp)
  ## partner (translocated) breakends: the other (non-`chr`) end of each boundary TRA
  partner_bp <- allbp[event %in% tra_ev & !(seqnames %in% chr)]
  partners <- unique(as.character(partner_bp$seqnames))
  out$tba <- "TRUE"
  out$n_boundary_tra <- length(tra_ev)
  out$tb_partner_chr <- paste(sort(partners), collapse = ",")

  ## --- confirmatory footprint (Lee et al., Nature 2023): a confident TB
  ## amplification shows the ASYMMETRIC bridge-arm LOH pattern -- a large LOH
  ## running to a landmark (telomere or centromere) on a "bridge" arm, on the
  ## amplicon's own arm and/or the translocated partner arm, while the amplicon
  ## chromosome's opposite ("non-bridge") arm is spared. The asymmetry is what
  ## distinguishes TB amplification from (symmetric) chromothripsis.
  if (!is.null(centromeres) && length(centromeres) &&
      !is.null(chrom_lengths) && chr %in% names(chrom_lengths)) {
    cendt <- gr2dt(centromeres)
    cn_all <- gr2dt(ctx$cnv); cn_all[, minorAlleleCopyNumber := as.numeric(minorAlleleCopyNumber)]
    arm_of <- function(cc, pos) {
      ce <- cendt[seqnames == cc]; if (!nrow(ce)) return(NA_character_)
      if (pos <= min(ce$start)) "p" else if (pos >= max(ce$end)) "q" else NA_character_
    }
    cen <- cendt[seqnames == chr]
    if (nrow(cen)) {
      cen_lo <- min(cen$start); cen_hi <- max(cen$end); L <- chrom_lengths[[chr]]
      arm <- if (max_c <= cen_lo) "p" else if (min_c >= cen_hi) "q" else NA_character_
      if (!is.na(arm)) {
        bridge_loh <- .arm_bridge_loh(cn_all, chr, arm, cen_lo, cen_hi, L, loh_max, bridge_loh_min_frac, bridge_anchor_tol)
        ## non-bridge (opposite) arm of the amplicon chromosome must be spared
        nb_arm <- if (arm == "p") "q" else "p"
        if (nb_arm == "p") { nlo <- 0; nhi <- cen_lo } else { nlo <- cen_hi; nhi <- L }
        nb <- cn_all[seqnames == chr & end >= nlo & start <= nhi]
        nb_w <- max(1, nhi - nlo)
        nb_loh <- nb[minorAlleleCopyNumber < loh_max]
        nb_loh_w <- if (nrow(nb_loh)) sum(pmax(0, pmin(nb_loh$end, nhi) - pmax(nb_loh$start, nlo))) else 0
        nb_spared <- (nb_loh_w / nb_w) <= nonbridge_loh_max_frac
        ## partner (translocated) bridge arm(s): the arm carrying each partner breakend
        partner_loh <- FALSE
        for (pc in partners) {
          if (!(pc %in% names(chrom_lengths))) next
          pce <- cendt[seqnames == pc]; if (!nrow(pce)) next
          ppos <- stats::median(partner_bp[seqnames == pc]$start)
          parm <- arm_of(pc, ppos); if (is.na(parm)) next
          if (.arm_bridge_loh(cn_all, pc, parm, min(pce$start), max(pce$end),
                              chrom_lengths[[pc]], loh_max, bridge_loh_min_frac, bridge_anchor_tol)) {
            partner_loh <- TRUE; break
          }
        }
        out$tb_bridge_arm_loh   <- if (bridge_loh) "TRUE" else "FALSE"
        out$tb_partner_arm_loh  <- if (partner_loh) "TRUE" else "FALSE"
        out$tb_nonbridge_spared <- if (nb_spared) "TRUE" else "FALSE"
        ## confident: LOH on either bridge arm (amplicon and/or partner) + non-bridge
        ## spared. high-confidence: the DUAL-LOH pattern -- LOH on BOTH the amplicon
        ## and partner bridge arms -- the strongest TB-amplification signature.
        if ((bridge_loh || partner_loh) && nb_spared) out$tb_confident <- "TRUE"
        if (bridge_loh && partner_loh && nb_spared)   out$tb_high_confidence <- "TRUE"
      }
    }
  }
  out
}

## ---- exported callers ------------------------------------------------------

#' Call breakage-replication/fusion (BRF) amplicons
#'
#' Standalone caller for the BRF annotation: an amplicon carrying *adjacent
#' parallel breakpoints*, the hallmark of breakage-replication/fusion (Zhang,
#' Mendez-Dorantes, Burns & Pellman, *Nat Genet* **58**, 88-99, 2026;
#' doi:10.1038/s41588-025-02434-5). Runs independently of the episomal call.
#'
#' A pair qualifies only when it clears all three of the source paper's tests:
#' the breakends are not part of a `+/-` insertion or overlap adjacency, they
#' share an orientation and sit within `max_dist` of one another, and their
#' estimated probability of independent origin --- the distance between them
#' divided by the distance out to the nearest opposite-orientation breakend ---
#' is at most `max_indep_p`. See [find_parallel_breakpoints()].
#'
#' The independence test is what keeps this specific. Without it the caller fires
#' on essentially any junction-dense amplicon, because same-orientation breakends
#' land within 20 kb of each other by chance once a footprint carries tens of
#' junctions; `validation/simulate_trajectories.R` shows the failure mode on
#' simulated chromothriptic amplicons that contain no BRF event by construction.
#'
#' @inheritParams call_simple_excision
#' @param max_dist,min_dist,max_indep_p,exclude_insertion_adjacency Passed to
#'   [find_parallel_breakpoints()]; the defaults follow the source paper.
#' @return A [data.table::data.table] of annotated breakpoints with `brf`
#'   (`"TRUE"`/`"FALSE"`), `n_parallel_pairs` and `min_indep_p` (the most
#'   confident pair's independence bound; `NA` when there is no pair).
#' @seealso [find_parallel_breakpoints()], [call_simple_excision()]
#' @export
call_brf <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                     ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                     max_dist = 2e4, min_dist = 1, max_indep_p = 0.05,
                     exclude_insertion_adjacency = TRUE, mc.cores = 1) {
  det <- function(this_chr, coords, ctx, ch, min_cn_ratio)
    .detect_brf(this_chr, coords, ctx, ch, min_cn_ratio,
                max_dist, min_dist, max_indep_p, exclude_insertion_adjacency)
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, det)
}

#' Call micronucleation + chromothripsis amplicons
#'
#' Standalone caller for the two-ecDNA micronucleation signature: an amplicon
#' joined to **another amplified locus on a non-homologous chromosome** by a
#' high-VF interchromosomal translocation, with both breakends amplified ("high
#' VF" = top quartile of the amplicon's own junctions). Fusing fragments of two
#' different chromosomes into one amplicon is only possible if both were present
#' together, so this is the signature of **two episomal ecDNAs co-encapsulated in
#' a micronucleus and recombined** after chromosome shattering (chromothripsis).
#'
#' Micronucleation need **not** involve a non-homologous chromosome: a single
#' chromosome -- or two homologous copies of one -- can be mis-segregated into a
#' micronucleus, shattered and rejoined, which presents as clustered
#' *intrachromosomal* rearrangements rather than an interchromosomal
#' translocation. This caller keys on the interchromosomal two-ecDNA-fusion
#' flavour; for the intrachromosomal shattering hallmarks (clustered breakpoints,
#' random fragment joins, copy-number oscillation) use [call_chromothripsis()].
#'
#' @inheritParams call_simple_excision
#' @return A [data.table::data.table] of annotated breakpoints with
#'   `micronucleation` (`"TRUE"`/`"FALSE"`).
#' @seealso [call_simple_excision()], [call_chromothripsis()]
#' @export
call_micronucleation <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                 ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                                 mc.cores = 1) {
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, .detect_micronucleation)
}

#' Call chromothripsis within (episomal) amplicons
#'
#' Standalone caller that scores each amplicon's INTERNAL structural variants for
#' chromothripsis, using the ShatterSeek hallmarks (Cortes-Ciriano et al., Nat
#' Genet 2020) restricted to the amplified footprint. It is designed to be run on
#' amplicons already called episomal (pass them as `ecdna_gr`) to flag the subset
#' that have since shattered -- the ecDNA -> micronucleus -> chromothripsis route
#' -- but works on any amplicon catalogue (or `NULL` to auto-detect seeds).
#'
#' A footprint is called chromothriptic when it carries at least `min_sv` distinct
#' internal SV events, its four intrachromosomal junction orientations are close to
#' equally represented (chi-squared goodness-of-fit p >= `join_p`, i.e. random
#' fragment joins), and its rounded copy-number profile changes direction at least
#' `min_oscillations` times across the footprint (oscillating copy number). All
#' three give `chromothripsis_conf = "high"`. Random fragment joins are required
#' for any positive call (they separate chromothripsis from orientation-biased
#' mechanisms such as BFB); prevalence with random joins but weak oscillation gives
#' `"low"`. A clean simple episome (one boundary junction, few internal SVs) fails
#' the prevalence test and is called `"FALSE"`.
#'
#' @inheritParams call_simple_excision
#' @param min_sv Minimum number of distinct internal SV events for the prevalence
#'   hallmark (default 6).
#' @param min_oscillations Minimum copy-number direction changes (turning points)
#'   across the footprint segments for the oscillation hallmark (default 3).
#' @param join_p Significance threshold for the fragment-join randomness test; the
#'   junction-orientation distribution must NOT differ from uniform at this level
#'   (default 0.05).
#' @return A [data.table::data.table] of annotated breakpoints with `chromothripsis`
#'   (`"TRUE"`/`"FALSE"`), `chromothripsis_conf` (`"high"`/`"low"`/`"none"`),
#'   `n_internal_sv`, `n_intrachr_sv`, `sv_type_pval`, `cn_oscillations` and
#'   `loh_interspersed`.
#' @seealso [call_simple_excision()], [call_micronucleation()], [call_bfb()]
#' @export
call_chromothripsis <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                                min_sv = 6L, min_oscillations = 3L, join_p = 0.05, mc.cores = 1) {
  det <- function(this_chr, coords, ctx, ch, min_cn_ratio)
    .detect_chromothripsis(this_chr, coords, ctx, ch, min_cn_ratio,
                           min_sv, min_oscillations, join_p)
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, det)
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

#' Call translocation-bridge amplification (TBA) amplicons
#'
#' Standalone caller for translocation-bridge (TB) amplification (Lee, Kim,
#' ..., Park, *Nature* 2023; doi:10.1038/s41586-023-06057-w): a focal amplicon
#' whose boundary is defined by an **inter-chromosomal translocation that is
#' itself amplified** (the "boundary translocation" carried within the amplicon).
#' In this model an oncogene neighbourhood is translocated in G1, creating a
#' dicentric chromosome that forms a chromosome bridge in mitosis; its breakage
#' and repair amplify the fragment (often as ecDNA), leaving an amplified
#' translocation demarcating the amplicon edge and co-amplifications on the
#' partner chromosome. Runs independently of the other mechanism callers.
#'
#' @inheritParams call_simple_excision
#' @param tb_edge_tol Integer; how close (bp) a translocation breakend must sit to
#'   an amplified amplicon edge to count as a boundary translocation (default
#'   `1e4`).
#' @param centromeres,chrom_lengths Optional [GenomicRanges::GRanges] of
#'   centromere spans (e.g. [load_centromeres()]) and named chromosome-length
#'   vector ([load_chrom_lengths()]). When both are supplied, each TBA amplicon is
#'   additionally scored for the paper's confirmatory footprint (see `tb_confident`
#'   in the return value); without them the confidence columns are `"NA"`.
#' @param loh_max Numeric; a copy-number segment is loss-of-heterozygosity (LOH)
#'   when its minor-allele copy number is below this (default `0.5`, i.e. ~0).
#' @param bridge_loh_min_frac Numeric; minimum fraction of a "bridge" arm that
#'   must be LOH, with an LOH block anchored near a landmark (telomere OR
#'   centromere), for that arm to count as a bridge arm (default `0.15`).
#' @param nonbridge_loh_max_frac Numeric; maximum LOH fraction allowed on the
#'   amplicon chromosome's opposite ("non-bridge") arm for it to count as spared
#'   (default `0.1`).
#' @param bridge_anchor_tol Integer; how close (bp) an LOH block must reach to a
#'   telomere or centromere to anchor the bridge-arm LOH (default `3e6`). The
#'   tolerance lets centromere-proximal LOH that starts a megabase or two past the
#'   centromere (a heterozygous peri-centromeric sliver in between) still anchor.
#' @return A [data.table::data.table] of annotated breakpoints with `tba`
#'   (`"TRUE"`/`"FALSE"`), `n_boundary_tra` (number of distinct boundary
#'   translocations) and `tb_partner_chr` (comma-separated partner chromosome(s)).
#'   When `centromeres` and `chrom_lengths` are supplied it also reports the
#'   confirmatory footprint of Lee et al. (*Nature* 2023): `tb_bridge_arm_loh`
#'   (bridge-arm LOH on the amplicon's own arm, to a telomere or centromere),
#'   `tb_partner_arm_loh` (the same on a translocated partner arm),
#'   `tb_nonbridge_spared` (the amplicon chromosome's opposite arm retains
#'   heterozygosity), `tb_confident` (`"TRUE"` when a bridge arm -- amplicon
#'   and/or partner -- shows LOH AND the non-bridge arm is spared: the asymmetric
#'   pattern distinguishing TB amplification from symmetric chromothripsis) and
#'   `tb_high_confidence` (`"TRUE"` only for the dual-LOH pattern: LOH on BOTH the
#'   amplicon and partner bridge arms with the non-bridge arm spared -- the
#'   strongest signature of a dicentric translocation bridge).
#' @seealso [call_simple_excision()], [call_bfb()], [call_micronucleation()]
#' @export
call_translocation_bridge_amp <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                          ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6,
                                          seed_min_width = 1e5, tb_edge_tol = 1e4,
                                          centromeres = NULL, chrom_lengths = NULL,
                                          loh_max = 0.5, bridge_loh_min_frac = 0.15,
                                          nonbridge_loh_max_frac = 0.1, bridge_anchor_tol = 3e6,
                                          mc.cores = 1) {
  det <- function(this_chr, coords, ctx, ch, min_cn_ratio)
    .detect_tba(this_chr, coords, ctx, ch, min_cn_ratio, tb_edge_tol,
                centromeres, chrom_lengths, loh_max, bridge_loh_min_frac,
                nonbridge_loh_max_frac, bridge_anchor_tol)
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, det)
}
