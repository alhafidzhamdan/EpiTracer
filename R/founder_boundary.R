## ---------------------------------------------------------------------------
## Founder-boundary nomination: VF-stratified test for a circularisation DUP.
##
## Episomal excision is defined by a boundary DUP (head-to-tail circularisation
## junction) that BOUNDS the amplicon. Because every copy of the episome carries
## that junction, it is an EARLY / high-copy event -- but not necessarily the
## single highest-VF junction (an internal fold-back inversion can out-copy it,
## as in inverted-duplication amplicons). So we VF-stratify the amplicon's
## junctions exactly as plot_sv_reconstruction() does (founder = max-VF stratum,
## then log-VF k-means clusters ordered by descending median VF) and look for a
## footprint-spanning DUP in the EARLY tier (founder OR the top cluster). We also
## record the class of the founder (max-VF) junction: a founder that is an
## inversion, not a DUP, is the inverted-duplication / fold-back signature and
## argues against simple excision.
## ---------------------------------------------------------------------------

## Per-junction VF stratum (1 = founder / max-VF; 2 = highest non-founder cluster;
## ...). Mirrors assign_strata() + founder isolation in plot_sv_reconstruction().
.assign_vf_strata <- function(vf, max_k = 4L) {
  n <- length(vf)
  if (n <= 1L || length(unique(vf)) == 1L) return(rep(1L, n))
  fmask <- vf == max(vf, na.rm = TRUE)            # founder = max VF (ties = founder)
  rest <- which(!fmask)
  if (!length(rest)) return(rep(1L, n))
  lv <- log(vf[rest]); lv[!is.finite(lv)] <- min(lv[is.finite(lv)], 0)
  uniq_n <- length(unique(round(lv, 8)))
  kk <- min(max_k, uniq_n)
  cl <- if (uniq_n <= 1L) rep(1L, length(rest))
        else if (requireNamespace("Ckmeans.1d.dp", quietly = TRUE))
          as.integer(Ckmeans.1d.dp::Ckmeans.1d.dp(lv, k = c(1L, kk))$cluster)
        else if (kk >= length(rest))          # too few points for k-means: rank into bins
          as.integer(factor(rank(lv, ties.method = "min")))
        else { set.seed(1L); as.integer(stats::kmeans(lv, centers = kk, nstart = 25L)$cluster) }
  ord <- names(sort(tapply(vf[rest], cl, stats::median), decreasing = TRUE))
  rest_stratum <- match(as.character(cl), ord)
  strat <- integer(n); strat[fmask] <- 1L; strat[rest] <- 1L + rest_stratum
  strat
}

.detect_founder_boundary <- function(this_chr, coords, ctx, chr, min_cn_ratio,
                                     max_k = 4L, span_frac = 0.8, early_strata = 2L) {
  out <- list(founder_class = NA_character_, founder_vf = NA_real_,
              early_boundary_dup = "FALSE", early_dup_vf = NA_real_,
              early_dup_stratum = NA_integer_, n_strata = 0L,
              interlocus_tra = "FALSE", max_tra_vf = NA_real_, n_ancestral_del = 0L)
  if (!isTRUE(coords$has_amp_region)) return(out)
  lo <- coords$min_amp_coord; hi <- coords$max_amp_coord; fw <- hi - lo
  ## per-event summary from breakends touching the footprint (both breakends of an
  ## intrachromosomal junction land in this_chr, so p1/p2 give its span)
  bp <- this_chr[start >= lo - 2e4 & start <= hi + 2e4]
  if (!nrow(bp)) return(out)
  ev <- bp[, .(vf = VF[1], svclass = svclass[1], p1 = min(start), p2 = max(start), nbp = .N),
           by = event]
  ev <- ev[is.finite(vf)]
  if (!nrow(ev)) return(out)
  ## Inter-locus translocations join this circle to ANOTHER amplicon; they are a
  ## fusion event, NOT part of this locus's founding excision. Assess this locus on
  ## its own (intrachromosomal) junctions and report the TRA separately -- two
  ## independently-episomal loci linked by a TRA are the two-ecDNA-recombination
  ## (micronuclear chromothripsis) hypothesis.
  tra <- ev[svclass == "TRA"]
  if (nrow(tra)) { out$interlocus_tra <- "TRUE"; out$max_tra_vf <- max(tra$vf) }
  ## Only copy-GAINING junctions can found an amplicon: a DUP (circularisation /
  ## tandem) or an INV (inverted duplication). A DELETION removes sequence, so a
  ## high-VF deletion is an ANCESTRAL event (present on the locus before/independent
  ## of the excision, high-copy only because every copy carries it), NOT a founder.
  ## Report it, but exclude it from the founding-mechanism call.
  out$n_ancestral_del <- nrow(ev[svclass == "DEL"])
  amp <- ev[svclass %in% c("DUP", "h2hINV", "t2tINV")]
  if (!nrow(amp)) return(out)
  amp[, stratum := .assign_vf_strata(vf, max_k = max_k)]
  out$n_strata <- max(amp$stratum)
  fi <- which.max(amp$vf)
  out$founder_class <- amp$svclass[fi]; out$founder_vf <- amp$vf[fi]
  ## early footprint-spanning DUP: DUP, two breakends, spanning >= span_frac of the
  ## footprint, in the founder or top cluster (stratum <= early_strata)
  dup <- amp[svclass == "DUP" & nbp >= 2L & stratum <= early_strata &
              p1 <= lo + (1 - span_frac) * fw & p2 >= hi - (1 - span_frac) * fw]
  if (nrow(dup)) {
    di <- which.max(dup$vf)
    out$early_boundary_dup <- "TRUE"; out$early_dup_vf <- dup$vf[di]
    out$early_dup_stratum <- dup$stratum[di]
  }
  out
}

#' Nominate episomes by a VF-early circularisation DUP
#'
#' Standalone caller that VF-stratifies each amplicon's junctions (as
#' [plot_sv_reconstruction()] does) and flags a footprint-spanning boundary DUP in
#' the EARLY tier -- the founder or the highest-VF cluster ("cluster 1"). The
#' circularisation DUP need not be the single highest-VF junction; an inverted-
#' duplication amplicon's fold-back can carry more copies, so requiring the boundary
#' DUP merely to be *early* rather than *maximal* is the correct test.
#'
#' The founding event of an amplicon must be a copy-GAINING junction: a DUP
#' (circularisation / tandem) or an inversion (inverted duplication). A DELETION
#' removes sequence and therefore cannot found an amplification -- a high-VF
#' deletion is an ANCESTRAL event (high-copy only because every copy carries it),
#' not a founder -- so deletions are excluded from the founding-mechanism call and
#' counted in `n_ancestral_del`. `founder_class` is thus DUP (simple excision) or
#' an inversion (inverted duplication, not simple excision).
#'
#' Each amplicon is assessed on its OWN (intrachromosomal) junctions:
#' inter-locus translocations are excluded from the founding-mechanism call and
#' reported separately (`interlocus_tra`, `max_tra_vf`), because a TRA joins this
#' circle to another amplicon rather than forming it. Two independently-episomal
#' (DUP-founder) loci linked by a TRA are the two-ecDNA-recombination (micronuclear
#' chromothripsis) hypothesis.
#'
#' @inheritParams call_simple_excision
#' @param max_k Maximum number of VF strata (default 4), as in
#'   [plot_sv_reconstruction()].
#' @param span_frac A boundary DUP must span at least this fraction of the
#'   amplicon footprint to count (default 0.8).
#' @param early_strata Highest stratum index still considered "early" (default 2:
#'   the founder and the top cluster).
#' @return A [data.table::data.table] of annotated breakpoints with
#'   `founder_class` (svclass of the max-VF copy-gaining junction: DUP or an
#'   inversion), `founder_vf`, `early_boundary_dup` (`"TRUE"`/`"FALSE"`),
#'   `early_dup_vf`, `early_dup_stratum`, `n_strata`, `interlocus_tra`
#'   (`"TRUE"`/`"FALSE"`; a translocation to another locus), `max_tra_vf` and
#'   `n_ancestral_del` (deletions excluded from the founder call).
#' @seealso [plot_sv_reconstruction()], [call_simple_excision()]
#' @export
call_founder_boundary <- function(ecdna_gr = NULL, breakpoints_gr, cnv_gr, cancer_genes_gr,
                                  ext = 1e7, min_cn_ratio = 3, seed_gap = 1e6, seed_min_width = 1e5,
                                  max_k = 4L, span_frac = 0.8, early_strata = 2L, mc.cores = 1) {
  det <- function(this_chr, coords, ctx, ch, min_cn_ratio)
    .detect_founder_boundary(this_chr, coords, ctx, ch, min_cn_ratio, max_k, span_frac, early_strata)
  .run_amplicon_detector(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                         ext, min_cn_ratio, seed_gap, seed_min_width, mc.cores, det)
}
