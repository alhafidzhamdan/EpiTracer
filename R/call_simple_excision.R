#' Classify a single ecDNA amplicon as episomal or not
#'
#' Internal worker called once per amplicon by [call_simple_excision()]. It
#' annotates the structural variant breakpoints falling within (or just
#' outside) an ecDNA amplicon, then applies the episome heuristic:
#' \enumerate{
#'   \item find duplication (DUP) breakpoints that sit at the amplicon
#'     boundaries and are themselves amplified;
#'   \item require the boundary DUP to carry the highest variant fraction (VF)
#'     of any DUP in the amplicon;
#'   \item require the chromosomal segments immediately flanking the boundaries
#'     to be non-gained (consistent with a circular episome excised from an
#'     otherwise diploid region);
#'   \item flag a shared flanking deletion as a candidate excision scar.
#' }
#'
#' @param this_amplicon_id Character scalar; a single value of `ecdna_gr$ID`.
#' @param ecdna_gr,breakpoints_gr,cnv_gr,cancer_genes_gr See
#'   [call_simple_excision()].
#' @param ext Integer; bp to extend the amplicon by when searching for
#'   boundary SVs.
#' @param verbose Logical; print progress/diagnostics.
#'
#' @return A [data.table::data.table] of annotated breakpoints for the
#'   amplicon (may have zero rows), or an empty data.table if the amplicon has
#'   no breakpoints in range.
#' @keywords internal
#' @importFrom data.table data.table rbindlist :=
#' @importFrom GenomicRanges trim
#' @importFrom dplyr filter select arrange rename mutate
classify_amplicon_episomal <- function(this_amplicon_id,
                                       ecdna_gr,
                                       breakpoints_gr,
                                       cnv_gr,
                                       cancer_genes_gr,
                                       ext = 1e7,
                                       flank_baseline = "chromosome",
                                       gain_ratio = 1.4,
                                       min_cn_ratio = 3,
                                       min_flank_width = 2000,
                                       bridge_gap = 1e6,
                                       founder_jcn = 30,
                                       centromeres = NULL,
                                       mh_min_homology = 2,
                                       hr_min_homology = 14,
                                       verbose = FALSE) {

  this_amplicon_gr <- ecdna_gr[ecdna_gr$ID %in% this_amplicon_id]

  ## Drop annotation columns that are re-derived below (if present):
  this_amplicon_gr$gene  <- NULL
  this_amplicon_gr$group <- NULL
  this_amplicon_gr$arm   <- NULL

  this_sample <- unique(this_amplicon_gr$WGS_ID)
  ## Base GRanges subsetting (avoids requiring plyranges for filter.GRanges):
  this_sample_cnv_gr <- cnv_gr[cnv_gr$sample == this_sample]

  this_sample_breakpoints <- breakpoints_gr[breakpoints_gr$WGS_ID == this_sample] %>%
    gr2dt()
  ## Optional per-breakend strand (for adjacent-parallel-breakpoint / BRF
  ## detection); absent in AmpliconArchitect-style inputs, so default to NA.
  if (!"bp_strand" %in% names(this_sample_breakpoints))
    this_sample_breakpoints$bp_strand <- NA_character_
  this_sample_breakpoints <- this_sample_breakpoints %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::select(seqnames, start, end, event, svclass, bp_strand,
                  AF = PURPLE_AF, JCN = PURPLE_JCN, VF, PURPLE_CN,
                  insLen, homLen = HOMLEN) %>%
    to_granges()

  ## Breakpoints at (or just outside) the amplicon boundaries:
  this_sample_breakpoints_ecdna_annotated <-
    (this_sample_breakpoints %$% GenomicRanges::trim((this_amplicon_gr + ext))) %>%
    gr2dt() %>%
    dplyr::filter(seqnames %in% as.character(this_amplicon_gr@seqnames)) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(ID != "") %>%
    dplyr::select(-c(strand, width))

  ## Add oncogene info:
  this_sample_breakpoints_ecdna_annotated <-
    ((this_sample_breakpoints_ecdna_annotated %>% to_granges()) %$% cancer_genes_gr) %>%
    gr2dt() %>%
    dplyr::select(-c(strand, width))

  ## Annotate with minor and major allele info:
  this_sample_breakpoints_ecdna_annotated <-
    ((this_sample_breakpoints_ecdna_annotated %>% to_granges()) %$% this_sample_cnv_gr) %>%
    gr2dt() %>%
    dplyr::select(-c(strand, width))

  if (nrow(this_sample_breakpoints_ecdna_annotated) > 0) {

    this_sample_breakpoints_ecdna_annotated$PURPLE_CN <-
      as.numeric(this_sample_breakpoints_ecdna_annotated$PURPLE_CN)
    this_sample_ploidy <- as.numeric(this_sample_breakpoints_ecdna_annotated$ploidy[1])

    ## Prep columns:
    this_sample_breakpoints_ecdna_annotated$duplication_at_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$duplication_at_boundary_has_highest_VF <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$episome_region <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$deletion_flanking_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$before_prox_boundary_not_gained <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$after_dist_boundary_not_gained <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_internal_sv_high_vf <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_internal_inversion <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_bridging_amplicon <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_micronucleus <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_chromosomal_bridge <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_no_boundary_sv <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_inv_at_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_tra_at_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$flag_tra_recombination <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$episomal_type2 <- "FALSE"
    ## Amplicon called by copy number alone (no supporting junction inside it),
    ## propagated from the SV-aware seed detector.
    this_sample_breakpoints_ecdna_annotated$flag_cn_only <-
      if (!is.null(this_amplicon_gr$cn_only) && length(this_amplicon_gr$cn_only) &&
          isTRUE(this_amplicon_gr$cn_only[1])) "TRUE" else "FALSE"
    this_sample_breakpoints_ecdna_annotated$episomal <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$has_excision_scar <- "FALSE"
    ## Circularisation-junction microhomology signature (Eugen-Olsen et al.,
    ## Nucleic Acids Res 2025, doi:10.1093/nar/gkaf122): the breakpoint homology
    ## length of the boundary DUP and the excision repair pathway it implies.
    this_sample_breakpoints_ecdna_annotated$boundary_homology <- NA_real_
    this_sample_breakpoints_ecdna_annotated$junction_homology_class <- NA_character_
    ## NOTE: the breakage-replication/fusion (BRF), micronucleation and
    ## breakage-fusion-bridge (BFB) annotations are computed by the standalone
    ## callers [call_brf()], [call_micronucleation()] and [call_bfb()], not here.

    ## The same amplicon can span >1 chromosome, so classify per chromosome:
    unique_chrs <- unique(this_sample_breakpoints_ecdna_annotated$seqnames)

    this_sample_breakpoints_ecdna_annotated <- lapply(seq_along(unique_chrs), function(i) {

      this_chr <- this_sample_breakpoints_ecdna_annotated %>%
        dplyr::filter(seqnames %in% unique_chrs[i])

      ## Min/max coords of amplified regions for this chromosome:
      has_amp_region <- nrow(this_sample_cnv_gr %>% gr2dt() %>%
                               dplyr::filter(seqnames %in% unique_chrs[i]) %>%
                               dplyr::filter(copyNumber > 3 * ploidy)) > 0

      if (has_amp_region) {

        min_amp_coord <- this_sample_cnv_gr %>% gr2dt() %>%
          dplyr::filter(seqnames %in% unique_chrs[i]) %>%
          dplyr::filter(copyNumber > 3 * ploidy) %>%
          dplyr::filter(start >= min(this_amplicon_gr[as.character(GenomeInfoDb::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@start)) %>%
          dplyr::filter(start == min(start)) %>% .$start

        max_amp_coord <- this_sample_cnv_gr %>% gr2dt() %>%
          dplyr::filter(seqnames %in% unique_chrs[i]) %>%
          dplyr::filter(copyNumber > 3 * ploidy) %>%
          dplyr::filter(end <= max(this_amplicon_gr[as.character(GenomeInfoDb::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@start +
                                     this_amplicon_gr[as.character(GenomeInfoDb::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@width)) %>%
          dplyr::filter(end == max(end)) %>% .$end

        ## No-boundary-SV flag: neither amplified edge (min/max_amp_coord) has a
        ## structural-variant breakend near it (within 10 kb), so the amplicon
        ## boundaries are pure copy-number transitions with no junction support --
        ## e.g. a focal box with only internal SVs (E28T1 EGFR). Informational.
        edge_bp <- gr2dt(this_sample_breakpoints)
        edge_bp <- edge_bp[seqnames %in% unique_chrs[i]]
        near_edge <- nrow(edge_bp) > 0 &&
          (any(abs(edge_bp$start - min_amp_coord) <= 1e4) ||
           any(abs(edge_bp$start - max_amp_coord) <= 1e4))
        if (!near_edge) this_chr$flag_no_boundary_sv <- "TRUE"

        ## Boundary INV / TRA exclusion: if an amplified edge is DEFINED by an
        ## inversion or translocation -- i.e. an INV/TRA breakend at the edge
        ## (within 10 kb, at amplicon-level copy number) whose junction copy number
        ## exceeds any duplication at that same edge -- the boundary is a fold-back
        ## / translocation junction, not a circularisation, so the amplicon is
        ## excluded. Requiring the INV/TRA to out-JCN the boundary DUP avoids
        ## flagging a genuine boundary-DUP episome that merely carries an
        ## incidental INV/TRA breakend nearby.
        for (edge_pos in c(min_amp_coord, max_amp_coord)) {
          e <- edge_bp[edge_bp$PURPLE_CN > min_cn_ratio * this_sample_ploidy &
                         abs(edge_bp$start - edge_pos) <= 1e4]
          if (!nrow(e)) next
          dup_jcn <- max(c(0, e$JCN[e$svclass == "DUP"]))
          inv_jcn <- max(c(0, e$JCN[e$svclass %in% c("h2hINV", "t2tINV")]))
          tra_jcn <- max(c(0, e$JCN[e$svclass == "TRA"]))
          if (inv_jcn > dup_jcn) this_chr$flag_inv_at_boundary <- "TRUE"
          if (tra_jcn > dup_jcn) this_chr$flag_tra_at_boundary <- "TRUE"
        }

        ## Two-ecDNA recombination flag: an amplified inter-chromosomal
        ## translocation joining this amplicon to ANOTHER amplified locus (both
        ## breakends at amplicon-level copy number) is the signature of two
        ## episomal ecDNAs recombining inside a micronucleus (chromothripsis). It
        ## does NOT disqualify this locus -- each locus is assessed on its own
        ## boundary DUP -- but is reported for the recombination hypothesis.
        allbp_tra <- gr2dt(this_sample_breakpoints)
        allbp_tra[, PURPLE_CN := as.numeric(PURPLE_CN)]
        tra_ev <- allbp_tra[svclass == "TRA", .(
          in_amp = any(seqnames %in% unique_chrs[i] & start >= min_amp_coord & start <= max_amp_coord),
          partner_amp = any(!(seqnames %in% unique_chrs[i]) &
                              PURPLE_CN > min_cn_ratio * this_sample_ploidy)), by = event]
        if (nrow(tra_ev[in_amp == TRUE & partner_amp == TRUE]))
          this_chr$flag_tra_recombination <- "TRUE"


        ## Chromosomal-bridge exclusion: if the amplified region's edge terminates
        ## WITHIN a centromere (rather than traversing across it), it is a
        ## dicentric / chromosomal-bridge event, not a self-contained episome.
        ## Flag it and exclude it from the episomal call.
        if (!is.null(centromeres) && length(centromeres)) {
          cen <- centromeres[as.character(GenomeInfoDb::seqnames(centromeres)) %in% unique_chrs[i]]
          if (length(cen)) {
            cen_lo <- min(GenomicRanges::start(cen)); cen_hi <- max(GenomicRanges::end(cen))
            edge_in_cen <- (min_amp_coord >= cen_lo & min_amp_coord <= cen_hi) ||
                           (max_amp_coord >= cen_lo & max_amp_coord <= cen_hi)
            traverses <- (min_amp_coord < cen_lo & max_amp_coord > cen_hi)
            if (edge_in_cen && !traverses) this_chr$flag_chromosomal_bridge <- "TRUE"
          }
        }

        ## Find DUP at boundary; ensure DUP has highest VF; check flanking CN;
        ## check for excision scar.

        ## A DUP counts as a boundary (circularisation) junction only if BOTH of
        ## its breakends are amplified (PURPLE_CN > min_cn_ratio * ploidy) -- i.e.
        ## both ends lie inside the amplicon. A DUP with one breakend in
        ## non-amplified sequence reaches OUT of the amplicon (a passenger tandem
        ## duplication, e.g. DO12952T1's JCN-1 DUP whose far end sits in diploid
        ## 56.7 Mb) and must not be mistaken for the boundary DUP at the edge. This
        ## is the same "both breakends amplified" rule used for SV-aware seed
        ## linking.
        dup_thr <- 3 * this_sample_ploidy
        intra_dup_events <- {
          dd <- this_chr %>% dplyr::filter(svclass == "DUP")
          if (!nrow(dd)) character(0) else {
            agg <- stats::aggregate(PURPLE_CN ~ event, data = dd,
                                    FUN = function(x) as.numeric(all(x > dup_thr) & length(x) >= 2))
            agg$event[agg$PURPLE_CN == 1]
          }
        }

        ## needs to be >1 DUP: a single DUP on a separate chr is unlikely to
        ## be involved in amplicon generation (outside of ecDNA)
        has_dup <- nrow(this_chr %>%
                          dplyr::filter(svclass == "DUP") %>%
                          dplyr::filter(event %in% intra_dup_events)) > 1

        if (has_dup) {

          ## Proximal / distal border DUPs must be ANCHORED at the amplicon edge
          ## (breakend within 10 kb of min_amp_coord / max_amp_coord), not merely
          ## somewhere past it. The breakpoint set is drawn from a +/- ext window
          ## around the amplicon, so a neighbouring amplified block (e.g. the
          ## 161/164 Mb chr6 blocks near DO12370T1's focal 155-156.8 Mb episome)
          ## contributes amplified DUP breakends far beyond the edge; an unbounded
          ## "past the edge" test would pick one of those as the border and split a
          ## clean single-boundary-DUP episome into a spurious "type 2".
          ## Proximal border:
          has_prox_border <- nrow(this_chr %>%
                                    dplyr::filter(svclass == "DUP") %>%
                                    dplyr::filter(event %in% intra_dup_events) %>%
                                    dplyr::filter(abs(start - min_amp_coord) <= 10000)) > 0
          if (has_prox_border) {
            prox_border <- this_chr %>%
              dplyr::filter(svclass == "DUP") %>%
              dplyr::filter(event %in% intra_dup_events) %>%
              dplyr::filter(abs(start - min_amp_coord) <= 10000) %>%
              dplyr::filter(start == min(start)) %>% .$event
          }

          ## Distal border:
          has_dist_border <- nrow(this_chr %>%
                                    dplyr::filter(svclass == "DUP") %>%
                                    dplyr::filter(event %in% intra_dup_events) %>%
                                    dplyr::filter(abs(end - max_amp_coord) <= 10000)) > 0
          if (has_dist_border) {
            dist_border <- this_chr %>%
              dplyr::filter(svclass == "DUP") %>%
              dplyr::filter(event %in% intra_dup_events) %>%
              dplyr::filter(abs(end - max_amp_coord) <= 10000) %>%
              dplyr::filter(start == max(start)) %>% .$event
          }

          if (has_prox_border & has_dist_border) {
            if (length(prox_border) == 1 & length(dist_border) == 1) {
              if (prox_border == dist_border) {
                ## The SAME DUP spans both boundaries -- proximal and distal
                ## breakends are the same SV. This is a connected circularisation
                ## junction: straightforward "type 1" episome.
                this_chr[event %in% dist_border]$duplication_at_boundary <- "TRUE"
              } else {
                ## Two DIFFERENT DUPs bound the amplicon (one reaches the proximal
                ## edge, another the distal edge) but they are NOT the same SV, so
                ## the boundary is not a single connected junction (the DUPs
                ## interleave / are unconnected, e.g. C3N-01334T1). Not called
                ## straightforward episomal -- flagged "episomal (type 2)" instead.
                this_chr$episomal_type2 <- "TRUE"
              }
            }
          }

          ## Any amplicon that carries a boundary DUP is taken through the flank
          ## test, regardless of whether an internal SV (DUP, inversion, deletion)
          ## out-VFs it -- internal rearrangements do not preclude a genuine
          ## episome, so the boundary DUP no longer has to be the highest-VF DUP.
          boundary_index <- which(this_chr$duplication_at_boundary == "TRUE")

          if (length(boundary_index) > 0) {

            ## max VF among DUPs, used only to ANNOTATE whether a boundary DUP
            ## carries the highest DUP VF (informational; no longer a gate).
            max_vf <- max(this_chr %>% dplyr::filter(svclass == "DUP") %>% .$VF)
            if (max(this_chr[boundary_index]$VF) == max_vf)
              this_chr[boundary_index]$duplication_at_boundary_has_highest_VF <- "TRUE"

            {

              ## Classify episome. The boundary DUP has two breakends; take the
              ## PROXIMAL/DISTAL boundary by genomic position, not by row order.
              ## `breakpoints_gr` is built by stacking all chrom1 then all chrom2
              ## breakends, so `boundary_index[1]` is not necessarily the lower
              ## coordinate -- reading them in row order silently swapped prox/dist
              ## and made the "before" flank land INSIDE the amplicon (gained),
              ## failing the flank test on otherwise clean episomes.
              boundary_starts <- sort(this_chr[boundary_index]$start)
              prox_boundary <- boundary_starts[1]
              dist_boundary <- boundary_starts[length(boundary_starts)]

              ## Chromosomal segments < prox_boundary and > dist_boundary should
              ## not be gained/amplified ABOVE the chromosome's background level.
              ## `flank_baseline = "chromosome"` compares the flanks to the local
              ## per-chromosome baseline (the width-weighted median copy number of
              ## its non-focally-amplified segments), so a focal episome excised
              ## from a POLYSOMIC chromosome (e.g. EGFR on the gained chr7 that is
              ## near-universal in glioblastoma) is not missed just because its
              ## flanks sit above the tumour ploidy. `"ploidy"` keeps the flanks
              ## compared to the global sample ploidy.
              chr_cn <- this_sample_cnv_gr %>% gr2dt() %>%
                dplyr::filter(seqnames %in% unique_chrs[i])
              if (identical(flank_baseline, "chromosome")) {
                bg <- chr_cn %>% dplyr::filter(copyNumber < min_cn_ratio * this_sample_ploidy)
                if (nrow(bg) == 0) bg <- chr_cn
                w <- pmax(1, bg$end - bg$start); o <- order(bg$copyNumber)
                chr_baseline <- max(1, bg$copyNumber[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]])
                ## Flank gain threshold (per-chromosome baseline mode). A flank
                ## disqualifies the episome only if it is itself AMPLIFIED
                ## (>= min_cn_ratio * ploidy) -- i.e. the amplification continues
                ## past the boundary rather than dropping to background -- OR it
                ## exceeds the local chromosome-gain level (gain_ratio * baseline,
                ## which binds on a polysomic chromosome whose baseline is above the
                ## amplified bar). Taking the max of the two avoids two failure
                ## modes: (i) on a chromosome with broad LOSS the width-weighted
                ## baseline falls below diploid, so a genuinely near-diploid flank
                ## would spuriously read as "gained" (DO12742T1 PDGFRA, chr4 largely
                ## at CN ~1.3); (ii) a modest flank a little above diploid abutting a
                ## very-high-CN focal circle would be rejected even though the copy
                ## number has plainly dropped back to background (C3L-02504T1 PDGFRA,
                ## distal flank CN ~3.4 beside a CN-30 amplicon). The boundary DUP
                ## remains the primary episome signal.
                flank_thresh <- max(gain_ratio * chr_baseline,
                                    min_cn_ratio * this_sample_ploidy)
              } else {
                ## Strict global-ploidy baseline (as-published): compare flanks to
                ## the sample ploidy only, with no amplified-bar floor.
                chr_baseline <- this_sample_ploidy
                flank_thresh <- gain_ratio * chr_baseline
              }

              ## Flank test: is the chromosome gained on the OUTSIDE of the
              ## amplified region? Anchor to the copy-number amplicon edges
              ## (`min_amp_coord` / `max_amp_coord`, the first/last CN>3*ploidy
              ## segment), NOT the boundary-DUP breakend positions. The breakend
              ## routinely lands 1-2 bp inside the flanking segment (the segment
              ## ends at, or just past, the breakend), so keying the flank off the
              ## breakend coordinate dropped the true flank on an off-by-one and
              ## failed genuine episomes whose flanks are plainly diploid.
              ##
              ## Copy-number callers also emit a tiny (sub-kb) intermediate
              ## segment right at a sharp amplicon edge (a "segmentation
              ## shoulder"); using that abutting sliver as THE flank reads the
              ## transition, not the true flanking level. Skip segments narrower
              ## than `min_flank_width` and test the nearest SUBSTANTIVE flank
              ## instead, falling back to the abutting segment only if none
              ## qualifies.
              prox_flank <- chr_cn %>% dplyr::filter(end < min_amp_coord)
              prox_wide  <- prox_flank %>% dplyr::filter((end - start) >= min_flank_width)
              if (nrow(prox_wide) > 0) prox_flank <- prox_wide
              before_prox_boundary_not_gained <- nrow(prox_flank %>%
                                                         dplyr::filter(end == max(end)) %>%
                                                         dplyr::filter(copyNumber < flank_thresh)) > 0

              dist_flank <- chr_cn %>% dplyr::filter(start > max_amp_coord)
              dist_wide  <- dist_flank %>% dplyr::filter((end - start) >= min_flank_width)
              if (nrow(dist_wide) > 0) dist_flank <- dist_wide
              after_dist_boundary_not_gained <- nrow(dist_flank %>%
                                                        dplyr::filter(start == min(start)) %>%
                                                        dplyr::filter(copyNumber < flank_thresh)) > 0

              ## --- Diagnostic flags (independent of the flank test) ------------
              ## Founder boundary-DUP support: a genuine episome is initiated by
              ## the boundary DUP, which remains the highest-VF junction.
              bd_vf <- max(this_chr[boundary_index]$VF)

              ## Circularisation-junction microhomology. The boundary DUP is the
              ## self-ligation junction of the episome; its breakpoint homology
              ## length discriminates the DSB-repair pathway that sealed the circle
              ## (Eugen-Olsen, Hariprakash, Oestergaard & Regenberg, Nucleic Acids
              ## Res 2025; doi:10.1093/nar/gkaf122): near-blunt junctions (< 2 bp)
              ## are NHEJ, short microhomology (2 to < 14 bp) is alt-EJ / MMEJ, and
              ## long homology (>= 14 bp, the human HR minimum) is homologous
              ## recombination. Both breakends of a junction share one HOMLEN, so
              ## take the (non-missing) boundary-DUP homology.
              bd_hom <- suppressWarnings(
                max(as.numeric(this_chr[boundary_index]$homLen), na.rm = TRUE))
              if (is.finite(bd_hom)) {
                this_chr$boundary_homology <- bd_hom
                this_chr$junction_homology_class <-
                  if (bd_hom < mh_min_homology) "NHEJ"
                  else if (bd_hom < hr_min_homology) "MMEJ"
                  else "HR"
              }

              ## Internal SVs: events whose breakends ALL lie strictly inside the
              ## boundary-DUP span (excluding the boundary DUP itself). An internal
              ## SV whose VF exceeds the founder boundary DUP is flagged. A
              ## fold-back INVERSION doing so is a break-fusion-bridge signature
              ## and disqualifies the episome; an internal DELETION doing so is an
              ## ancestral (pre-circularisation) deletion -- kept, flag is
              ## informational only.
              int_ev <- this_chr[duplication_at_boundary != "TRUE",
                                 .(all_inside = all(start > prox_boundary & start < dist_boundary),
                                   cls = svclass[1], vf = max(VF)), by = event]
              int_ev <- int_ev[all_inside == TRUE]
              if (nrow(int_ev) && max(int_ev$vf) > bd_vf) {
                this_chr$flag_internal_sv_high_vf <- "TRUE"
                inv_ev <- int_ev[cls %in% c("h2hINV", "t2tINV")]
                if (nrow(inv_ev) && max(inv_ev$vf) > bd_vf)
                  this_chr$flag_internal_inversion <- "TRUE"
              }

              ## Bridging: an amplicon-connecting junction (a same-chromosome SV
              ## whose BOTH breakends sit at amplicon-level copy number) that spans
              ## a contiguous non-amplified gap wider than `bridge_gap` fuses two
              ## SEPARATE amplicons (more than one seed-gap apart) rather than
              ## bounding one contiguous episome. Genuine internal deletions leave
              ## only sub-`bridge_gap` gaps. Use the FULL sample breakpoints (not
              ## the `ext`-limited amplicon set): a bridge can reach a separate
              ## amplicon tens of Mb away, beyond `ext`, so its far breakend is
              ## absent from `this_chr`. Keep junctions with at least one breakend
              ## inside this amplicon's boundary span.
              chr_bp <- gr2dt(this_sample_breakpoints)
              chr_bp <- chr_bp[seqnames %in% unique_chrs[i]]
              chr_bp[, PURPLE_CN := as.numeric(PURPLE_CN)]
              jun <- chr_bp[, .(p1 = min(start), p2 = max(start), n = .N,
                                both_amp = all(PURPLE_CN > min_cn_ratio * this_sample_ploidy),
                                hits_amp = any(start >= prox_boundary & start <= dist_boundary)),
                            by = event]
              jun <- jun[n >= 2 & both_amp == TRUE & hits_amp == TRUE]
              if (nrow(jun)) {
                gapw <- vapply(seq_len(nrow(jun)), function(k) {
                  seg <- chr_cn %>% dplyr::filter(end >= jun$p1[k] & start <= jun$p2[k])
                  if (!nrow(seg)) return(0)
                  seg <- seg[order(seg$start), ]
                  w <- pmin(seg$end, jun$p2[k]) - pmax(seg$start, jun$p1[k])
                  amp <- seg$copyNumber > min_cn_ratio * this_sample_ploidy
                  runs <- rle(amp); e <- cumsum(runs$lengths); s0 <- e - runs$lengths + 1L
                  mx <- 0
                  for (m in seq_along(runs$values))
                    if (!runs$values[m]) mx <- max(mx, sum(w[s0[m]:e[m]]))
                  mx
                }, numeric(1))
                if (any(gapw > bridge_gap)) this_chr$flag_bridging_amplicon <- "TRUE"
              }

              ## Founder-boundary micronucleus test: among HIGH junction-copy-
              ## number (clonal, high-copy) intra-chromosomal SVs overlapping this
              ## amplicon, does any INVERSION begin INSIDE a DUP's span? Interleaved
              ## high-JCN DUP+INV -- an inversion whose proximal breakend is
              ## contained within a duplication's segment -- is the topological
              ## signature of interactions WITHIN a micronucleus (chromothripsis),
              ## not a clean circularisation, so the amplicon is not episomal. Only
              ## founder-level (high-JCN) junctions count; low-JCN subclonal
              ## interleaving (later internal rearrangement of a genuine episome) is
              ## ignored.
              hj <- chr_bp[JCN >= founder_jcn]
              if (nrow(hj)) {
                hj_ev <- hj[, .(lo = min(start), hi = max(start), cls = svclass[1], n = .N,
                                in_amp = any(start >= min_amp_coord & start <= max_amp_coord)),
                            by = event]
                hj_ev <- hj_ev[n >= 2 & in_amp == TRUE]
                dd <- hj_ev[cls == "DUP"]; ii <- hj_ev[cls %in% c("h2hINV", "t2tINV")]
                if (nrow(dd) && nrow(ii)) {
                  ## The interleaving DUP must itself span a non-amplified gap
                  ## (inter-locus / bridging): a local fold-back inversion inside
                  ## one contiguous amplicon is an internal rearrangement of a
                  ## genuine episome, NOT micronucleus. Only a DUP that reaches
                  ## across a gap between two loci, crossed by an inversion, is the
                  ## micronucleus signature.
                  dgap <- vapply(seq_len(nrow(dd)), function(a) {
                    seg <- chr_cn %>% dplyr::filter(end >= dd$lo[a] & start <= dd$hi[a])
                    if (!nrow(seg)) return(0)
                    seg <- seg[order(seg$start), ]
                    w <- pmin(seg$end, dd$hi[a]) - pmax(seg$start, dd$lo[a])
                    amp <- seg$copyNumber > min_cn_ratio * this_sample_ploidy
                    runs <- rle(amp); e <- cumsum(runs$lengths); s0 <- e - runs$lengths + 1L
                    mx <- 0
                    for (m in seq_along(runs$values))
                      if (!runs$values[m]) mx <- max(mx, sum(w[s0[m]:e[m]]))
                    mx
                  }, numeric(1))
                  micro <- FALSE
                  for (a in seq_len(nrow(dd))) if (dgap[a] > bridge_gap)
                    for (b in seq_len(nrow(ii)))
                      if (dd$lo[a] < ii$lo[b] && ii$lo[b] < dd$hi[a]) { micro <- TRUE; break }
                  if (micro) this_chr$flag_micronucleus <- "TRUE"
                }
              }

              ## NOTE: the original script tested
              ##   `after_dist_boundary_not_gained == after_dist_boundary_not_gained`
              ## which is a tautology (always TRUE). Corrected here to require
              ## BOTH flanks to be non-gained before calling an episome region.
              if (before_prox_boundary_not_gained & after_dist_boundary_not_gained) {
                this_chr$before_prox_boundary_not_gained <- "TRUE"
                this_chr$after_dist_boundary_not_gained <- "TRUE"
                this_chr[start >= prox_boundary & start <= dist_boundary,
                         episome_region := "TRUE"]
              }

              ## Define excision scar:
              prox_boundary_sv <- this_chr[boundary_index[1] - 1]$event
              dist_boundary_sv <- this_chr[boundary_index[2] + 1]$event

              if (length(prox_boundary_sv) > 0 & length(dist_boundary_sv) > 0) {
                if (!is.na(prox_boundary_sv) & !is.na(dist_boundary_sv)) {
                  if (prox_boundary_sv == dist_boundary_sv) {
                    if (unique(this_chr[event == prox_boundary_sv]$svclass == "DEL")) {
                      this_chr[event == prox_boundary_sv]$deletion_flanking_boundary <- "TRUE"
                    }
                  }
                }
              }
            }
          }
        }
      }

      ## Classify the whole amplicon (per chromosome):
      if (nrow(this_chr %>% dplyr::filter(episome_region == "TRUE")) > 0) {
        this_chr$episomal <- "TRUE"
      }
      ## The bridging and internal-high-VF-SV flags are informational only -- an
      ## internal DELETION that out-copies the boundary DUP is an ANCESTRAL event
      ## (a smaller-circle / pre-excision deletion), not a rival formation
      ## mechanism, so it does not override the episome call. Two signatures ARE
      ## disqualifying:
      ##  (i) a MICRONUCLEUS (interleaved high-JCN founder DUP+INV across a gap) --
      ##      intra-micronucleus chromothripsis, not a circularised episome; and
      ## (ii) an internal INVERSION that out-copies the boundary DUP -- the
      ##      inverted-duplication signature: the amplicon was founded by an
      ##      inverted duplication (a copy-gaining mechanism), not simple excision.
      if (any(this_chr$flag_micronucleus == "TRUE")) {
        this_chr$episomal <- "FALSE"
      }
      if (any(this_chr$flag_internal_inversion == "TRUE")) {
        this_chr$episomal <- "FALSE"
        this_chr$episomal_type2 <- "FALSE"
      }
      ## A chromosomal-bridge amplicon (edge terminates in a centromere) is
      ## excluded from the episomal call (and from the type-2 candidate set).
      if (any(this_chr$flag_chromosomal_bridge == "TRUE")) {
        this_chr$episomal <- "FALSE"
        this_chr$episomal_type2 <- "FALSE"
      }
      ## A boundary defined by an inversion or translocation (fold-back / complex
      ## junction) is not a circularisation -- excluded from episomal calling.
      if (any(this_chr$flag_inv_at_boundary == "TRUE") ||
          any(this_chr$flag_tra_at_boundary == "TRUE")) {
        this_chr$episomal <- "FALSE"
        this_chr$episomal_type2 <- "FALSE"
      }
      if (nrow(this_chr %>% dplyr::filter(deletion_flanking_boundary == "TRUE")) == 2) {
        this_chr$has_excision_scar <- "TRUE"
      }

      this_chr
    }) %>% data.table::rbindlist()
  }

  this_sample_breakpoints_ecdna_annotated
}


#' Call episomal extrachromosomal DNA from WGS structural variants
#'
#' Detects ecDNA amplicons whose structure is consistent with the *episome*
#' (breakage-independent) model of formation: a circular amplicon bounded by a
#' duplication breakpoint, arising from an otherwise non-amplified chromosomal
#' region, and often leaving a deletion "excision scar" at the origin locus.
#'
#' Each amplicon called by AmpliconArchitect (supplied in `ecdna_gr`, a
#' compulsory input) is processed independently (optionally in parallel). Within
#' each amplicon the function annotates the flanking structural variant
#' breakpoints with oncogene and allele-specific copy-number context, then
#' applies the heuristic described in [classify_amplicon_episomal()].
#'
#' @param ecdna_gr A [GenomicRanges::GRanges] of ecDNA amplicon regions with
#'   metadata columns `ID` (unique amplicon identifier) and `WGS_ID` (sample
#'   identifier) — typically the AmpliconArchitect amplicon catalogue. If `NULL`
#'   (the default), focal-amplicon seeds are detected from `cnv_gr` with
#'   [detect_amplicon_seeds()], so EpiTracer can run without AmpliconArchitect.
#' @param breakpoints_gr A [GenomicRanges::GRanges] of PURPLE (HMF pipeline) SV
#'   breakpoints with metadata columns `WGS_ID`, `event`, `svclass` (e.g. "DUP",
#'   "DEL"), `PURPLE_AF`, `PURPLE_JCN`, `VF`, `PURPLE_CN`, `insLen`, `HOMLEN`.
#' @param cnv_gr A [GenomicRanges::GRanges] of PURPLE (HMF pipeline)
#'   allele-specific copy-number segments with metadata columns `sample`,
#'   `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`.
#' @param cancer_genes_gr A [GenomicRanges::GRanges] of cancer gene loci with a
#'   `gene` (gene symbol) metadata column, used to annotate breakpoints with the
#'   overlapping oncogene.
#' @param ext Integer; number of base pairs to extend each amplicon by when
#'   searching for boundary breakpoints (default `1e7`).
#' @param mc.cores Integer; number of cores for [parallel::mclapply()]
#'   (default `1`). Values > 1 are ignored on Windows.
#' @param verbose Logical; print per-amplicon progress (default `FALSE`).
#' @param min_cn_ratio,seed_gap,seed_min_width Passed to
#'   [detect_amplicon_seeds()] when `ecdna_gr` is `NULL` (copy-number amplicon
#'   threshold `copyNumber > min_cn_ratio * ploidy`, gap to merge across, and
#'   minimum seed width). Ignored when `ecdna_gr` is supplied.
#' @param flank_baseline How the "flanks not gained" test is calibrated:
#'   `"chromosome"` (default) compares the amplicon flanks to the local
#'   per-chromosome baseline (the width-weighted median copy number of that
#'   chromosome's non-focally-amplified segments), `"ploidy"` compares them to
#'   the global sample ploidy. Use `"chromosome"` for focal episomes on a
#'   polysomic chromosome (e.g. EGFR on the gained chr7 in glioblastoma), which
#'   `"ploidy"` misses because the flanks sit above tumour ploidy.
#' @param gain_ratio Numeric; a flank is "gained" when its copy number is at or
#'   above `gain_ratio` times the baseline (default `1.4`).
#' @param min_flank_width Integer; ignore copy-number segments narrower than this
#'   (in bp) when reading the flanking copy number, so a sub-kb "segmentation
#'   shoulder" emitted at a sharp amplicon edge is not mistaken for the true
#'   flank (default `2000`). The nearest segment at least this wide is tested;
#'   if none qualifies, the abutting segment is used as a fallback.
#' @param bridge_gap Numeric; the maximum non-amplified gap (in bp) an
#'   amplicon-connecting junction may span before the amplicon is flagged as
#'   `bridging` and rejected (default: `seed_gap`, i.e. the same distance that
#'   separates distinct amplicon seeds). A junction whose both breakends sit at
#'   amplicon-level copy number but that spans a wider gap fuses two separate
#'   amplicons rather than bounding one episome; a genuine internal deletion
#'   leaves only a sub-`bridge_gap` gap.
#' @param founder_jcn Numeric; the junction-copy-number above which an SV counts
#'   as a clonal, high-copy "founder" junction for the micronucleus test (default
#'   `30`). When a founder-level inversion begins inside a founder-level
#'   duplication's span (interleaved DUP+INV), the amplicon is a micronucleus /
#'   chromothripsis product and is called non-episomal. Low-JCN (subclonal)
#'   interleaving -- later internal rearrangement of a genuine episome -- is
#'   ignored.
#' @param centromeres Optional [GenomicRanges::GRanges] of centromere spans (one
#'   per chromosome, e.g. from [load_centromeres()]). When supplied, an amplicon
#'   whose amplified edge terminates *within* a centromere (rather than
#'   traversing across it) is flagged `flag_chromosomal_bridge` and excluded from
#'   the episomal call as a dicentric / chromosomal-bridge event. `NULL` (default)
#'   skips the check.
#' @param mh_min_homology,hr_min_homology Numeric; breakpoint-homology (bp)
#'   thresholds used to classify the circularisation (boundary-DUP) junction into
#'   an inferred DSB-repair pathway, following Eugen-Olsen et al. (Nucleic Acids
#'   Res 2025; doi:10.1093/nar/gkaf122). A junction with homology
#'   `< mh_min_homology` is called `"NHEJ"` (near-blunt), `>= mh_min_homology` and
#'   `< hr_min_homology` is `"MMEJ"` (short microhomology / alternative
#'   end-joining), and `>= hr_min_homology` is `"HR"` (long homology / homologous
#'   recombination). Defaults `2` and `14` are the human microhomology and HR
#'   minimum-homology lengths reported in that review.
#'
#' @return A [data.table::data.table] combining the annotated breakpoints of all
#'   amplicons, with per-breakpoint classification columns including
#'   `duplication_at_boundary`, `duplication_at_boundary_has_highest_VF`,
#'   `episome_region`, `deletion_flanking_boundary`, `episomal`,
#'   `has_excision_scar`, and the diagnostic flags `flag_internal_sv_high_vf`
#'   (an internal SV out-VFs the boundary DUP; informational), the disqualifying
#'   `flag_internal_inversion` (a fold-back inversion does so -- break-fusion-
#'   bridge), and `flag_bridging_amplicon` (a junction fuses two separate
#'   amplicons) -- all character `"TRUE"`/`"FALSE"`. An amplicon carrying either
#'   disqualifying flag is set `episomal = "FALSE"`. For an amplicon with a
#'   boundary DUP it also reports the circularisation-junction microhomology as
#'   `boundary_homology` (numeric bp) and `junction_homology_class`
#'   (`"NHEJ"`/`"MMEJ"`/`"HR"`; `NA` when no boundary DUP is found), the inferred
#'   DSB-repair pathway that sealed the circle.
#'
#'   The other amplicon-formation mechanisms are computed by dedicated standalone
#'   callers, not here: [call_brf()] (breakage-replication/fusion / adjacent
#'   parallel breakpoints), [call_micronucleation()] (high-VF non-homologous
#'   translocation), [call_bfb()] (breakage-fusion-bridge) and
#'   [call_translocation_bridge_amp()] (translocation-bridge amplification). Join
#'   their output to this one by `WGS_ID` + `ID` to assemble a combined mechanism
#'   table.
#'
#' @examples
#' \dontrun{
#' ecdna_gr        <- readRDS("ecDNA_amplicon_regions.rds")
#' breakpoints_gr  <- readRDS("SV_catalogue.rds")
#' cnv_gr          <- readRDS("CN_segments.rds")
#' cancer_genes_gr <- readRDS("cancer_genes.rds")
#'
#' episomal <- call_simple_excision(
#'   ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
#'   ext = 1e7, mc.cores = 4
#' )
#' }
#'
#' @seealso [classify_amplicon_episomal()], [plot_sv_linear()]
#' @export
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
call_simple_excision <- function(ecdna_gr = NULL,
                                breakpoints_gr,
                                cnv_gr,
                                cancer_genes_gr,
                                ext = 1e7,
                                mc.cores = 1L,
                                verbose = FALSE,
                                min_cn_ratio = 3,
                                seed_gap = 1e6,
                                seed_min_width = 1e5,
                                flank_baseline = c("chromosome", "ploidy"),
                                gain_ratio = 1.4,
                                min_flank_width = 2000,
                                bridge_gap = seed_gap,
                                founder_jcn = 30,
                                centromeres = NULL,
                                mh_min_homology = 2,
                                hr_min_homology = 14) {

  flank_baseline <- match.arg(flank_baseline)
  stopifnot(
    methods::is(breakpoints_gr, "GRanges"),
    methods::is(cnv_gr, "GRanges"),
    methods::is(cancer_genes_gr, "GRanges")
  )

  ## Caller-agnostic column checks. EpiTracer works with SV/CN calls from any
  ## source (not only PURPLE / AmpliconArchitect) once coerced to these columns;
  ## fail early with an explicit message naming what is missing.
  .need <- function(gr, cols, what) {
    miss <- setdiff(cols, names(S4Vectors::mcols(gr)))
    if (length(miss))
      stop(sprintf("%s is missing required metadata column(s): %s",
                   what, paste(miss, collapse = ", ")), call. = FALSE)
  }
  .need(breakpoints_gr,
        c("WGS_ID", "event", "svclass", "PURPLE_AF", "PURPLE_JCN",
          "VF", "PURPLE_CN", "insLen", "HOMLEN"), "breakpoints_gr")
  .need(cnv_gr,
        c("sample", "copyNumber", "ploidy",
          "majorAlleleCopyNumber", "minorAlleleCopyNumber"), "cnv_gr")
  .need(cancer_genes_gr, "gene", "cancer_genes_gr")

  ## Amplicon seeds: use the supplied AmpliconArchitect catalogue, or (when
  ## ecdna_gr is NULL) detect focal amplicons from copy number alone.
  if (is.null(ecdna_gr)) {
    if (verbose) message("No ecdna_gr supplied; detecting focal-amplicon seeds from cnv_gr ...")
    ecdna_gr <- detect_amplicon_seeds(cnv_gr, min_cn_ratio = min_cn_ratio,
                                      gap = seed_gap, min_width = seed_min_width,
                                      breakpoints = breakpoints_gr)
    if (length(ecdna_gr) == 0L)
      stop("No focal amplicons detected in cnv_gr ",
           "(need copyNumber > ", min_cn_ratio, " x ploidy).", call. = FALSE)
  } else {
    stopifnot(methods::is(ecdna_gr, "GRanges"),
              !is.null(ecdna_gr$ID), !is.null(ecdna_gr$WGS_ID))
  }

  amplicon_ids <- unique(ecdna_gr$ID)
  if (verbose) message("Classifying ", length(amplicon_ids), " amplicons ...")

  results <- parallel::mclapply(seq_along(amplicon_ids), function(x) {
    if (verbose) message("[", x, "/", length(amplicon_ids), "] ", amplicon_ids[x])
    classify_amplicon_episomal(
      this_amplicon_id = amplicon_ids[x],
      ecdna_gr         = ecdna_gr,
      breakpoints_gr   = breakpoints_gr,
      cnv_gr           = cnv_gr,
      cancer_genes_gr  = cancer_genes_gr,
      ext              = ext,
      flank_baseline   = flank_baseline,
      gain_ratio       = gain_ratio,
      min_cn_ratio     = min_cn_ratio,
      min_flank_width  = min_flank_width,
      bridge_gap       = bridge_gap,
      founder_jcn      = founder_jcn,
      centromeres      = centromeres,
      mh_min_homology  = mh_min_homology,
      hr_min_homology  = hr_min_homology,
      verbose          = verbose
    )
  }, mc.cores = mc.cores)

  data.table::rbindlist(results, fill = TRUE)
}
