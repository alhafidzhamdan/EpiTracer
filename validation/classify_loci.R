## ---------------------------------------------------------------------------
## Per-locus ecDNA mechanism classification (refined EpiTracer heuristic)
##
## Reads one AmpliconArchitect amplicon graph and classifies EACH major
## amplified locus by mechanism, rather than giving one call for the whole
## (possibly chimeric) amplicon. Calibrated on the CCLE episomal calls against
## expert ground truth (5637, KATO III, NCI-H2170, NCI-H526).
##
## A locus is `episomal` iff:
##   1. a boundary DUP spans it and is its highest-VF DUP (self-ligated circle),
##   2. its flanks are diploid (relative to the per-chromosome baseline),
##   3. it is NOT breakage-fusion-bridge: < `foldback_min` near-self inverted
##      (fold-back) junctions, and
##   4. no inter-chromosomal (TRA) junction reaches the boundary DUP's copy
##      number (a dominant TRA = translocation-formed, not a circle).
## Otherwise it is tagged by mechanism:
##   BFB                  - fold-back inversions, no amplified translocation.
##   translocation-bridge - an amplified inter-chromosomal translocation makes a
##                          dicentric bridge that amplifies, few fold-backs
##                          (Lee et al., Nature 2023).
##   LTA                  - an amplified translocation PLUS fold-backs (BFB cycles);
##                          `arm_loss` flags the sub-baseline loss on the amplicon
##                          chromosome (Espejo Valle-Inclan et al., Cell 2024).
##   complex              - none of the above (no clean boundary DUP / too many junctions).
## Inter-locus TRAs BELOW the per-locus boundary DUPs are episome fusion: the
## episomal loci of such an amplicon are flagged `fused = TRUE` (e.g. 5637).
##
## Depends on read_aa_graph()/aa_svclass() from aa_to_epitracer.R.
## ---------------------------------------------------------------------------

## classify every major locus of one amplicon graph; returns a data.frame
classify_loci <- function(graph_path, sample_id, amplicon = NA_character_,
                          cnv_bed = NULL,
                          min_locus_width = 2e5, min_cn = 6,
                          boundary_tol_frac = 0.05, foldback_dist = 5e4, foldback_min = 3L,
                          gain_mult = 1.4, min_break_reads = 5,
                          min_bdup_cn = 2, complexity_max = 3L, hi_junc_frac = 0.25,
                          centromeres = NULL, centromere_pad = 5e5,
                          tra_frac = 0.25, loss_frac = 0.6, loss_min_width = 1e6) {
  g   <- read_aa_graph(graph_path)
  seg <- g$segments
  br  <- g$breaks[g$breaks$type == "discordant" & g$breaks$n_reads >= min_break_reads, , drop = FALSE]
  if (nrow(br)) br$svclass <- vapply(seq_len(nrow(br)), function(i)
    aa_svclass(br$chr1[i], br$pos1[i], br$strand1[i], br$chr2[i], br$pos2[i], br$strand2[i]),
    character(1)) else br$svclass <- character(0)

  ## genome-wide CN (for flanks + per-chr baseline); else fall back to AA segments
  if (!is.null(cnv_bed) && file.exists(cnv_bed)) {
    bed <- utils::read.delim(cnv_bed, header = FALSE, stringsAsFactors = FALSE)
    gw <- data.frame(chr = bed[[1]], start = as.numeric(bed[[2]]), end = as.numeric(bed[[3]]),
                     cn = as.numeric(bed[[ncol(bed)]]), stringsAsFactors = FALSE)
  } else {
    gw <- data.frame(chr = seg$chr, start = seg$start, end = seg$end, cn = seg$cn)
  }
  baseline <- function(ch) {
    v <- gw[gw$chr == ch & gw$cn < min_cn, ]
    if (!nrow(v)) return(2)
    max(2, round(stats::median(rep(v$cn, pmax(1, round((v$end - v$start) / 1e4))))))
  }

  ## major loci: per chromosome, merge CN>min_cn AA segments; keep width > thresh
  loci <- list()
  for (ch in unique(seg$chr)) {
    s <- seg[seg$chr == ch & seg$cn > min_cn, , drop = FALSE]
    if (!nrow(s) || sum(s$end - s$start) < min_locus_width) next
    loci[[ch]] <- c(lo = min(s$start), up = max(s$end),
                    maxcn = max(s$cn), width = sum(s$end - s$start))
  }
  if (!length(loci)) return(NULL)

  out <- lapply(names(loci), function(ch) {
    L <- loci[[ch]]; lo <- L["lo"]; up <- L["up"]; span <- up - lo; base <- baseline(ch)
    intra <- br[br$chr1 == ch & br$chr2 == ch &
                pmin(br$pos1, br$pos2) <= up & pmax(br$pos1, br$pos2) >= lo, , drop = FALSE]
    ## boundary DUP: intra-chr DUP spanning most of the locus, near both ends
    ## boundary DUP: a DUP whose BOTH ends sit at the amplified-region edges (a
    ## true self-ligated circle), not one that stops short of the boundary or
    ## overshoots far past it into flanking sequence.
    bd_cn <- NA_real_; bd_vf <- NA_real_
    btol <- max(2e5, boundary_tol_frac * span)
    d <- intra[intra$svclass == "DUP", , drop = FALSE]
    if (nrow(d)) {
      a <- pmin(d$pos1, d$pos2); b <- pmax(d$pos1, d$pos2)
      ok <- abs(a - lo) <= btol & abs(b - up) <= btol & d$edge_cn >= min_bdup_cn
      if (any(ok)) { j <- which(ok)[which.max(d$edge_cn[ok])]; bd_cn <- d$edge_cn[j]; bd_vf <- d$n_reads[j] }
    }
    ## fold-back (BFB) content: near-self inverted junctions
    inv <- intra[intra$svclass %in% c("h2hINV", "t2tINV"), , drop = FALSE]
    fb_near <- if (nrow(inv)) sum(abs(inv$pos1 - inv$pos2) < foldback_dist) else 0L
    ## junction complexity: competing high-copy junctions at the locus. A simple
    ## episome has one dominant circularisation (+ maybe a scar); many high-CN
    ## junctions mean internal rearrangement (complex/BFB), not a clean circle.
    n_hi_junc <- sum(intra$edge_cn >= hi_junc_frac * L["maxcn"])
    ## flanks: nearest genome-wide segment on each side outside the locus
    left  <- gw[gw$chr == ch & gw$end   <= lo, ]; right <- gw[gw$chr == ch & gw$start >= up, ]
    lflank <- if (nrow(left))  left$cn[which.max(left$end)]    else base
    rflank <- if (nrow(right)) right$cn[which.min(right$start)] else base
    diploid_flanks <- lflank < gain_mult * base & rflank < gain_mult * base
    ## dominant inter-chromosomal junction touching this chr
    tra <- br[br$svclass == "TRA" & (br$chr1 == ch | br$chr2 == ch), , drop = FALSE]
    max_tra_cn <- if (nrow(tra)) max(tra$edge_cn) else 0
    substantial_tra <- max_tra_cn >= tra_frac * L["maxcn"]   # an amplified boundary translocation
    bfb <- fb_near >= foldback_min                            # fold-back inversions = BFB cycles
    ## centromere adjacency: an amplicon abutting the centromere is more likely
    ## pericentromeric/BFB than an excised episome.
    cen_adj <- FALSE
    if (!is.null(centromeres)) {
      cc <- centromeres[centromeres$chr == ch, , drop = FALSE]
      if (nrow(cc)) cen_adj <- any(up >= cc$start - centromere_pad & lo <= cc$end + centromere_pad)
    }
    ## a sizeable sub-baseline (lost) segment on the amplicon's chromosome — the
    ## "loss" limb of loss-translocation-amplification.
    armseg <- gw[gw$chr == ch, , drop = FALSE]
    has_arm_loss <- any(armseg$cn < loss_frac * base & (armseg$end - armseg$start) > loss_min_width)

    ## Mechanism. Both LTA and translocation-bridge (TB) are driven by an amplified
    ## inter-chromosomal translocation; they are distinguished by fold-back content
    ## (LTA proceeds through BFB cycles -> fold-backs; TB's translocation makes the
    ## dicentric directly -> few fold-backs). See Espejo Valle-Inclan et al. Cell
    ## 2024 (LTA) and Lee et al. Nature 2023 (TB).
    mech <- if (bfb && substantial_tra) "LTA"
    else if (is.na(bd_cn)) {
      if (substantial_tra) "translocation-bridge"
      else if (bfb) "BFB" else "complex"
    } else if (bfb) "BFB"
    else if (n_hi_junc > complexity_max) "complex"       # thicket of junctions, not a clean circle
    else if (!diploid_flanks) "complex"
    else if (cen_adj) "complex"                          # pericentromeric, not a clean episome
    else if (max_tra_cn > bd_cn) "translocation-bridge"
    else "episomal"

    data.frame(sample = sample_id, amplicon = amplicon, chr = ch,
               start = unname(lo), end = unname(up), max_cn = unname(L["maxcn"]),
               baseline = base, boundary_dup_cn = bd_cn, boundary_dup_vf = bd_vf,
               foldbacks = fb_near, n_hi_junc = n_hi_junc, diploid_flanks = diploid_flanks,
               centromere_adjacent = cen_adj, max_tra_cn = max_tra_cn, arm_loss = has_arm_loss,
               mechanism = mech, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  ## fusion: an amplicon whose episomal loci are joined by (subclonal) TRAs
  n_epi <- sum(res$mechanism == "episomal")
  has_tra <- any(br$svclass == "TRA")
  res$fused <- res$mechanism == "episomal" & (n_epi > 1 | has_tra)
  res
}
