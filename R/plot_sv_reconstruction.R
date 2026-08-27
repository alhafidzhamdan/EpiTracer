#' Temporally-stratified rearrangement plot (SV+CN grouped by read support)
#'
#' A companion to [plot_sv_linear()] that decomposes a locus's structural
#' rearrangements into read-support (variant-fraction, `VF`) strata and draws one
#' stacked panel per stratum. Junctions with the highest read support -- taken as
#' a proxy for the earliest / highest-copy events -- are plotted at the top, then
#' successively lower-`VF` strata below, so the figure reads as a rough temporal
#' decomposition of amplicon evolution.
#'
#' Strata are found by one-dimensional clustering of `log(VF)` (k-means; the
#' optimal 1-D solver \pkg{Ckmeans.1d.dp} is used when installed, otherwise
#' [stats::kmeans]). `VF` here is the junction read-support count (as emitted by
#' PURPLE), which is strongly right-skewed, so clustering is done on the log scale
#' -- equal-width `VF` bins would collapse almost every junction into the lowest
#' bracket. All panels share the same concatenated-locus x-axis and the same
#' copy-number / read-support y-scaling, so arc heights are comparable across
#' strata; each panel draws only the copy-number segments lying near that
#' stratum's breakpoints, emphasising the CN change each set of events bounds.
#'
#' Loci are resolved exactly as in [plot_sv_linear()] (explicit `loci`;
#' `chromosome` + `chromosome_range`; `chromosome` alone; or auto-detection of all
#' amplified loci).
#'
#' @inheritParams plot_sv_linear
#' @param vf_col Name of the read-support column in `sv_data` (default `"VF"`).
#' @param k Number of `VF` strata. `"auto"` (default) chooses the number of
#'   clusters automatically (via \pkg{Ckmeans.1d.dp} BIC when available, otherwise
#'   a within-cluster-sum-of-squares elbow capped at `max_k`); an integer forces
#'   that many. Reduced automatically when a sample has fewer distinct `VF` values
#'   than `k`.
#' @param max_k Upper bound on the number of strata when `k = "auto"` (default
#'   `4`).
#' @param vf_breaks Optional numeric vector of explicit `VF` cut points (upper
#'   edges), e.g. `c(200, 80)` for strata `>200`, `80-200`, `<80`. Overrides `k`.
#' @param min_vf Junctions with `VF < min_vf` are dropped before clustering
#'   (default `1`).
#' @param isolate_founder Logical (default `TRUE`); pull the single highest-`VF`
#'   junction (the defining excision) out of the clustering and give it its own
#'   top panel, headed "Max VF". Its arc is drawn at the same weight as every
#'   other junction. It is excluded from the k-means / breaks step (so one
#'   dominant outlier no longer skews the clusters), and the remaining strata are
#'   headed "Cluster 1..X" (highest `VF` first). Ties at the maximum `VF` are
#'   treated as one group. See `founder_offscale` for how it relates to the axis.
#' @param founder_offscale Logical (default `FALSE`); controls the shared
#'   read-support axis when a founder is isolated. `FALSE`: the axis reflects the
#'   true maximum `VF` (founder included), so the founder arc's height is literal.
#'   `TRUE`: the axis is driven by the highest non-founder junction (so the lower
#'   strata are not compressed by a large founder outlier) and the founder arc is
#'   clamped to the axis top -- its true `VF` is still labelled, and every panel's
#'   read-support axis title is marked "(scaled)" to flag that the axis no longer
#'   maps 1:1 to `VF`. Preferred when the founder `VF` greatly exceeds the next
#'   highest.
#' @param vf_scale How arc heights (and the right-hand read-support axis) are
#'   scaled. `"shared"` (default) puts every panel on one global
#'   read-support-to-copy-number scale, so arc heights are directly comparable
#'   across strata -- best for multi-locus / hub amplicons. `"per_panel"` scales
#'   each panel to its own stratum's `VF` range, so a single very-high-`VF`
#'   junction does not crush the dynamic range of the lower strata -- usually the
#'   better choice for a dense single-chromosome amplicon. The left copy-number
#'   axis and the x-axis stay shared either way, so panels remain aligned.
#' @param cn_near_flank Numeric; in `cn_display = "actual"`, a copy-number
#'   segment is drawn in a panel if it lies within this many bp of any of that
#'   stratum's breakpoints (and, for intra-locus junctions, if it lies between
#'   the two breakpoints). Default `1e5`.
#' @param cn_display How the major-allele copy number is drawn per panel.
#'   `"reconstruct"` (default) shows an *iterative* reconstruction of how the
#'   copy-number profile is built up wave by wave. The breakpoints of the
#'   junctions introduced so far (strata \eqn{1..k}) partition each locus into
#'   intervals, and every interval is drawn flat at the minimum of the two
#'   copy-number segments immediately adjacent to (just inside) its bounding
#'   breakpoints. The founder panel therefore shows a single flat baseline over
#'   the founder span; each lower panel introduces more breakpoints, so the
#'   intervals subdivide and copy-number structure emerges; and the bottom panel
#'   is the full observed bulk copy number that the reconstruction evolves
#'   toward. This is a visual model, not a formal amplicon deconvolution.
#'   `"actual"` instead draws the real observed CN segments near each stratum's
#'   breakpoints in every panel (see `cn_near_flank`). The minor-allele CN is
#'   always drawn in full regardless of this setting. In `"reconstruct"` the
#'   allele-loss bars (LOH where minor CN `< loh_thresh`, and homozygous
#'   deletions where total CN `< homdel_thresh`) are placed on the same timeline:
#'   a loss segment is revealed in the wave of the earliest intra-locus deletion
#'   junction whose deleted span covers it, or -- if no deletion covers it -- from
#'   the founder wave down (treated as pre-existing).
#' @param prior_sv_alpha Numeric in `[0, 1]`; opacity at which the junctions of
#'   earlier (higher-`VF`) strata are re-drawn faintly in each lower panel, so
#'   the accumulation of rearrangements is visible going down. The current
#'   stratum is always drawn at full opacity. Set to `0` to show only the current
#'   stratum's junctions in each panel (default `0.15`).
#' @param founder_alpha Numeric in `[0, 1]`; opacity at which the founder
#'   junction (the defining, highest-`VF` event) is re-drawn as a prior in the
#'   lower panels. Kept higher than `prior_sv_alpha` so the founder stays
#'   trackable all the way down the reconstruction, while later waves fade more.
#'   The founder arc is also drawn with a slightly bolder line than the other
#'   junctions. Has no effect on the founder's own (top) panel, where it is
#'   always at full opacity (default `0.5`).
#' @param cn_border_lines Logical (default `TRUE`); draw vertical dashed guide
#'   lines at the borders of the final (fully reconstructed) copy-number
#'   segments. Because the panels share the x-axis, the guides line up across the
#'   whole stack, marking where copy-number transitions occur.
#' @param cn_border_min_step Numeric minimum copy-number change for an *internal*
#'   border line to be drawn (each locus's outer edges are always drawn). If
#'   `NULL` (default) a proportional threshold (`max(0.25, 0.015 * cn_axis_max)`)
#'   is used, so a guide is drawn at essentially every perceptible copy-number
#'   step (only flat, equal-level segment splits are skipped). Set larger to show
#'   only the biggest transitions, smaller to also mark sub-copy noise.
#' @param drop_empty_strata Logical; drop strata that end up with no drawable
#'   in-locus junction (default `TRUE`).
#' @param panel_rel_height Numeric relative height of the (label-bearing) bottom
#'   panel versus the others (default `1.4`).
#'
#' @return A \pkg{patchwork} object stacking the per-stratum panels (invisibly
#'   when a file is written); the written path is attached as attribute `"path"`,
#'   and the per-junction stratum assignment as attribute `"strata"`.
#'
#' @examples
#' \dontrun{
#' plot_sv_reconstruction("DUMC12T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
#'               k = "auto", outdir = "plots")
#'
#' # Fixed read-support brackets:
#' plot_sv_reconstruction("DUMC12T1", cnv, sv, karyotype = K, gene_coord = G,
#'               vf_breaks = c(200, 80))
#' }
#' @seealso [plot_sv_linear()], [call_simple_excision()]
#' @export
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom data.table rbindlist
#' @importFrom grDevices pdf dev.off
#' @importFrom stats kmeans median
#' @importFrom utils read.table head
plot_sv_reconstruction <- function(sample,
                          cnv_data,
                          sv_data,
                          wgd_data = NULL,
                          genome = c("hg38", "hg19", "mm10"),
                          karyotype = NULL,
                          gene_coord = NULL,
                          chromosome = NULL,
                          chromosome_range = NULL,
                          loci = NULL,
                          events = "amp",
                          cluster_gap = 5e6,
                          flank_pct = 10,
                          min_cn_ratio = 3,
                          gain_ratio = 1.4,
                          loh_thresh = 0.5,
                          homdel_thresh = 0.5,
                          min_amp_width = 1e5,
                          gap_frac = 0.06,
                          vf_col = "VF",
                          k = "auto",
                          max_k = 4,
                          vf_breaks = NULL,
                          min_vf = 1,
                          isolate_founder = TRUE,
                          founder_offscale = FALSE,
                          vf_scale = c("shared", "per_panel"),
                          cn_near_flank = 1e5,
                          cn_display = c("reconstruct", "actual"),
                          prior_sv_alpha = 0.15,
                          founder_alpha = 0.5,
                          cn_border_lines = TRUE,
                          cn_border_min_step = NULL,
                          drop_empty_strata = TRUE,
                          panel_rel_height = 1.4,
                          genes_to_highlight = NULL,
                          gene_label_angle = NULL,
                          repel_labels = TRUE,
                          displayExon = FALSE,
                          cds_gr = NULL,
                          cn_max = NULL,
                          offset_gene = 1.15,
                          ymax_highlight_ratio = 1.08,
                          karyotype_rel_size = 0.048,
                          loh_position_ratio = 0.5,
                          highlight_amp = TRUE,
                          highlight_hom_del = TRUE,
                          highlight_events = NULL,
                          highlight_id_col = NULL,
                          highlight_colour = "#d95f0e",
                          dim_unhighlighted = FALSE,
                          wgd_sample_col = NULL,
                          outdir = NULL,
                          save = TRUE,
                          plot_width_custom = NULL,
                          plot_height_custom = NULL,
                          verbose = FALSE) {

  vf_scale <- match.arg(vf_scale)
  cn_display <- match.arg(cn_display)

  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("plot_sv_reconstruction() needs the 'patchwork' package to stack panels. ",
         "Install it with install.packages('patchwork').")

  this_sample <- sample

  ## ---- Reference data (identical handling to plot_sv_linear) --------------
  ## `genome` selects the bundled karyotype + oncogene panel (hg38 / hg19 /
  ## mm10); an explicit `karyotype` / `gene_coord` overrides.
  genome <- match.arg(genome)
  if (is.null(karyotype))
    karyotype <- system.file("extdata", paste0("chr_info_", genome, ".rds"), package = "EpiTracer")
  if (is.null(gene_coord))
    gene_coord <- system.file("extdata", paste0("oncogene_coord_", genome, ".bed"), package = "EpiTracer")
  if (is.character(karyotype)) {
    if (!nzchar(karyotype) || !file.exists(karyotype))
      stop("Karyotype reference not found. Pass `karyotype` (a data.frame or .rds path).")
    karyotype <- readRDS(karyotype)
  }
  karyotype <- as.data.frame(karyotype)
  names(karyotype)[1] <- "chr"
  karyotype$color[karyotype$gieStain == "gneg"]    <- "white"
  karyotype$color[karyotype$gieStain == "gpos25"]  <- "grey75"
  karyotype$color[karyotype$gieStain == "gpos50"]  <- "grey50"
  karyotype$color[karyotype$gieStain == "gpos75"]  <- "grey25"
  karyotype$color[karyotype$gieStain == "gpos100"] <- "grey0"
  karyotype$color[karyotype$gieStain == "acen"]    <- "red"
  chr_len <- tapply(karyotype$end, karyotype$chr, max)

  if (is.character(gene_coord)) {
    if (!nzchar(gene_coord) || !file.exists(gene_coord))
      stop("Gene-coordinate reference not found. Pass `gene_coord` (a data.frame or BED path).")
    gene_coord <- utils::read.table(gene_coord, header = FALSE, sep = "\t")
  }
  gene_coord <- as.data.frame(gene_coord)
  names(gene_coord) <- c("chr", "start", "end", "strand", "gene")
  ## default to the whole loaded panel, so labels adapt to `genome`
  genes <- if (is.null(genes_to_highlight)) unique(gene_coord$gene) else genes_to_highlight

  ## ---- Sample data --------------------------------------------------------
  cnv <- cnv_data[cnv_data$sample == this_sample, , drop = FALSE]
  cnv$chr <- as.character(cnv$seqnames)
  if (nrow(cnv) == 0) stop("No copy-number data for sample ", this_sample, ".")

  wgd_status <- NULL
  if (!is.null(wgd_data)) {
    wgd_data <- as.data.frame(wgd_data)
    wgd_col <- if (!is.null(wgd_sample_col)) wgd_sample_col else
      if ("sample" %in% names(wgd_data)) "sample" else "WGS_ID"
    wrow <- wgd_data[wgd_data[[wgd_col]] == this_sample, , drop = FALSE]
    wgd_status <- if (nrow(wrow) >= 1 && wrow$Polyploidy[1] == "No") "Diploid" else "WGD"
  }

  ## ---- Resolve loci (identical to plot_sv_linear) -------------------------
  events <- match.arg(events, c("amp", "gain", "loh", "homdel"), several.ok = TRUE)
  is_event <- function(d, ploidy) {
    keep <- rep(FALSE, nrow(d))
    if ("amp"    %in% events) keep <- keep | (d$copyNumber > min_cn_ratio * ploidy)
    if ("gain"   %in% events) keep <- keep | (d$copyNumber > gain_ratio * ploidy &
                                              d$copyNumber <= min_cn_ratio * ploidy)
    if ("loh"    %in% events) keep <- keep | (d$minorAlleleCopyNumber < loh_thresh)
    if ("homdel" %in% events) keep <- keep | (d$copyNumber < homdel_thresh)
    keep
  }
  detect_loci <- function(chroms = NULL) {
    ploidy <- cnv$ploidy[1]
    ev <- cnv[is_event(cnv, ploidy), , drop = FALSE]
    clusters <- list()
    for (cc in unique(ev$chr)) {
      d <- ev[ev$chr == cc, , drop = FALSE]
      d <- d[order(d$start), , drop = FALSE]
      prev_max_end <- cummax(d$end)
      newgrp <- c(TRUE, (d$start[-1] - prev_max_end[-nrow(d)]) > cluster_gap)
      grp <- cumsum(newgrp)
      for (g in unique(grp)) {
        dd <- d[grp == g, , drop = FALSE]
        clusters[[length(clusters) + 1L]] <- data.frame(
          chr = cc, start = min(dd$start), end = max(dd$end),
          amp_bp = sum(dd$end - dd$start), stringsAsFactors = FALSE)
      }
    }
    agg <- if (length(clusters)) do.call(rbind, clusters) else
      data.frame(chr = character(), start = numeric(), end = numeric(), amp_bp = numeric())
    agg <- agg[agg$amp_bp >= min_amp_width, , drop = FALSE]
    if (nrow(agg) > 0) {
      w <- agg$end - agg$start
      flank <- (flank_pct / 100) * w
      agg$start <- pmax(0, round(agg$start - flank))
      agg$end   <- pmin(chr_len[agg$chr], round(agg$end + flank))
    }
    if (is.null(chroms)) {
      if (nrow(agg) == 0) stop("No ", paste(events, collapse = "/"), " loci found for ",
                               this_sample, " (>= ", min_amp_width,
                               " bp). Specify `chromosome` or `loci`.")
      res <- agg[order(suppressWarnings(as.integer(gsub("chr", "", agg$chr))), agg$chr, agg$start),
                 c("chr", "start", "end")]
    } else {
      res <- do.call(rbind, lapply(chroms, function(cc) {
        if (cc %in% agg$chr) agg[agg$chr == cc, c("chr", "start", "end"), drop = FALSE]
        else data.frame(chr = cc, start = 0, end = unname(chr_len[cc]), stringsAsFactors = FALSE)
      }))
    }
    res
  }

  if (!is.null(loci)) {
    if (is.character(loci)) {
      m <- regmatches(loci, regexec("^(chr[^:]+):([0-9.eE+]+)-([0-9.eE+]+)$", loci))
      loci <- do.call(rbind, lapply(m, function(x)
        data.frame(chr = x[2], start = as.numeric(x[3]), end = as.numeric(x[4]),
                   stringsAsFactors = FALSE)))
    }
    loci <- as.data.frame(loci)[, c("chr", "start", "end")]
  } else if (!is.null(chromosome_range)) {
    cr <- chromosome_range
    if (is.null(dim(cr))) cr <- matrix(cr, ncol = 2, byrow = TRUE)
    loci <- data.frame(chr = chromosome, start = cr[, 1], end = cr[, 2],
                       stringsAsFactors = FALSE)
  } else if (!is.null(chromosome)) {
    loci <- detect_loci(chromosome)
  } else {
    loci <- detect_loci(NULL)
  }
  loci$width <- loci$end - loci$start
  if (verbose) {
    message("Loci:")
    for (i in seq_len(nrow(loci)))
      message(sprintf("  %s:%.0f-%.0f (%.1f Mb)", loci$chr[i], loci$start[i],
                      loci$end[i], loci$width[i] / 1e6))
  }

  ## ---- Concatenated layout (shared across every panel) --------------------
  total_w <- sum(loci$width)
  gap <- total_w * gap_frac
  loci$offset <- utils::head(c(0, cumsum(loci$width + gap)), nrow(loci))
  loci$gx_start <- loci$offset
  loci$gx_end   <- loci$offset + loci$width
  loci$gx_mid   <- (loci$gx_start + loci$gx_end) / 2

  locus_idx <- function(chr, pos) {
    idx <- rep(NA_integer_, length(pos))
    for (i in seq_len(nrow(loci))) {
      hit <- chr == loci$chr[i] & pos >= loci$start[i] & pos <= loci$end[i]
      idx[hit] <- i
    }
    idx
  }
  ## genomic position -> concatenated x
  gx_of <- function(li, pos) loci$offset[li] + (pos - loci$start[li])

  ## ---- CN segments in plot coords (full backbone; panels subset it) -------
  cn_plot <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    d <- cnv[cnv$chr == loci$chr[i] & cnv$end >= loci$start[i] & cnv$start <= loci$end[i], , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$start <- pmax(d$start, loci$start[i]); d$end <- pmin(d$end, loci$end[i])
    d$gx_start <- loci$offset[i] + (d$start - loci$start[i])
    d$gx_end   <- loci$offset[i] + (d$end   - loci$start[i])
    d$li <- i
    d
  }))
  if (is.null(cn_plot) || nrow(cn_plot) == 0) stop("No copy-number data in the selected loci.")

  ploidy <- cnv$ploidy[1]
  max.cn <- max(cn_plot$majorAlleleCopyNumber)

  ## `nice_ceiling()` lives in R/utils-axis.R (shared with plot_sv_linear).
  ## Floor the axis top at 2 (matching plot_sv_linear) so a near-diploid region
  ## uses a full 0-2 axis; harmless for amplicons where max.cn is far above 2.
  cn_axis_max <- if (!is.null(cn_max)) cn_max else max(2, nice_ceiling(max.cn))

  ## ---- SVs ----------------------------------------------------------------
  sv <- sv_data[sv_data$sample == this_sample, , drop = FALSE]
  if (nrow(sv) == 0) stop("No structural-variant data for sample ", this_sample, ".")
  if (!vf_col %in% names(sv))
    stop("Read-support column '", vf_col, "' not found in `sv_data`. ",
         "Set `vf_col` to one of: ", paste(names(sv), collapse = ", "), ".")
  sv$VF  <- suppressWarnings(as.numeric(sv[[vf_col]]))
  .chrp <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))  # accept 'chr'-prefixed or bare
  sv$chr1 <- .chrp(sv$chrom1); sv$chr2 <- .chrp(sv$chrom2)
  sv$pos1 <- as.integer(sv$start1);    sv$pos2 <- as.integer(sv$start2)
  sv$strands <- paste0(sv$strand1, sv$strand2)
  sv$strands <- ifelse(sv$svclass == "TRA", "TRA", sv$strands)
  sv$l1 <- locus_idx(sv$chr1, sv$pos1)
  sv$l2 <- locus_idx(sv$chr2, sv$pos2)
  ## keep junctions with at least one end in a locus and a usable VF
  sv <- sv[!(is.na(sv$l1) & is.na(sv$l2)) & is.finite(sv$VF) & sv$VF >= min_vf, , drop = FALSE]
  if (nrow(sv) == 0)
    stop("No structural variants with VF >= ", min_vf, " fall inside the selected loci for ",
         this_sample, ".")

  DEL_colour <- "#bde0fe"; DUP_colour <- "#c1447e"; h2hINV_colour <- "#a5a6ae"
  t2tINV_colour <- "#384351"; TRA_colour <- "#fac881"
  sv_col <- function(s) ifelse(s %in% c("+-", "DEL"), DEL_colour,
                        ifelse(s %in% c("++", "h2hINV"), h2hINV_colour,
                        ifelse(s %in% c("--", "t2tINV"), t2tINV_colour,
                        ifelse(s %in% c("-+", "DUP"), DUP_colour, TRA_colour))))
  sv$colour <- sv_col(sv$strands)

  ## Optional event highlight (shared with plot_sv_linear / plot_sv_circos):
  ## matched SVs are recoloured to `highlight_colour` and drawn bolder (their
  ## per-row `hl` feeds the line-width multiplier below); `dim_unhighlighted`
  ## greys the rest. Carried on `sv` so it survives subsetting into VF strata.
  sv$hl <- .resolve_highlight(sv, highlight_events, highlight_id_col)
  if (any(sv$hl)) {
    if (isTRUE(dim_unhighlighted)) sv$colour[!sv$hl] <- "grey85"
    sv$colour[sv$hl] <- highlight_colour
  }
  highlight_lwd_mult <- 3

  ## ---- Cluster junctions into VF (read-support) strata --------------------
  ## Right-skewed counts -> cluster on log scale. Returns integer stratum id per
  ## row of `sv`, then order strata by descending median VF (earliest on top).
  assign_strata <- function(vf) {
    lv <- log(vf)
    uniq_n <- length(unique(round(lv, 8)))
    if (!is.null(vf_breaks)) {
      edges <- sort(unique(c(-Inf, vf_breaks, Inf)))
      return(as.integer(cut(vf, breaks = edges, labels = FALSE, include.lowest = TRUE)))
    }
    if (uniq_n == 1L) return(rep(1L, length(vf)))
    if (identical(k, "auto")) {
      kk <- min(max_k, uniq_n)
      if (requireNamespace("Ckmeans.1d.dp", quietly = TRUE)) {
        cl <- Ckmeans.1d.dp::Ckmeans.1d.dp(lv, k = c(1L, kk))$cluster
        return(as.integer(cl))
      }
      ## too few points for k-means (centres would meet/exceed the point count,
      ## which stats::kmeans rejects): rank into bins instead
      if (kk >= length(lv))
        return(as.integer(factor(rank(lv, ties.method = "min"))))
      ## WSS elbow fallback
      wss <- vapply(seq_len(kk), function(kc) {
        if (kc == 1L) return(sum((lv - mean(lv))^2))
        sum(km_seeded(lv, kc)$withinss)
      }, numeric(1))
      drops <- if (length(wss) > 1) -diff(wss) / pmax(wss[-length(wss)], .Machine$double.eps) else numeric(0)
      kbest <- 1L
      for (i in seq_along(drops)) if (drops[i] >= 0.15) kbest <- i + 1L
      cl <- if (kbest == 1L) rep(1L, length(lv)) else km_seeded(lv, kbest)$cluster
      return(as.integer(cl))
    }
    ki <- min(as.integer(k), uniq_n)
    if (ki <= 1L) return(rep(1L, length(vf)))
    if (ki >= length(lv)) return(as.integer(factor(rank(lv, ties.method = "min"))))
    as.integer(km_seeded(lv, ki)$cluster)
  }
  ## kmeans with a fixed seed that does not disturb the caller's RNG stream
  km_seeded <- function(x, centers) {
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(1L)
    stats::kmeans(x, centers = centers, nstart = 25L)
  }

  ## Founder isolation: the single highest-VF junction (the founder / defining
  ## excision) is pulled out of the clustering and always given the top panel,
  ## so it neither distorts the k-means solution nor the shared read-support
  ## scale. Ties at the maximum VF are treated as one founder group.
  sv$is_founder <- FALSE
  founder_mask <- rep(FALSE, nrow(sv))
  if (isolate_founder && length(unique(sv$VF)) > 1L) {
    founder_mask <- sv$VF == max(sv$VF, na.rm = TRUE)
    sv$is_founder <- founder_mask
  }
  rest <- which(!founder_mask)

  ## cluster only the non-founder junctions
  rest_raw <- assign_strata(sv$VF[rest])
  ord <- names(sort(tapply(sv$VF[rest], rest_raw, stats::median), decreasing = TRUE))
  rest_stratum <- match(as.character(rest_raw), ord)

  sv$stratum <- NA_integer_
  founder_offset <- if (any(founder_mask)) 1L else 0L
  if (any(founder_mask)) sv$stratum[founder_mask] <- 1L
  sv$stratum[rest] <- founder_offset + rest_stratum
  n_strata <- founder_offset + length(ord)

  strata_meta <- do.call(rbind, lapply(seq_len(n_strata), function(s) {
    v <- sv$VF[sv$stratum == s]
    data.frame(stratum = s, n = length(v),
               vf_min = min(v), vf_max = max(v), vf_med = stats::median(v),
               is_founder = any(sv$is_founder[sv$stratum == s]),
               stringsAsFactors = FALSE)
  }))
  if (verbose) {
    message("VF strata (highest read support first):")
    for (s in seq_len(n_strata))
      message(sprintf("  [%d]%s VF %.0f-%.0f (median %.0f): %d junctions",
                      s, if (strata_meta$is_founder[s]) " [founder]" else "",
                      strata_meta$vf_min[s], strata_meta$vf_max[s],
                      strata_meta$vf_med[s], strata_meta$n[s]))
  }

  ## ---- Global read-support -> CN-unit scaling (shared by all panels) ------
  ## By default the axis top reflects the true maximum VF (founder included), so
  ## the founder arc's height is honest. With `founder_offscale = TRUE` the axis
  ## is instead driven by the highest NON-founder junction (giving the lower
  ## strata the full dynamic range), and the founder arc is drawn clamped to the
  ## axis top with its true VF labelled -- useful when the founder is a large
  ## outlier that would otherwise compress every other stratum.
  max_vf <- max(sv$VF, na.rm = TRUE)
  max_vf_scale <- if (founder_offscale && length(rest)) max(sv$VF[rest], na.rm = TRUE) else max_vf
  rs_axis_max <- if (max_vf_scale > 0) nice_ceiling(max_vf_scale) else 0
  coeff <- if (rs_axis_max > 0) rs_axis_max / cn_axis_max else 1

  ## When the founder is clamped off-scale the whole figure's read-support axis
  ## no longer maps 1:1 to VF, so mark EVERY panel's axis title "(scaled)".
  founder_vf <- if (any(sv$is_founder)) max(sv$VF[sv$is_founder], na.rm = TRUE) else 0
  rs_axis_name <- if (isTRUE(founder_offscale) && founder_vf > rs_axis_max)
    "Read support (scaled)" else "Read support"

  ## ---- Shared colours / sizes / geometry (mirrors plot_sv_linear) ---------
  color_minor_cn <- "#3a9387"; color_major_cn <- "#d92a05"
  color_homdel <- "#0077b6"; color_loh <- "#bde0fe"
  cn_size <- 0.9; size_sv_line <- 0.2; size_interchr_line <- 0.3
  founder_lwd_mult <- 1.6   # founder junction drawn slightly bolder than the rest
  size_gene_label <- 3; size_text <- 8
  curv_intra <- 0.18; curv_inter <- -0.18

  if (max.cn > 250) { minorAllele_offset <- 5; upper_limit_karyotype <- -max.cn * 0.14
  } else if (max.cn > 110) { minorAllele_offset <- 3; upper_limit_karyotype <- -max.cn * 0.16
  } else if (max.cn > 40) { minorAllele_offset <- 1; upper_limit_karyotype <- -max.cn * 0.18
  } else if (max.cn < 6) { minorAllele_offset <- 0.05; upper_limit_karyotype <- -max.cn * 0.25
  } else { minorAllele_offset <- 0.1; upper_limit_karyotype <- -max.cn * 0.10 }
  yend_outside_range <- max.cn * 0.03
  lower_limit_karyotype <- upper_limit_karyotype - (cn_axis_max * karyotype_rel_size)
  loh_bar_y <- upper_limit_karyotype * loh_position_ratio

  ## shared y-limits so panels align and arc heights compare 1:1
  y_lo <- lower_limit_karyotype * 1.05
  y_hi <- max.cn * offset_gene * 1.22

  ## bottom-axis label geometry (only drawn on the bottom panel)
  unit_y      <- cn_axis_max * 0.07
  tick_y0     <- lower_limit_karyotype
  tick_y1     <- lower_limit_karyotype - 0.7 * unit_y
  mb_label_y  <- lower_limit_karyotype - 1.6 * unit_y
  chr_label_y <- lower_limit_karyotype - 3.0 * unit_y

  plot_width <- if (!is.null(plot_width_custom)) plot_width_custom else max(5, 2.2 * nrow(loci) + 1)
  min_tick_gap <- 0.62 * total_w / plot_width
  tick_df <- data.frame(); mb_df <- data.frame()
  for (i in seq_len(nrow(loci))) {
    ticks <- pretty(c(loci$start[i], loci$end[i]), n = 3)
    ticks <- ticks[ticks >= loci$start[i] & ticks <= loci$end[i]]
    ticks <- ticks[ticks != 0]
    if (length(ticks) == 0) next
    gxt <- loci$offset[i] + (ticks - loci$start[i])
    keep <- rep(FALSE, length(gxt)); last <- -Inf
    for (j in seq_along(gxt)) if (gxt[j] - last >= min_tick_gap) { keep[j] <- TRUE; last <- gxt[j] }
    gxt <- gxt[keep]; ticks <- ticks[keep]
    tick_df <- rbind(tick_df, data.frame(x = gxt))
    mb_df <- rbind(mb_df, data.frame(x = gxt, lab = formatC(ticks / 1e6, format = "f", digits = 1)))
  }
  chr_lab <- data.frame(gx = loci$gx_mid, chr = gsub("chr", "Chr ", loci$chr))
  ideo <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    d <- karyotype[karyotype$chr == loci$chr[i] & karyotype$end >= loci$start[i] & karyotype$start <= loci$end[i], , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$start <- pmax(d$start, loci$start[i]); d$end <- pmin(d$end, loci$end[i])
    d$gx_start <- loci$offset[i] + (d$start - loci$start[i])
    d$gx_end   <- loci$offset[i] + (d$end   - loci$start[i])
    d
  }))
  cn_breaks <- unique(c(0, cn_axis_max / 2, cn_axis_max))
  xlim <- c(loci$gx_start[1], loci$gx_end[nrow(loci)])

  ## which CN segments lie near a stratum's breakpoints
  near_cn <- function(sv_sub) {
    keep <- rep(FALSE, nrow(cn_plot))
    for (i in seq_len(nrow(sv_sub))) {
      ends <- list()
      if (!is.na(sv_sub$l1[i])) ends[[length(ends) + 1L]] <- list(li = sv_sub$l1[i], pos = sv_sub$pos1[i])
      if (!is.na(sv_sub$l2[i])) ends[[length(ends) + 1L]] <- list(li = sv_sub$l2[i], pos = sv_sub$pos2[i])
      for (e in ends) {
        keep <- keep | (cn_plot$li == e$li &
                        cn_plot$end   >= e$pos - cn_near_flank &
                        cn_plot$start <= e$pos + cn_near_flank)
      }
      ## intra-locus junction: also the segments it bounds
      if (!is.na(sv_sub$l1[i]) && !is.na(sv_sub$l2[i]) && sv_sub$l1[i] == sv_sub$l2[i]) {
        lo <- min(sv_sub$pos1[i], sv_sub$pos2[i]); hi <- max(sv_sub$pos1[i], sv_sub$pos2[i])
        keep <- keep | (cn_plot$li == sv_sub$l1[i] & cn_plot$end >= lo & cn_plot$start <= hi)
      }
    }
    cn_plot[keep, , drop = FALSE]
  }

  ## ---- Iterative CN reconstruction (cn_display = "reconstruct") ------------
  ## The breakpoints of the junctions introduced so far (strata 1..k) partition
  ## each locus into intervals; every interval is drawn flat at the MINIMUM of
  ## the two copy-number segments immediately adjacent to (just inside) its
  ## bounding breakpoints. The founder alone gives one flat baseline over its
  ## span; each lower panel introduces more breakpoints, so the intervals
  ## subdivide and copy-number structure emerges; the bottom panel is the full
  ## observed bulk copy number (the reconstruction target). A visual model of
  ## how the amplicon is built up wave by wave, not a formal deconvolution.
  reconstruct_cn <- function(max_stratum) {
    ## bottom panel: the full observed bulk CN that the reconstruction evolves to
    if (max_stratum >= n_strata) {
      return(data.frame(li = cn_plot$li, start = cn_plot$start, end = cn_plot$end,
                        y = cn_plot$majorAlleleCopyNumber,
                        gx_start = cn_plot$gx_start, gx_end = cn_plot$gx_end))
    }
    jj <- sv[sv$stratum <= max_stratum, , drop = FALSE]
    out <- list()
    for (li in seq_len(nrow(loci))) {
      cn_li <- cn_plot[cn_plot$li == li, , drop = FALSE]
      if (nrow(cn_li) == 0) next
      cn_li <- cn_li[order(cn_li$start), , drop = FALSE]
      ## breakpoints introduced in this locus so far, plus the locus edges
      bps <- c(loci$start[li], loci$end[li],
               jj$pos1[which(jj$l1 == li)], jj$pos2[which(jj$l2 == li)])
      bps <- sort(unique(pmin(pmax(bps, loci$start[li]), loci$end[li])))
      if (length(bps) < 2) next
      for (i in seq_len(length(bps) - 1L)) {
        a <- bps[i]; b <- bps[i + 1L]; if (b <= a) next
        ## the CN segments immediately inside the interval's edges, found by
        ## segment boundary (robust to the few-bp offset between an SV breakpoint
        ## and the copy-number segment boundary -- an epsilon-position lookup can
        ## otherwise land in the flanking segment and collapse the level).
        left  <- cn_li[cn_li$start >= a & cn_li$start < b, , drop = FALSE]
        right <- cn_li[cn_li$end   <= b & cn_li$end   > a, , drop = FALSE]
        cand <- c(if (nrow(left))  left$majorAlleleCopyNumber[which.min(left$start)],
                  if (nrow(right)) right$majorAlleleCopyNumber[which.max(right$end)])
        if (length(cand) == 0) {                         # interval sits within one segment
          s <- cn_li[cn_li$start <= a & cn_li$end >= b, , drop = FALSE]
          if (nrow(s)) cand <- s$majorAlleleCopyNumber[1]
        }
        if (length(cand) == 0) next
        out[[length(out) + 1L]] <- data.frame(li = li, start = a, end = b, y = min(cand))
      }
    }
    if (length(out) == 0) return(NULL)
    res <- as.data.frame(data.table::rbindlist(out))
    res$gx_start <- gx_of(res$li, res$start)
    res$gx_end   <- gx_of(res$li, res$end)
    res
  }

  ## Borders of the final (fully reconstructed) CN segments -> vertical dashed
  ## guides repeated in every panel (shared x -> they align across the stack).
  ## Each locus's outer edges plus every internal CN step larger than
  ## `border_step` copies are drawn, so a guide lands on essentially every
  ## copy-number transition (flat equal-level splits and sub-pixel duplicates are
  ## the only ones dropped).
  border_step <- if (!is.null(cn_border_min_step)) cn_border_min_step else max(0.25, 0.015 * cn_axis_max)
  cn_border_gx <- numeric(0)
  if (cn_border_lines) {
    recF <- reconstruct_cn(n_strata)
    if (!is.null(recF) && nrow(recF) > 0) {
      for (lc in unique(recF$li)) {
        d <- recF[recF$li == lc, , drop = FALSE]
        d <- d[order(d$start), , drop = FALSE]
        cn_border_gx <- c(cn_border_gx, d$gx_start[1], d$gx_end[nrow(d)])   # outer edges
        if (nrow(d) > 1) {
          chg <- which(abs(d$y[-1] - d$y[-nrow(d)]) > border_step)          # major level changes
          if (length(chg)) cn_border_gx <- c(cn_border_gx, d$gx_end[chg])
        }
      }
      cn_border_gx <- sort(unique(round(cn_border_gx, 3)))
      ## collapse only borders that are essentially coincident (< ~0.3% of the
      ## plotted width, i.e. within a line-width of each other) so every distinct
      ## copy-number transition keeps its own guide, while sub-pixel duplicates
      ## from the segmentation do not stack into a single fat line.
      if (length(cn_border_gx) > 1) {
        min_gap <- total_w * 0.003
        keep <- rep(TRUE, length(cn_border_gx)); last <- cn_border_gx[1]
        for (j in 2:length(cn_border_gx)) {
          if (cn_border_gx[j] - last < min_gap) keep[j] <- FALSE else last <- cn_border_gx[j]
        }
        cn_border_gx <- cn_border_gx[keep]
      }
    }
  }

  ## ---- LOH / hom-del reconstruction ---------------------------------------
  ## Place allele-loss segments on the same timeline as the major allele: a
  ## segment is "introduced" by the earliest (highest-VF = lowest stratum) intra-
  ## locus deletion junction whose deleted span covers the segment midpoint; a
  ## segment covered by no deletion is treated as pre-existing (assigned to the
  ## founder wave, i.e. present from the top panel down). In `reconstruct` each
  ## bar is then revealed only in the wave that creates it and every wave below,
  ## so e.g. a large deletion and the LOH it causes appear in the same panel.
  loh_seg <- cn_plot[cn_plot$minorAlleleCopyNumber < loh_thresh    & cn_plot$end > cn_plot$start, , drop = FALSE]
  hd_seg  <- cn_plot[cn_plot$copyNumber            < homdel_thresh & cn_plot$end > cn_plot$start, , drop = FALSE]
  del_j <- sv[(sv$svclass == "DEL" | sv$strands %in% c("+-", "DEL")) &
              !is.na(sv$l1) & !is.na(sv$l2) & sv$l1 == sv$l2, , drop = FALSE]
  wave_of_seg <- function(seg) {
    if (nrow(seg) == 0) return(integer(0))
    w <- rep(1L, nrow(seg))                        # default: pre-founder background
    if (nrow(del_j) > 0) {
      lo <- pmin(del_j$pos1, del_j$pos2); hi <- pmax(del_j$pos1, del_j$pos2)
      mid <- (seg$start + seg$end) / 2
      for (i in seq_len(nrow(seg))) {
        hit <- del_j$l1 == seg$li[i] & lo <= mid[i] & hi >= mid[i]
        if (any(hit)) w[i] <- min(del_j$stratum[hit])
      }
    }
    w
  }
  loh_seg$wave <- wave_of_seg(loh_seg)
  hd_seg$wave  <- wave_of_seg(hd_seg)
  if (verbose) {
    message("Deletion junctions per stratum: ",
            paste(sprintf("s%d=%d", sort(unique(del_j$stratum)),
                          as.integer(table(del_j$stratum)[as.character(sort(unique(del_j$stratum)))])),
                  collapse = ", "))
    message("LOH segments introduced per wave (1 = pre-founder): ",
            paste(sprintf("w%s=%d", names(table(loh_seg$wave)), as.integer(table(loh_seg$wave))), collapse = ", "))
    if (nrow(hd_seg) > 0)
      message("Hom-del segments introduced per wave: ",
              paste(sprintf("w%s=%d", names(table(hd_seg$wave)), as.integer(table(hd_seg$wave))), collapse = ", "))
  }

  ## small geom helpers shared by the reconstruct / actual panel branches: a
  ## full-height highlight rectangle, and an allele-loss bar at `loh_bar_y`. Both
  ## no-op on an empty data frame so call sites stay branch-free.
  add_shade <- function(p, df, fill, alpha) if (nrow(df) > 0)
    p + geom_rect(data = df, aes(xmin = gx_start, xmax = gx_end,
        ymin = 0, ymax = max.cn * ymax_highlight_ratio), fill = fill, alpha = alpha) else p
  add_bar <- function(p, df, colour) if (nrow(df) > 0)
    p + geom_segment(data = df, aes(x = gx_start, xend = gx_end,
        y = loh_bar_y, yend = loh_bar_y), colour = colour, linewidth = cn_size) else p

  ## ---- Per-stratum panel builder ------------------------------------------
  build_panel <- function(s, is_bottom, is_top) {
    sv_s <- sv[sv$stratum == s, , drop = FALSE]
    cn_s <- near_cn(sv_s)
    panel_vf_max <- if (nrow(sv_s) > 0) max(sv_s$VF, na.rm = TRUE) else 0

    ## Read-support -> CN-unit scaling. "shared" (default): one global scale so
    ## arc heights compare across strata. "per_panel": each panel fills its own
    ## VF range. The left CN axis and x-axis stay shared, so panels remain aligned.
    if (vf_scale == "per_panel" && panel_vf_max > 0) {
      rs_max_s <- nice_ceiling(panel_vf_max)
      coeff_s  <- rs_max_s / cn_axis_max
    } else {
      rs_max_s <- rs_axis_max
      coeff_s  <- coeff
    }
    ## `rs_axis_name` ("Read support" or "Read support (scaled)") is set once for
    ## the whole figure, so every panel's axis title reads the same.
    rs_name <- rs_axis_name

    p <- ggplot()

    ## karyotype ideogram + bottom labels only on the bottom panel
    if (is_bottom && !is.null(ideo) && nrow(ideo) > 0)
      p <- p + geom_rect(data = ideo, aes(xmin = gx_start, xmax = gx_end,
               ymin = lower_limit_karyotype, ymax = upper_limit_karyotype),
               fill = ideo$color, colour = "black", linewidth = 0.2)

    ## vertical dashed CN-border guides (behind the data). Each panel extends its
    ## lines beyond its own plotting area (clip = "off") so they bridge the gaps
    ## between panels and read as continuous verticals down the whole figure; the
    ## figure's outer top and bottom are left clean.
    if (length(cn_border_gx) > 0) {
      ext <- (y_hi - y_lo) * 0.75
      yb <- if (is_bottom) y_lo else y_lo - ext
      yt <- if (is_top) cn_axis_max else y_hi + ext
      p <- p + geom_segment(data = data.frame(x = cn_border_gx),
               aes(x = x, xend = x, y = yb, yend = yt),
               linetype = "dashed", colour = "grey78", linewidth = 0.2, alpha = 0.6)
    }

    ## Major-allele CN: either the cumulative reconstruction through this wave
    ## ("reconstruct") or the observed near-breakpoint segments ("actual").
    rec <- if (cn_display == "reconstruct") reconstruct_cn(s) else NULL

    ## amplification / hom-del shading
    if (cn_display == "reconstruct") {
      if (highlight_amp && !is.null(rec)) p <- add_shade(p, rec[rec$y > 3 * ploidy, , drop = FALSE], "#d92a05", 0.1)
      if (highlight_hom_del)             p <- add_shade(p, hd_seg[hd_seg$wave <= s, , drop = FALSE], color_homdel, 0.05)
    } else if (nrow(cn_s) > 0) {
      if (highlight_amp)     p <- add_shade(p, cn_s[cn_s$copyNumber > 3 * cn_s$ploidy, , drop = FALSE], "#d92a05", 0.1)
      if (highlight_hom_del) p <- add_shade(p, cn_s[cn_s$copyNumber < homdel_thresh, , drop = FALSE], color_homdel, 0.05)
    }

    ## Minor-allele CN is drawn across the FULL locus in every panel -- it is the
    ## stable background allele and serves as a constant reference, so it must not
    ## disappear when the near-breakpoint filter trims the (amplified) major
    ## allele in the lower strata.
    p <- p +
      geom_segment(data = cn_plot, aes(x = gx_start, xend = gx_end,
                   y = minorAlleleCopyNumber - minorAllele_offset,
                   yend = minorAlleleCopyNumber - minorAllele_offset),
                   colour = color_minor_cn, linewidth = cn_size)

    ## Major-allele CN, then allele-loss bars. In reconstruct the bars are
    ## revealed by the wave that introduces them (or, for pre-existing loss, from
    ## the founder wave down); in actual they are the observed near-breakpoint
    ## LOH / hom-del segments.
    if (cn_display == "reconstruct") {
      if (!is.null(rec) && nrow(rec) > 0)
        p <- p + geom_segment(data = rec, aes(x = gx_start, xend = gx_end, y = y, yend = y),
                              colour = color_major_cn, linewidth = cn_size)
      p <- add_bar(p, loh_seg[loh_seg$wave <= s, , drop = FALSE], color_loh)
      p <- add_bar(p, hd_seg[hd_seg$wave  <= s, , drop = FALSE], color_homdel)
    } else if (nrow(cn_s) > 0) {
      p <- p +
        geom_segment(data = cn_s, aes(x = gx_start, xend = gx_end,
                     y = majorAlleleCopyNumber, yend = majorAlleleCopyNumber),
                     colour = color_major_cn, linewidth = cn_size)
      p <- add_bar(p, cn_s[cn_s$minorAlleleCopyNumber < loh_thresh,    , drop = FALSE], color_loh)
      p <- add_bar(p, cn_s[cn_s$copyNumber            < homdel_thresh, , drop = FALSE], color_homdel)
      amp_seg <- cn_s[cn_s$copyNumber > 3 * cn_s$ploidy, , drop = FALSE]
      if (nrow(amp_seg) > 0) p <- p + geom_segment(data = amp_seg,
        aes(x = gx_start, xend = gx_end, y = majorAlleleCopyNumber, yend = majorAlleleCopyNumber),
        colour = color_major_cn, linewidth = cn_size - 0.3)
    }

    ## SV arcs + breakpoint drop-lines (arc height = VF / coeff). Earlier
    ## (higher-VF) strata are re-drawn faintly in each lower panel so the
    ## accumulating rearrangements stay visible; the current stratum is at full
    ## opacity. Every junction is drawn at the same weight (including the max-VF
    ## one); the max-VF junction is distinguished only by having its own top panel.
    y_cap <- cn_axis_max
    add_sv_layer <- function(p, svset, alpha) {
      if (nrow(svset) == 0) return(p)
      arc_rows <- list(); seg_rows <- list()
      for (i in seq_len(nrow(svset))) {
        yv <- min(svset$VF[i] / coeff_s, y_cap); col <- svset$colour[i]
        ## the founder junction is drawn with a slightly bolder line so it stays
        ## legible even where it spans the whole amplicon or is re-drawn faintly.
        fmul <- if (isTRUE(svset$is_founder[i])) founder_lwd_mult else 1
        if (isTRUE(svset$hl[i])) fmul <- max(fmul, highlight_lwd_mult)   # highlighted events bolder
        lwd_seg <- size_sv_line * fmul
        in1 <- !is.na(svset$l1[i]); in2 <- !is.na(svset$l2[i])
        gx1 <- if (in1) gx_of(svset$l1[i], svset$pos1[i]) else NA
        gx2 <- if (in2) gx_of(svset$l2[i], svset$pos2[i]) else NA
        if (in1 && in2) {
          same <- svset$l1[i] == svset$l2[i]
          cv <- if (same) {
            base <- if (abs(gx2 - gx1) <= total_w * 0.15) curv_intra * 1.4 else curv_intra
            if (svset$strands[i] %in% c("DEL", "h2hINV", "+-", "--")) base else -base
          } else curv_inter
          arc_rows[[length(arc_rows) + 1L]] <- data.frame(x = gx1, xend = gx2, y = yv, yend = yv,
            curvature = cv, colour = col,
            lwd = (if (same) size_sv_line else size_interchr_line) * fmul)
          seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx1, xend = gx1, y = 0, yend = yv, colour = col, lwd = lwd_seg)
          seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx2, xend = gx2, y = 0, yend = yv, colour = col, lwd = lwd_seg)
        } else {
          gx <- if (in1) gx1 else gx2
          seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx, xend = gx, y = 0, yend = yv, colour = col, lwd = lwd_seg)
          seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx, xend = gx - total_w * 0.006, y = 0, yend = -yend_outside_range, colour = col, lwd = lwd_seg)
        }
      }
      if (length(arc_rows) > 0) {
        arc_df <- as.data.frame(data.table::rbindlist(arc_rows))
        arc_df$grp <- paste(arc_df$curvature, arc_df$lwd, sep = "_")
        for (g in unique(arc_df$grp)) {
          d <- arc_df[arc_df$grp == g, , drop = FALSE]
          p <- p + geom_curve(data = d, aes(x = x, xend = xend, y = y, yend = yend, colour = colour),
                              curvature = d$curvature[1], linewidth = d$lwd[1], alpha = alpha)
        }
      }
      if (length(seg_rows) > 0) {
        seg_df <- as.data.frame(data.table::rbindlist(seg_rows))
        for (w in unique(seg_df$lwd)) {
          d <- seg_df[seg_df$lwd == w, , drop = FALSE]
          p <- p + geom_segment(data = d, aes(x = x, xend = xend, y = y, yend = yend, colour = colour),
                                linewidth = w, alpha = alpha)
        }
      }
      p
    }
    ## Prior (already-placed) junctions are faded, but the founder is held at a
    ## higher opacity so it can still be followed all the way down the panels.
    prior <- sv[sv$stratum < s, , drop = FALSE]
    prior_founder <- prior[prior$is_founder, , drop = FALSE]
    prior_rest    <- prior[!prior$is_founder, , drop = FALSE]
    if (prior_sv_alpha > 0 && nrow(prior_rest) > 0)
      p <- add_sv_layer(p, prior_rest, alpha = prior_sv_alpha)
    if (nrow(prior_founder) > 0)
      p <- add_sv_layer(p, prior_founder, alpha = max(prior_sv_alpha, founder_alpha))
    p <- add_sv_layer(p, sv_s, alpha = 1)
    p <- p + scale_colour_identity()

    ## Exon models (top panel): draw each gene's exons as boxes on a body line;
    ## exons that fall entirely within a deletion junction are coloured red, so
    ## an intragenic loss (e.g. EGFR exons 2-7 = EGFRvIII) is immediately visible.
    if (is_top && displayExon && !is.null(cds_gr)) {
      cds <- as.data.frame(cds_gr)
      cds$chr <- if (all(grepl("^chr", cds$seqnames))) as.character(cds$seqnames) else paste0("chr", cds$seqnames)
      yb <- max.cn * offset_gene
      min_ex_w <- total_w * 0.004
      dels <- sv[sv$svclass == "DEL" | sv$strands %in% c("+-", "DEL"), , drop = FALSE]
      ex_df <- data.frame(); body_df <- data.frame(); lab_df <- data.frame()
      for (i in seq_len(nrow(loci))) {
        for (gn in intersect(genes, cds$gene_name[cds$chr == loci$chr[i]])) {
          ex <- cds[cds$gene_name == gn & cds$chr == loci$chr[i] &
                    cds$end >= loci$start[i] & cds$start <= loci$end[i], , drop = FALSE]
          if (nrow(ex) == 0) next
          del <- rep(FALSE, nrow(ex))
          if (nrow(dels) > 0) {
            onchr <- dels$chr1 == loci$chr[i] & dels$chr2 == loci$chr[i]
            lo <- pmin(dels$pos1, dels$pos2); hi <- pmax(dels$pos1, dels$pos2)
            for (k in seq_len(nrow(ex)))
              del[k] <- any(onchr & lo <= ex$start[k] & hi >= ex$end[k])
          }
          gxs <- loci$offset[i] + (pmax(ex$start, loci$start[i]) - loci$start[i])
          gxe <- loci$offset[i] + (pmin(ex$end, loci$end[i]) - loci$start[i])
          short <- (gxe - gxs) < min_ex_w; mid <- (gxs + gxe) / 2
          gxs[short] <- mid[short] - min_ex_w / 2; gxe[short] <- mid[short] + min_ex_w / 2
          ex_df <- rbind(ex_df, data.frame(gxs = gxs, gxe = gxe, deleted = del))
          body_df <- rbind(body_df, data.frame(gxs = min(gxs), gxe = max(gxe)))
          lab_df <- rbind(lab_df, data.frame(gx = (min(gxs) + max(gxe)) / 2, gene = gn))
        }
      }
      if (nrow(ex_df) > 0) {
        p <- p + geom_segment(data = body_df, aes(x = gxs, xend = gxe, y = yb * 0.945, yend = yb * 0.945),
                              colour = "grey40", linewidth = 0.4)
        en <- ex_df[!ex_df$deleted, , drop = FALSE]; ed <- ex_df[ex_df$deleted, , drop = FALSE]
        if (nrow(en) > 0) p <- p + geom_rect(data = en, aes(xmin = gxs, xmax = gxe, ymin = yb * 0.90, ymax = yb * 0.99),
                              fill = "grey25", colour = "grey15", linewidth = 0.1)
        if (nrow(ed) > 0) p <- p + geom_rect(data = ed, aes(xmin = gxs, xmax = gxe, ymin = yb * 0.90, ymax = yb * 0.99),
                              fill = "#d92a05", colour = "grey15", linewidth = 0.1)
        p <- p + geom_text(data = lab_df, aes(x = gx, y = yb * 1.08, label = gene),
                           size = size_gene_label, fontface = "italic")
      }
    } else if (is_top) {
      gl <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
        g <- gene_coord[gene_coord$chr == loci$chr[i] &
                        gene_coord$start >= loci$start[i] & gene_coord$end <= loci$end[i] &
                        gene_coord$gene %in% genes, , drop = FALSE]
        if (nrow(g) == 0) return(NULL)
        g$gx <- loci$offset[i] + ((g$start + g$end) / 2 - loci$start[i])
        g
      }))
      if (!is.null(gl) && nrow(gl) > 0) {
        gl$y <- max.cn * offset_gene
        ang <- gene_label_angle
        if (is.null(ang)) {
          gp <- sort(gl$gx); ang <- if (length(gp) > 1 && min(diff(gp)) < 0.05 * total_w) 45 else 0
        }
        p <- p + geom_point(data = gl, aes(x = gx, y = y * 0.94), shape = 16, size = 1, colour = "black")
        if (repel_labels && requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(data = gl, aes(x = gx, y = y, label = gene),
            size = size_gene_label, fontface = "italic", angle = ang, hjust = 0,
            direction = "both", nudge_y = max.cn * 0.12, ylim = c(max.cn * 1.02, NA),
            segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0,
            box.padding = 0.25, max.overlaps = Inf, seed = 1L)
        } else {
          p <- p + geom_text(data = gl, aes(x = gx, y = y, label = gene),
                             size = size_gene_label, fontface = "italic", angle = ang)
        }
      }
    }

    ## CN axis frame (left/right verticals), bottom labels on bottom panel only
    p <- p +
      geom_segment(data = data.frame(x = c(loci$gx_start[1], loci$gx_end[nrow(loci)])),
                   aes(x = x, xend = x, y = 0, yend = cn_axis_max), colour = "black", linewidth = 0.4)
    if (is_bottom) {
      p <- p +
        geom_segment(data = tick_df, aes(x = x, xend = x, y = tick_y0, yend = tick_y1),
                     colour = "black", linewidth = 0.4) +
        geom_text(data = mb_df, aes(x = x, y = mb_label_y, label = lab), size = size_text / 2.1, vjust = 1) +
        geom_text(data = chr_lab, aes(x = gx, y = chr_label_y, label = chr),
                  size = size_text / 1.78, vjust = 1, fontface = "bold")
    }

    ## per-panel read-support stratum header (sits above the panel, so it never
    ## collides with in-panel gene labels or arcs). The isolated max-VF junction
    ## reads "Max VF"; the remaining VF clusters are numbered 1..X (highest first).
    tag <- if (isTRUE(strata_meta$is_founder[s]))
      sprintf("Max VF: %.0f  (n = %d)", strata_meta$vf_max[s], strata_meta$n[s])
    else
      sprintf("Cluster %d: VF %.0f-%.0f  (n = %d)", s - founder_offset,
              strata_meta$vf_min[s], strata_meta$vf_max[s], strata_meta$n[s])

    p <- p +
      ggtitle(tag) +
      labs(x = NULL, y = "Allele specific\ncopy number") +
      coord_cartesian(xlim = xlim, ylim = c(y_lo, y_hi), clip = "off", expand = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.01))
    if (max_vf > 0) {
      rs_breaks_s <- unique(c(0, rs_max_s / 2, rs_max_s))
      p <- p + scale_y_continuous(breaks = cn_breaks,
        sec.axis = sec_axis(trans = ~ . * coeff_s, breaks = rs_breaks_s, name = rs_name))
    } else {
      p <- p + scale_y_continuous(breaks = cn_breaks)
    }
    p <- p + theme(
      text = element_text(size = size_text, colour = "black"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_text(size = size_text + 1, colour = "black"),
      axis.text.y.right = element_text(size = size_text + 1, colour = "black"),
      axis.title.y = element_text(size = size_text + 1, colour = "black"),
      axis.title.y.right = element_text(size = size_text + 1, colour = "black"),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y.right = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length.y = unit(0.2, "cm"),
      panel.background = element_blank(), plot.background = element_blank(), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0, size = size_text, face = "bold", colour = "grey20",
                                margin = margin(b = 1)),
      plot.margin = unit(c(if (is_top) .15 else .12, .5, if (is_bottom) 1.3 else .08, .2), "cm")
    )
    p
  }

  panels <- lapply(seq_len(n_strata), function(s)
    build_panel(s, is_bottom = (s == n_strata), is_top = (s == 1L)))

  heights <- c(rep(1, max(0, n_strata - 1)), panel_rel_height)[seq_len(n_strata)]
  title_txt <- if (is.null(wgd_status)) this_sample else paste0(this_sample, " (", wgd_status, ")")
  combined <- patchwork::wrap_plots(panels, ncol = 1, heights = heights) +
    patchwork::plot_annotation(title = title_txt,
      theme = theme(plot.title = element_text(hjust = 0.5, size = size_text + 4)))

  ## ---- Size + save --------------------------------------------------------
  per_panel_h <- 1.6
  plot_height <- if (!is.null(plot_height_custom)) plot_height_custom else
    per_panel_h * n_strata + 0.9

  outfile <- NULL
  if (!is.null(outdir) && save) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    region_tag <- if (nrow(loci) == 1)
      sprintf("_%.0f-%.0f", round(loci$start[1]), round(loci$end[1])) else ""
    outfile <- file.path(outdir, paste0(this_sample, "_", paste(loci$chr, collapse = "_"),
                                        region_tag, "_by_VF.pdf"))
    grDevices::pdf(outfile, width = plot_width, height = plot_height, useDingbats = FALSE)
    print(combined); grDevices::dev.off()
    if (verbose) message("Wrote ", outfile)
  }
  attr(combined, "path") <- outfile
  attr(combined, "strata") <- sv[, c("chr1", "pos1", "chr2", "pos2", "svclass", "VF", "stratum", "is_founder")]
  if (!is.null(outfile)) invisible(combined) else combined
}
