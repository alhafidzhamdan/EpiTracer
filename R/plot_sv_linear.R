#' Linear copy-number and structural-variant "recon" plot
#'
#' Draws allele-specific copy number, structural-variant arcs, karyotype
#' ideograms and gene labels for one or more loci laid out side-by-side on a
#' single concatenated x-axis. Because all loci share one coordinate system, the
#' structural variants that interconnect separate amplicons (e.g. the junctions
#' of a multi-fragment / hub ecDNA) are drawn as arcs spanning the loci.
#'
#' The loci to display are resolved in this order:
#' \enumerate{
#'   \item `loci` if supplied (explicit windows);
#'   \item `chromosome` + `chromosome_range` (explicit windows, one per
#'     chromosome);
#'   \item `chromosome` alone -- the amplified region on each named chromosome is
#'     auto-detected and padded by `margin`; a chromosome with no amplification
#'     falls back to its whole length;
#'   \item nothing -- every amplified locus in the sample is detected
#'     automatically.
#' }
#' A locus is "amplified" where `copyNumber > min_cn_ratio * ploidy`; its
#' extent is the min-max of amplified segments on the chromosome, padded by
#' `margin`, with single-segment artefacts below `min_amp_width` dropped.
#'
#' @param sample Character scalar; sample identifier (matched in `cnv_data`,
#'   `sv_data`, `wgd_data`).
#' @param cnv_data A data.frame of copy-number segments with columns `sample`,
#'   `seqnames`, `start`, `end`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`,
#'   `minorAlleleCopyNumber`.
#' @param sv_data A data.frame of SVs with columns `chrom1`, `start1`, `chrom2`,
#'   `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`.
#' @param wgd_data A data.frame with a sample-identifier column (see
#'   `wgd_sample_col`) and a `Polyploidy` column (`"No"` = diploid, otherwise
#'   WGD) used to label the plot title.
#' @param karyotype A data.frame of ideogram bands (UCSC `cytoBand` / `chr_info`
#'   style) with columns including chromosome, `start`, `end`, and `gieStain`;
#'   the first column is treated as the chromosome name. A path to an `.rds`
#'   holding such a data.frame is also accepted.
#' @param gene_coord A data.frame of gene coordinates with columns
#'   `chr`,`start`,`end`,`strand`,`gene`. A path to a headerless tab-separated
#'   BED-like file with those five columns is also accepted.
#' @param chromosome Optional character vector of chromosomes to display, e.g.
#'   `c("chr7", "chr12")`.
#' @param chromosome_range Optional two-column matrix/data.frame of `start`,`end`
#'   window limits, one row per entry in `chromosome`.
#' @param loci Optional explicit loci: a `data.frame` with columns
#'   `chr`,`start`,`end`, or a character vector of `"chr:start-end"` strings.
#'   Takes precedence over `chromosome`/`chromosome_range`.
#' @param margin Numeric; fraction of each auto-detected amplicon's width to pad
#'   on both sides (default `0.15`).
#' @param min_cn_ratio Numeric; a segment is "amplified" when
#'   `copyNumber > min_cn_ratio * ploidy` (default `3`).
#' @param min_amp_width Numeric; drop auto-detected loci whose total amplified
#'   span is below this many bp (default `1e5`).
#' @param gap_frac Numeric; gap between loci as a fraction of the total plotted
#'   width (default `0.06`).
#' @param genes_to_highlight Optional character vector of gene symbols. If `NULL`,
#'   a default oncogene panel is used; only genes falling inside a locus window
#'   are drawn.
#' @param gene_label_angle Optional numeric label rotation in degrees. If `NULL`
#'   (default) labels are horizontal, switching to 45 degrees automatically when
#'   genes are crowded.
#' @param repel_labels Logical; use \pkg{ggrepel} to de-collide gene labels
#'   (default `TRUE`; falls back to plain labels if ggrepel is absent).
#' @param cn_max Optional numeric override for the copy-number axis top; if
#'   `NULL` it is rounded up from the data to two significant figures.
#' @param offset_gene,ymax_highlight_ratio,karyotype_rel_size,loh_position_ratio
#'   Layout tuning parameters.
#' @param highlight_amp,highlight_hom_del Logical; shade amplified / homozygously
#'   deleted segments.
#' @param wgd_sample_col Optional name of the sample column in `wgd_data`
#'   (default: `sample`, falling back to `WGS_ID`).
#' @param outdir Optional directory in which to write the plot. If `NULL` no file
#'   is written and only the plot object is returned.
#' @param save Logical; if `TRUE` (default) and `outdir` is supplied, write the
#'   plot.
#' @param format Character vector of output formats: any of `"pdf"` and `"png"`
#'   (default both). The PNG is rendered from the plot object at `dpi`.
#' @param dpi Numeric PNG resolution (default `300`).
#' @param plot_width_custom,plot_height_custom Optional numeric overrides for the
#'   output dimensions (inches).
#' @param verbose Logical; print progress/diagnostic messages.
#'
#' @return The assembled [ggplot2::ggplot] object (invisibly when a file is
#'   written); written paths are attached as attribute `"path"`.
#'
#' @examples
#' \dontrun{
#' # Single focused locus:
#' plot_sv_linear("DO11441T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
#'                chromosome = "chr7",
#'                chromosome_range = matrix(c(52e6, 56e6), nrow = 1),
#'                outdir = "plots")
#'
#' # All amplified loci, interconnected (multi-fragment ecDNA hub):
#' plot_sv_linear("DUMC12T1", cnv, sv, wgd, karyotype = K, gene_coord = G,
#'                outdir = "plots")
#'
#' # Explicit loci:
#' plot_sv_linear("S1", cnv, sv, wgd, karyotype = K, gene_coord = G,
#'                loci = c("chr4:50e6-64e6", "chr12:57e6-59e6"))
#' }
#' @seealso [call_episomal_ecdna()]
#' @export
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom data.table rbindlist
#' @importFrom grDevices pdf dev.off
#' @importFrom utils read.table head
plot_sv_linear <- function(sample,
                           cnv_data,
                           sv_data,
                           wgd_data,
                           karyotype,
                           gene_coord,
                           chromosome = NULL,
                           chromosome_range = NULL,
                           loci = NULL,
                           margin = 0.15,
                           min_cn_ratio = 3,
                           min_amp_width = 1e5,
                           gap_frac = 0.06,
                           genes_to_highlight = NULL,
                           gene_label_angle = NULL,
                           repel_labels = TRUE,
                           cn_max = NULL,
                           offset_gene = 1.15,
                           ymax_highlight_ratio = 1.08,
                           karyotype_rel_size = 0.06,
                           loh_position_ratio = 0.5,
                           highlight_amp = TRUE,
                           highlight_hom_del = TRUE,
                           wgd_sample_col = NULL,
                           outdir = NULL,
                           save = TRUE,
                           format = c("pdf", "png"),
                           dpi = 300,
                           plot_width_custom = NULL,
                           plot_height_custom = NULL,
                           verbose = FALSE) {

  this_sample <- sample

  ## ---- Reference data -----------------------------------------------------
  if (is.character(karyotype)) karyotype <- readRDS(karyotype)
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
    gene_coord <- utils::read.table(gene_coord, header = FALSE, sep = "\t")
  }
  gene_coord <- as.data.frame(gene_coord)
  names(gene_coord) <- c("chr", "start", "end", "strand", "gene")
  if (is.null(genes_to_highlight)) {
    genes <- c("EGFR", "MDM4", "CDK4", "MDM2", "CCND2", "MYC", "FOXO1", "MET",
               "PRDM2", "AKT3", "EXO1", "PABPC1", "NTRK1", "GMPS", "SOX2",
               "RECQL5", "H3F3B", "KLF4", "RAD23B", "TAL2", "CDK6", "HSP90AB1",
               "NFKBIE", "NF1", "TERT", "PTEN", "PDGFRA", "TP53", "PIK3CA",
               "CDKN2A")
  } else {
    genes <- genes_to_highlight
  }

  ## ---- Sample data --------------------------------------------------------
  cnv <- cnv_data[cnv_data$sample == this_sample, , drop = FALSE]
  cnv$chr <- as.character(cnv$seqnames)
  if (nrow(cnv) == 0) stop("No copy-number data for sample ", this_sample, ".")

  wgd_data <- as.data.frame(wgd_data)
  wgd_col <- if (!is.null(wgd_sample_col)) wgd_sample_col else
    if ("sample" %in% names(wgd_data)) "sample" else "WGS_ID"
  wrow <- wgd_data[wgd_data[[wgd_col]] == this_sample, , drop = FALSE]
  wgd_status <- if (nrow(wrow) >= 1 && wrow$Polyploidy[1] == "No") "Diploid" else "WGD"

  ## ---- Resolve loci -------------------------------------------------------
  detect_loci <- function(chroms = NULL) {
    ploidy <- cnv$ploidy[1]
    amp <- cnv[cnv$copyNumber > min_cn_ratio * ploidy, , drop = FALSE]
    agg <- if (nrow(amp) > 0) do.call(rbind, lapply(split(amp, amp$chr), function(d)
      data.frame(chr = d$chr[1], start = min(d$start), end = max(d$end),
                 amp_bp = sum(d$end - d$start), stringsAsFactors = FALSE))) else
      data.frame(chr = character(), start = numeric(), end = numeric(), amp_bp = numeric())
    agg <- agg[agg$amp_bp >= min_amp_width, , drop = FALSE]
    if (nrow(agg) > 0) {
      w <- agg$end - agg$start
      agg$start <- pmax(0, round(agg$start - margin * w))
      agg$end   <- pmin(chr_len[agg$chr], round(agg$end + margin * w))
    }
    if (is.null(chroms)) {
      if (nrow(agg) == 0) stop("No amplified loci found for ", this_sample,
                               " (copyNumber > ", min_cn_ratio, " x ploidy, >= ",
                               min_amp_width, " bp). Specify `chromosome` or `loci`.")
      res <- agg[, c("chr", "start", "end")]
      res <- res[order(suppressWarnings(as.integer(gsub("chr", "", res$chr))), res$chr), ]
    } else {
      res <- do.call(rbind, lapply(chroms, function(cc) {
        if (cc %in% agg$chr) agg[agg$chr == cc, c("chr", "start", "end")]
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

  ## ---- Concatenated layout ------------------------------------------------
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

  ## ---- CN segments in plot coords -----------------------------------------
  cn_plot <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    d <- cnv[cnv$chr == loci$chr[i] & cnv$end >= loci$start[i] & cnv$start <= loci$end[i], , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$start <- pmax(d$start, loci$start[i]); d$end <- pmin(d$end, loci$end[i])
    d$gx_start <- loci$offset[i] + (d$start - loci$start[i])
    d$gx_end   <- loci$offset[i] + (d$end   - loci$start[i])
    d
  }))
  if (is.null(cn_plot) || nrow(cn_plot) == 0) stop("No copy-number data in the selected loci.")

  ploidy <- cnv$ploidy[1]
  max.cn <- max(cn_plot$majorAlleleCopyNumber)

  ## ---- Axis scaling (shared CN + read-support, 2-sig-fig tops) ------------
  nice_step <- function(span, n = 5) {
    if (!is.finite(span) || span <= 0) return(1)
    raw <- span / n; mag <- 10^floor(log10(raw)); frac <- raw / mag
    (if (frac < 1.5) 1 else if (frac < 3) 2 else if (frac < 7) 5 else 10) * mag
  }
  nice_ceiling <- function(x, digits = 2) {
    if (!is.finite(x) || x <= 0) return(1)
    mag <- 10^(floor(log10(x)) - (digits - 1)); ceiling(x / mag) * mag
  }
  cn_axis_max <- if (!is.null(cn_max)) cn_max else nice_ceiling(max.cn)

  ## ---- SVs ----------------------------------------------------------------
  sv <- sv_data[sv_data$sample == this_sample, , drop = FALSE]
  sv$chr1 <- paste0("chr", sv$chrom1); sv$chr2 <- paste0("chr", sv$chrom2)
  sv$pos1 <- as.integer(sv$start1);    sv$pos2 <- as.integer(sv$start2)
  sv$strands <- paste0(sv$strand1, sv$strand2)
  sv$strands <- ifelse(sv$svclass == "TRA", "TRA", sv$strands)
  sv$l1 <- locus_idx(sv$chr1, sv$pos1)
  sv$l2 <- locus_idx(sv$chr2, sv$pos2)
  sv <- sv[!(is.na(sv$l1) & is.na(sv$l2)), , drop = FALSE]

  DEL_colour <- "#bde0fe"; DUP_colour <- "#c1447e"; h2hINV_colour <- "#a5a6ae"
  t2tINV_colour <- "#384351"; TRA_colour <- "#fac881"
  sv_col <- function(s) ifelse(s %in% c("+-", "DEL"), DEL_colour,
                        ifelse(s %in% c("++", "h2hINV"), h2hINV_colour,
                        ifelse(s %in% c("--", "t2tINV"), t2tINV_colour,
                        ifelse(s %in% c("-+", "DUP"), DUP_colour, TRA_colour))))
  sv$colour <- sv_col(sv$strands)

  max_vf <- if (nrow(sv) > 0) max(sv$VF, na.rm = TRUE) else 0
  rs_axis_max <- if (max_vf > 0) nice_ceiling(max_vf) else 0
  coeff <- if (max_vf > 0) rs_axis_max / cn_axis_max else 1

  ## ---- Colours / sizes ----------------------------------------------------
  color_minor_cn <- "#3a9387"; color_major_cn <- "#d92a05"
  color_homdel <- "#0077b6"; color_loh <- "#bde0fe"
  cn_size <- 0.9; size_sv_line <- 0.2; size_interchr_line <- 0.3
  size_gene_label <- 3; size_text <- 8
  curv_intra <- 0.18; curv_inter <- -0.18

  ## ---- Karyotype / gap geometry -------------------------------------------
  if (max.cn > 250) { minorAllele_offset <- 5; upper_limit_karyotype <- -max.cn * 0.14
  } else if (max.cn > 110) { minorAllele_offset <- 3; upper_limit_karyotype <- -max.cn * 0.16
  } else if (max.cn > 40) { minorAllele_offset <- 1; upper_limit_karyotype <- -max.cn * 0.18
  } else if (max.cn < 6) { minorAllele_offset <- 0.05; upper_limit_karyotype <- -max.cn * 0.25
  } else { minorAllele_offset <- 0.1; upper_limit_karyotype <- -max.cn * 0.10 }
  yend_outside_range <- max.cn * 0.03
  lower_limit_karyotype <- upper_limit_karyotype - (cn_axis_max * karyotype_rel_size)
  loh_bar_y <- upper_limit_karyotype * loh_position_ratio

  ## ---- Build plot ---------------------------------------------------------
  p <- ggplot()

  ideo <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    d <- karyotype[karyotype$chr == loci$chr[i] & karyotype$end >= loci$start[i] & karyotype$start <= loci$end[i], , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$start <- pmax(d$start, loci$start[i]); d$end <- pmin(d$end, loci$end[i])
    d$gx_start <- loci$offset[i] + (d$start - loci$start[i])
    d$gx_end   <- loci$offset[i] + (d$end   - loci$start[i])
    d
  }))
  if (!is.null(ideo) && nrow(ideo) > 0)
    p <- p + geom_rect(data = ideo, aes(xmin = gx_start, xmax = gx_end,
             ymin = lower_limit_karyotype, ymax = upper_limit_karyotype),
             fill = ideo$color, colour = "black", linewidth = 0.2)

  if (highlight_amp) {
    d <- cn_plot[cn_plot$copyNumber > 3 * cn_plot$ploidy, , drop = FALSE]
    if (nrow(d) > 0) p <- p + geom_rect(data = d, aes(xmin = gx_start, xmax = gx_end,
      ymin = 0, ymax = max.cn * ymax_highlight_ratio), fill = "#d92a05", alpha = 0.1)
  }
  if (highlight_hom_del) {
    d <- cn_plot[cn_plot$copyNumber < 0.5, , drop = FALSE]
    if (nrow(d) > 0) p <- p + geom_rect(data = d, aes(xmin = gx_start, xmax = gx_end,
      ymin = 0, ymax = max.cn * ymax_highlight_ratio), fill = color_homdel, alpha = 0.05)
  }

  p <- p +
    geom_segment(data = cn_plot, aes(x = gx_start, xend = gx_end,
                 y = minorAlleleCopyNumber - minorAllele_offset,
                 yend = minorAlleleCopyNumber - minorAllele_offset),
                 colour = color_minor_cn, linewidth = cn_size) +
    geom_segment(data = cn_plot, aes(x = gx_start, xend = gx_end,
                 y = majorAlleleCopyNumber, yend = majorAlleleCopyNumber),
                 colour = color_major_cn, linewidth = cn_size)
  loh <- cn_plot[cn_plot$minorAlleleCopyNumber < 0.5, , drop = FALSE]
  if (nrow(loh) > 0) p <- p + geom_segment(data = loh,
    aes(x = gx_start, xend = gx_end, y = loh_bar_y, yend = loh_bar_y),
    colour = color_loh, linewidth = cn_size)
  hd <- cn_plot[cn_plot$copyNumber < 0.5, , drop = FALSE]
  if (nrow(hd) > 0) p <- p + geom_segment(data = hd,
    aes(x = gx_start, xend = gx_end, y = loh_bar_y, yend = loh_bar_y),
    colour = color_homdel, linewidth = cn_size)

  ## SV arcs + breakpoint lines (batched by curvature/linewidth for speed):
  arc_rows <- list(); seg_rows <- list()
  for (i in seq_len(nrow(sv))) {
    yv <- sv$VF[i] / coeff; col <- sv$colour[i]
    in1 <- !is.na(sv$l1[i]); in2 <- !is.na(sv$l2[i])
    gx1 <- if (in1) loci$offset[sv$l1[i]] + (sv$pos1[i] - loci$start[sv$l1[i]]) else NA
    gx2 <- if (in2) loci$offset[sv$l2[i]] + (sv$pos2[i] - loci$start[sv$l2[i]]) else NA
    if (in1 && in2) {
      same <- sv$l1[i] == sv$l2[i]
      cv <- if (same) {
        base <- if (abs(gx2 - gx1) <= total_w * 0.15) curv_intra * 1.4 else curv_intra
        if (sv$strands[i] %in% c("DEL", "h2hINV", "+-", "--")) base else -base
      } else curv_inter
      arc_rows[[length(arc_rows) + 1L]] <- data.frame(x = gx1, xend = gx2, y = yv, yend = yv,
        curvature = cv, colour = col, lwd = if (same) size_sv_line else size_interchr_line)
      seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx1, xend = gx1, y = 0, yend = yv, colour = col)
      seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx2, xend = gx2, y = 0, yend = yv, colour = col)
    } else {
      gx <- if (in1) gx1 else gx2
      seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx, xend = gx, y = 0, yend = yv, colour = col)
      seg_rows[[length(seg_rows) + 1L]] <- data.frame(x = gx, xend = gx - total_w * 0.006, y = 0, yend = -yend_outside_range, colour = col)
    }
  }
  if (length(arc_rows) > 0) {
    arc_df <- as.data.frame(data.table::rbindlist(arc_rows))
    arc_df$grp <- paste(arc_df$curvature, arc_df$lwd, sep = "_")
    for (g in unique(arc_df$grp)) {
      d <- arc_df[arc_df$grp == g, , drop = FALSE]
      p <- p + geom_curve(data = d, aes(x = x, xend = xend, y = y, yend = yend, colour = colour),
                          curvature = d$curvature[1], linewidth = d$lwd[1])
    }
  }
  if (length(seg_rows) > 0) {
    seg_df <- as.data.frame(data.table::rbindlist(seg_rows))
    p <- p + geom_segment(data = seg_df, aes(x = x, xend = xend, y = y, yend = yend, colour = colour),
                          linewidth = size_sv_line)
  }
  p <- p + scale_colour_identity()

  amp_seg <- cn_plot[cn_plot$copyNumber > 3 * cn_plot$ploidy, , drop = FALSE]
  if (nrow(amp_seg) > 0) p <- p + geom_segment(data = amp_seg,
    aes(x = gx_start, xend = gx_end, y = majorAlleleCopyNumber, yend = majorAlleleCopyNumber),
    colour = color_major_cn, linewidth = cn_size - 0.3)

  ## Gene labels (in-window), de-collided with ggrepel:
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
    if (is.null(gene_label_angle)) {
      gp <- sort(gl$gx); gene_label_angle <- if (length(gp) > 1 && min(diff(gp)) < 0.05 * total_w) 45 else 0
    }
    p <- p + geom_point(data = gl, aes(x = gx, y = y * 0.94), shape = 16, size = 1, colour = "black")
    if (repel_labels && requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = gl, aes(x = gx, y = y, label = gene),
        size = size_gene_label, fontface = "italic", angle = gene_label_angle, hjust = 0,
        direction = "both", nudge_y = max.cn * 0.12, ylim = c(max.cn * 1.02, NA),
        segment.size = 0.2, segment.colour = "grey60", min.segment.length = 0,
        box.padding = 0.25, max.overlaps = Inf, seed = 1L)
    } else {
      p <- p + geom_text(data = gl, aes(x = gx, y = y, label = gene),
                         size = size_gene_label, fontface = "italic", angle = gene_label_angle)
    }
  }

  ## ---- Bottom labelling (Mb ticks + Mb labels + chromosome names) ---------
  unit_y      <- cn_axis_max * 0.07
  tick_y0     <- lower_limit_karyotype
  tick_y1     <- lower_limit_karyotype - 0.45 * unit_y
  mb_label_y  <- lower_limit_karyotype - 1.25 * unit_y
  chr_label_y <- lower_limit_karyotype - 2.7 * unit_y

  tick_df <- data.frame(); mb_df <- data.frame()
  for (i in seq_len(nrow(loci))) {
    ticks <- pretty(c(loci$start[i], loci$end[i]), n = 3)
    ticks <- ticks[ticks >= loci$start[i] & ticks <= loci$end[i]]
    gxt <- loci$offset[i] + (ticks - loci$start[i])
    tick_df <- rbind(tick_df, data.frame(x = gxt))
    mb_df <- rbind(mb_df, data.frame(x = gxt, lab = formatC(ticks / 1e6, format = "f", digits = 1)))
  }
  chr_lab <- data.frame(gx = loci$gx_mid, chr = gsub("chr", "Chr ", loci$chr))
  cn_breaks <- unique(c(0, cn_axis_max / 2, cn_axis_max))

  p <- p +
    geom_segment(data = data.frame(x = loci$gx_start[1]),
                 aes(x = x, xend = x, y = 0, yend = cn_axis_max), colour = "black", linewidth = 0.4) +
    geom_segment(data = data.frame(x = loci$gx_end[nrow(loci)]),
                 aes(x = x, xend = x, y = 0, yend = cn_axis_max), colour = "black", linewidth = 0.4) +
    geom_segment(data = tick_df, aes(x = x, xend = x, y = tick_y0, yend = tick_y1),
                 colour = "black", linewidth = 0.3) +
    geom_text(data = mb_df, aes(x = x, y = mb_label_y, label = lab), size = size_text / 2.6, vjust = 1) +
    geom_text(data = chr_lab, aes(x = gx, y = chr_label_y, label = chr),
              size = size_text / 2.0, vjust = 1, fontface = "bold") +
    ggtitle(paste0(this_sample, " (", wgd_status, ")")) +
    labs(x = NULL, y = "Allele specific\ncopy number") +
    coord_cartesian(clip = "off", expand = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    theme(
      text = element_text(size = size_text, colour = "black"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_text(size = size_text, colour = "black"),
      axis.title.y = element_text(size = size_text + 2, colour = "black"),
      axis.ticks.y.right = element_line(colour = "black", linewidth = 0.4),
      panel.background = element_blank(), plot.background = element_blank(), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = size_text + 3),
      plot.margin = unit(c(.3, .5, 1.2, .2), "cm")
    )

  if (max_vf > 0) {
    rs_breaks <- unique(c(0, rs_axis_max / 2, rs_axis_max))
    p <- p + scale_y_continuous(breaks = cn_breaks,
      sec.axis = sec_axis(trans = ~ . * coeff, breaks = rs_breaks, name = "Read support"))
  } else {
    p <- p + scale_y_continuous(breaks = cn_breaks)
  }

  ## ---- Size + save --------------------------------------------------------
  plot_width  <- if (!is.null(plot_width_custom)) plot_width_custom else max(5, 2.2 * nrow(loci) + 1)
  plot_height <- if (!is.null(plot_height_custom)) plot_height_custom else 3.0

  outfile <- NULL
  if (!is.null(outdir) && save) {
    format <- match.arg(format, several.ok = TRUE)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    region_tag <- if (nrow(loci) == 1)
      sprintf("_%.0f-%.0f", round(loci$start[1]), round(loci$end[1])) else ""
    stem <- file.path(outdir, paste0(this_sample, "_", paste(loci$chr, collapse = "_"),
                                     region_tag, "_linear_plot"))
    written <- character(0)
    if ("pdf" %in% format) {
      pp <- paste0(stem, ".pdf")
      grDevices::pdf(pp, width = plot_width, height = plot_height, useDingbats = FALSE)
      print(p); grDevices::dev.off(); written <- c(written, pp)
    }
    if ("png" %in% format) {
      pp <- paste0(stem, ".png")
      ggplot2::ggsave(pp, p, width = plot_width, height = plot_height, dpi = dpi, bg = "white")
      written <- c(written, pp)
    }
    outfile <- written
    if (verbose) message("Wrote ", paste(written, collapse = ", "))
  }
  attr(p, "path") <- outfile
  if (!is.null(outfile)) invisible(p) else p
}
