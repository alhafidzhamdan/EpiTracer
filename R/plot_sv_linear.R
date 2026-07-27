#' Draw a linear copy-number and structural-variant "recon" plot
#'
#' Produces a linear, allele-specific copy-number track for one or more
#' chromosomes (or sub-chromosomal windows) of a single sample, overlaid with
#' arcs for intra- and inter-chromosomal structural variants (SVs), a karyotype
#' ideogram, LOH / homozygous-deletion bars, and gene labels. The assembled
#' `ggplot` object is returned, and optionally written to a PDF in `outdir`.
#'
#' Reference annotation that was previously hard-coded is now supplied via the
#' `karyotype`, `gene_coord`, and `cds_gr` arguments so the function is
#' self-contained.
#'
#' @param sample Character scalar; the sample identifier. Rows with
#'   `sample == sample` are selected from `wgd_data`, `cnv_data`, and `sv_data`.
#' @param chromosome Character vector of chromosomes to display, e.g.
#'   `c("chr7", "chr12")`.
#' @param chromosome_range Optional two-column matrix/data.frame of
#'   `start`,`end` window limits, one row per entry in `chromosome`. If `NULL`
#'   (default), the full chromosome (or the SV-spanning window) is shown.
#' @param genes_to_highlight Optional character vector of gene symbols to label.
#'   If `NULL`, a default oncogene panel is used.
#' @param karyotype A data.frame of ideogram bands (as from the UCSC
#'   `cytoBand`/`chr_info` table) with columns including chromosome, `start`,
#'   `end`, and `gieStain`; the first column is treated as the chromosome name.
#'   A file path to an `.rds` holding such a data.frame is also accepted.
#' @param gene_coord A data.frame of gene coordinates with columns
#'   `chr`,`start`,`end`,`strand`,`gene`. A path to a headerless tab-separated
#'   BED-like file with those five columns is also accepted.
#' @param wgd_data A data.frame with a sample-identifier column (see
#'   `wgd_sample_col`) and a `Polyploidy` column (`"No"` = diploid, otherwise
#'   WGD) used to label the plot title.
#' @param cnv_data A data.frame (or object coercible via `as.data.frame`) of
#'   copy-number segments with columns `sample`, `seqnames`, `start`, `end`,
#'   `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`.
#' @param sv_data A data.frame of SVs with columns `chrom1`, `start1`, `chrom2`,
#'   `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`.
#' @param cds_gr Optional [GenomicRanges::GRanges] of CDS/exon ranges (metadata
#'   column `gene_name`) used only when `displayExon = TRUE` to draw exon models.
#' @param wgd_sample_col Optional name of the sample-identifier column in
#'   `wgd_data`. If `NULL` (default) the function uses `sample`, falling back to
#'   `WGS_ID` if present.
#' @param sec_axis_adj,offset_gene,ymax_highlight_ratio Numeric layout tuning
#'   parameters (see Details of the original implementation).
#' @param displayExon Logical; if `TRUE`, draw exon models (requires `cds_gr`).
#' @param scale_ticks Numeric x-axis tick spacing in bp. If `NULL` (default) a
#'   sensible spacing is derived from the width of the widest panel.
#' @param yend_left,yend_right,yend_right_line,loh_position_ratio Numeric layout
#'   tuning parameters for axes and LOH bars. `yend_left`/`yend_right` may be
#'   `NULL` (default) to auto-scale the copy-number axis to the data.
#' @param repel_labels Logical; if `TRUE` (default) use \pkg{ggrepel} to keep
#'   gene labels from overlapping (falls back to plain labels if ggrepel is not
#'   installed).
#' @param plot_height_custom,plot_width_custom Optional numeric overrides for the
#'   auto-selected PDF dimensions (inches).
#' @param highlight_amp,highlight_hom_del Logical; shade amplified /
#'   homozygously deleted segments.
#' @param karyotype_rel_size Numeric; ideogram height relative to the CN axis.
#' @param outdir Optional directory in which to write the PDF. If `NULL`
#'   (default) no file is written and only the plot object is returned.
#' @param save Logical; if `TRUE` (default) and `outdir` is supplied, write the
#'   PDF. Set `FALSE` to build the plot without writing a file.
#' @param verbose Logical; if `TRUE`, print progress/diagnostic messages.
#'
#' @return The assembled [ggplot2::ggplot] object (invisibly when a PDF is
#'   written). The path of any written PDF is attached as attribute
#'   `"path"`.
#'
#' @examples
#' \dontrun{
#' p <- plot_sv_linear(
#'   sample     = "DO12742T1",
#'   chromosome = c("chr7", "chr12"),
#'   karyotype  = "chr_info_hg38.rds",
#'   gene_coord = "gene.coord_strand_name.bed",
#'   wgd_data   = wgd_df,
#'   cnv_data   = cnv_df,
#'   sv_data    = sv_df,
#'   outdir     = "plots"
#' )
#' p + ggplot2::labs(subtitle = "EGFR amplicon")   # compose further
#' }
#'
#' @seealso [call_episomal_ecdna()]
#' @export
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom scales number_format
#' @importFrom data.table rbindlist
#' @importFrom dplyr filter select rename mutate distinct
#' @importFrom grDevices pdf dev.off
#' @importFrom utils read.table
plot_sv_linear <- function(sample,
                           chromosome,
                           chromosome_range = NULL,
                           genes_to_highlight = NULL,
                           karyotype,
                           gene_coord,
                           wgd_data,
                           cnv_data,
                           sv_data,
                           cds_gr = NULL,
                           wgd_sample_col = NULL,
                           sec_axis_adj = 10,
                           offset_gene = 1.2,
                           ymax_highlight_ratio = 1.08,
                           displayExon = FALSE,
                           scale_ticks = NULL,
                           yend_left = NULL,
                           yend_right = 100,
                           yend_right_line = 2,
                           loh_position_ratio = 0.3,
                           repel_labels = TRUE,
                           plot_height_custom = NULL,
                           plot_width_custom = NULL,
                           highlight_amp = TRUE,
                           highlight_hom_del = TRUE,
                           karyotype_rel_size = 0.05,
                           outdir = NULL,
                           save = TRUE,
                           verbose = FALSE) {

  #### Base functions ####
  chromosome_labeller <- function(chr, value) {
    chrm_name <- gsub("chr", "", chr)
    return(paste0("Chr ", chrm_name))
  }
  ## Tick breaks that span the actual plotted window (respecting both ends),
  ## so zoomed panels still get labelled ticks.
  scaler <- function(step) {
    function(y) {
      lo <- floor(min(y) / step) * step
      hi <- ceiling(max(y) / step) * step
      seq(lo, hi, by = step)
    }
  }
  ## Choose a "nice" tick spacing (1/2/5 x 10^k) targeting ~n intervals.
  nice_step <- function(span, n = 5) {
    if (!is.finite(span) || span <= 0) return(10e6)
    raw <- span / n
    mag <- 10^floor(log10(raw))
    frac <- raw / mag
    mult <- if (frac < 1.5) 1 else if (frac < 3) 2 else if (frac < 7) 5 else 10
    mult * mag
  }
  ## Round a value up to a "nice" multiple, for a tidy copy-number axis top.
  nice_ceiling <- function(x) {
    if (!is.finite(x) || x <= 0) return(1)
    s <- nice_step(x, 2)
    ceiling(x / s) * s
  }

  #### Build karyotype ####
  if (is.character(karyotype)) karyotype <- readRDS(karyotype)
  karyotype_data <- as.data.frame(karyotype)
  names(karyotype_data)[1] <- "chr"
  karyotype_data$color[karyotype_data$gieStain == "gneg"]    <- "white"
  karyotype_data$color[karyotype_data$gieStain == "gpos25"]  <- "grey75"
  karyotype_data$color[karyotype_data$gieStain == "gpos50"]  <- "grey50"
  karyotype_data$color[karyotype_data$gieStain == "gpos75"]  <- "grey25"
  karyotype_data$color[karyotype_data$gieStain == "gpos100"] <- "grey0"
  karyotype_data$color[karyotype_data$gieStain == "acen"]    <- "red"

  #### Get gene coordinates ####
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

  #### Get sample and chromosomes info ####
  this_sample <- sample
  chrs <- chromosome
  if (is.null(chromosome_range)) {
    chr_selection <- data.frame(chr = chrs, start = rep(0, length(chrs)),
                                end = rep(250000000, length(chrs)))
  } else {
    ## NB: column must be `chr` (not `chrs`) to match every downstream
    ## `chr_selection$chr` reference; the original script mis-named it, which
    ## silently broke the `chromosome_range` (zoom) code path.
    chr_selection <- data.frame(chr = chrs, start = chromosome_range[, 1],
                                end = chromosome_range[, 2])
  }

  if (verbose) {
    message("sample: ", this_sample)
    message("chromosomes: ", paste(chrs, collapse = ", "))
  }

  #### Get wgd status ####
  ## Accept a WGD table keyed by `sample` (default) or `WGS_ID`, or a
  ## user-specified column name.
  wgd_data <- as.data.frame(wgd_data)
  wgd_col <- wgd_sample_col
  if (is.null(wgd_col)) {
    wgd_col <- if ("sample" %in% names(wgd_data)) "sample"
               else if ("WGS_ID" %in% names(wgd_data)) "WGS_ID"
               else stop("wgd_data must contain a 'sample' or 'WGS_ID' column, ",
                         "or set `wgd_sample_col`.")
  }
  wgd_row <- wgd_data[wgd_data[[wgd_col]] == this_sample, , drop = FALSE]
  if (nrow(wgd_row) >= 1 && wgd_row$Polyploidy[1] == "No") {
    wgd_status_this_sample <- "Diploid"
  } else {
    wgd_status_this_sample <- "WGD"
  }
  if (verbose) message("WGD status: ", wgd_status_this_sample)

  #### Define parameters ####
  cn_size <- 0.9
  if (length(chrs) > 6) {
    plot_height <- 3.5; plot_width <- 9
  } else if (length(chrs) > 3) {
    plot_height <- 3; plot_width <- 7.5
  } else if (length(chrs) > 2) {
    plot_height <- 2.7; plot_width <- 6
  } else {
    plot_height <- 2.5; plot_width <- 4.5
  }
  if (!is.null(plot_width_custom))  plot_width  <- plot_width_custom
  if (!is.null(plot_height_custom)) plot_height <- plot_height_custom

  window <- 10000000
  xscale <- 1e6
  npc_now <- .00625 * 1.5   ## gap between chromosome panels; keep < 3
  curvature_intrachr_SVs <- -0.15
  curvature_interchr_SVs <- -0.15
  label_interchr_SV <- FALSE

  ## Colours and text sizes:
  TRA_colour <- "#fac881"
  DEL_colour <- "#bde0fe"
  DUP_colour <- "#c1447e"
  h2hINV_colour <- "#a5a6ae"
  t2tINV_colour <- "#384351"
  color_minor_cn <- "#3a9387"
  color_major_cn <- "#d92a05"
  color_homdel <- "#0077b6"
  color_loh <- "#bde0fe"
  size_text <- 8
  size_title <- 10
  size_chr_labels <- 10
  size_sv_line <- 0.2
  size_interchr_line <- 0.3
  size_gene_label <- 3

  #### Load CNV, SV and karyotype data: ####
  cnv <- cnv_data %>% dplyr::filter(sample == this_sample) %>%
    dplyr::rename(chr = seqnames) %>% as.data.frame()
  if (verbose) message("CN segments for sample: ", nrow(cnv))

  sv <- sv_data %>%
    dplyr::mutate(chrom1 = paste0("chr", chrom1), chrom2 = paste0("chr", chrom2)) %>%
    dplyr::mutate(strands = paste0(strand1, strand2)) %>%
    dplyr::mutate(strands = ifelse(svclass == "TRA", svclass, strands)) %>%
    dplyr::select(chr1 = chrom1, pos1 = start1, chr2 = chrom2, pos2 = start2,
                  strands, VF, JCN, sample) %>%
    dplyr::filter(sample == this_sample) %>%
    as.data.frame()
  if (verbose) message("Translocations (TRA) for sample: ",
                       nrow(sv[sv$strands == "TRA", ]))

  ## Validate SVs with no width info (assume intrachromosomal):
  idx <- which(sv$pos1 == sv$pos2)
  if (length(idx) > 0) sv$pos2[idx] <- sv$pos2[idx] + 1

  ## Karyotype info:
  karyotype_data_now <- karyotype_data[karyotype_data$chr %in% chr_selection$chr, ]
  karyotype_data_now$y <- rep(1, nrow(karyotype_data_now))
  if (nrow(karyotype_data_now) < 7) {
    karyotype_data_annot <- karyotype_data_now
  } else {
    karyotype_data_annot <- karyotype_data_now[seq(3, (nrow(karyotype_data_now) - 3), 3), ]
  }
  karyotype_data_now$chr <- factor(karyotype_data_now$chr, levels = unique(chr_selection$chr))

  #### Process data: ####
  sv$pos1 <- as.integer(sv$pos1)
  sv$pos2 <- as.integer(sv$pos2)
  sv$chr1 <- as.character(sv$chr1)
  sv$chr2 <- as.character(sv$chr2)

  ## Arc colours and curvatures:
  if (nrow(sv) >= 1) {
    sv$colour <- apply(sv, 1, function(x) {
      if (x["strands"] == "+-" | x["strands"] == "DEL") {
        c <- DEL_colour
      } else if (x["strands"] == "++" | x["strands"] == "h2hINV") {
        c <- h2hINV_colour
      } else if (x["strands"] == "--" | x["strands"] == "t2tINV") {
        c <- t2tINV_colour
      } else if (x["strands"] == "-+" | x["strands"] == "DUP") {
        c <- DUP_colour
      } else if (x["strands"] == "TRA") {
        c <- TRA_colour
      }
      return(c)
    })
    sv$curve <- apply(sv, 1, function(x) {
      if (as.character(x["chr1"]) == as.character(x["chr2"])) {
        if (abs(as.integer(x["pos2"]) - as.integer(x["pos1"])) <= 10000000) {
          c <- 2 * curvature_intrachr_SVs
        } else {
          c <- curvature_intrachr_SVs
        }
      } else {
        c <- curvature_interchr_SVs
      }
      c
    })
  }

  ## Intrachromosomal SVs (of the chromosomes of interest):
  intraSV <- lapply(seq_along(chr_selection$chr), function(i) {
    sv[sv$chr1 == chr_selection$chr[i] & sv$chr2 == chr_selection$chr[i], ]
  }) %>% data.table::rbindlist()

  ## Flag interchromosomal SVs touching the chromosomes of interest:
  inter <- sv[sv$chr1 != sv$chr2, ]
  if ((sum(inter$chr1 %in% chr_selection$chr) + sum(inter$chr2 %in% chr_selection$chr)) > 0) {
    interFlag <- TRUE
  } else {
    interFlag <- FALSE
  }

  ## X-axis limits for each chromosome of interest:
  for (i in 1:nrow(chr_selection)) {
    chr.now <- chr_selection$chr[i]
    start.now <- chr_selection$start[i]
    end.now <- chr_selection$end[i]
    if (is.na(start.now) | is.na(end.now)) {
      pos1 <- sv$pos1[sv$chr1 == chr.now]
      if (length(pos1) > 0) pos1 <- min(pos1)
      pos2 <- sv$pos2[sv$chr2 == chr.now]
      if (length(pos2) > 0) pos2 <- min(pos2)
      start.now <- c(pos1, pos2)
      idx <- unlist(sapply(start.now, is.numeric))
      if (length(idx) > 0) {
        start.now <- min(start.now[idx])
      } else {
        start.now <- 0
      }
      chr_selection$start[i] <- start.now - window
      if (chr_selection$start[i] < 0) chr_selection$start[i] <- 0

      end1 <- sv$pos1[sv$chr1 == chr.now]
      if (length(end1) > 0) end1 <- max(end1)
      end2 <- sv$pos2[sv$chr2 == chr.now]
      if (length(end2) > 0) end2 <- max(end2)
      end.now <- c(end1, end2)
      idx <- unlist(sapply(end.now, is.numeric))
      if (length(idx) > 0) {
        end.now <- max(end.now[idx])
      } else {
        end.now <- 250000000
      }
      chr_selection$end[i] <- end.now + window
    }
    if (chr_selection$end[i] > max(karyotype_data_now$end[karyotype_data_now$chr %in% chr_selection$chr])) {
      chr_selection$end[i] <- max(karyotype_data_now$end[karyotype_data_now$chr %in% chr_selection$chr])
    }
  }

  ## Auto x-axis tick spacing from the widest panel, unless the user set it.
  if (is.null(scale_ticks)) {
    widest <- max(chr_selection$end - chr_selection$start, na.rm = TRUE)
    scale_ticks <- nice_step(widest)
    if (verbose) message("Auto x-axis tick spacing: ", scale_ticks, " bp")
  }

  ## TRAs + CNV/karyotype setup:
  if (interFlag) {
    idx1 <- which(sv$chr1 != sv$chr2 & sv$chr1 %in% chr_selection$chr)
    idx2 <- which(sv$chr1 != sv$chr2 & sv$chr2 %in% chr_selection$chr)
    if (length(idx1) > 0 | length(idx2) > 0) {
      interSV <- sv[unique(c(idx1, idx2)), ]
      interSV$pos1 <- as.integer(interSV$pos1)
      interSV$pos2 <- as.integer(interSV$pos2)
    }

    idx1 <- which(sv$chr1 != sv$chr2 & sv$chr1 %in% chr_selection$chr & !(sv$chr2 %in% chr_selection$chr))
    idx2 <- which(sv$chr1 != sv$chr2 & sv$chr2 %in% chr_selection$chr & !(sv$chr1 %in% chr_selection$chr))
    if (length(idx1) > 0 | length(idx2) > 0) {
      interSV_other_chrs <- sv[unique(c(idx1, idx2)), ]
      interSV_other_chrs$pos1 <- as.integer(interSV_other_chrs$pos1)
      interSV_other_chrs$pos2 <- as.integer(interSV_other_chrs$pos2)
    }
  }

  ## Process CNV data (same regardless of interFlag):
  cnv.plot <- lapply(1:nrow(chr_selection), function(i) {
    main.cnv <- subset(cnv, chr == chr_selection$chr[i] &
                         ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                            (start <= chr_selection$start[i] & end >= chr_selection$start[i]) |
                            (start >= chr_selection$start[i] & start <= chr_selection$end[i])))
    main.cnv[main.cnv$start < chr_selection$start[i], "start"] <- chr_selection$start[i]
    main.cnv[main.cnv$end > chr_selection$end[i], "end"] <- chr_selection$end[i]
    main.cnv
  }) %>% data.table::rbindlist() %>% as.data.frame()

  ## Prune karyotype data to updated chr_selection windows:
  karyotype_data_now <- lapply(1:nrow(chr_selection), function(i) {
    karyo_subset <- subset(karyotype_data_now, chr == chr_selection$chr[i] &
                             ((start >= chr_selection$start[i] & end <= chr_selection$end[i]) |
                                (start <= chr_selection$start[i] & end >= chr_selection$end[i]) |
                                (start >= chr_selection$start[i] & start <= chr_selection$end[i]) |
                                (end >= chr_selection$start[i] & end <= chr_selection$end[i])))
    karyo_subset[karyo_subset$start < chr_selection$start[i], "start"] <- chr_selection$start[i]
    karyo_subset[karyo_subset$end > chr_selection$end[i], "end"] <- chr_selection$end[i]
    karyo_subset
  }) %>% data.table::rbindlist() %>% as.data.frame()

  ## High CN handling:
  max.cn <- max(cnv.plot$majorAlleleCopyNumber)
  max_y <- max.cn
  max_y_rectangle <- max_y

  #### Plot: ####
  p <- ggplot()

  ## Karyotype scaling by max major CN:
  max_major <- max(cnv.plot[cnv.plot$chr %in% chr_selection$chr, ]$majorAlleleCopyNumber)
  if (max_major > 250) {
    minorAllele_offset <- 5;    upper_limit_karyotype <- -100; yend_outside_range <- 5
  } else if (max_major > 110) {
    minorAllele_offset <- 3;    upper_limit_karyotype <- -50;  yend_outside_range <- 5
  } else if (max_major > 40) {
    minorAllele_offset <- 1;    upper_limit_karyotype <- -22;  yend_outside_range <- 5
  } else if (max_major < 6) {
    minorAllele_offset <- 0.05; upper_limit_karyotype <- -1.4; yend_outside_range <- 0.3
  } else {
    minorAllele_offset <- 0.1;  upper_limit_karyotype <- -5;   karyotype_rel_size <- 0.05
    yend_outside_range <- 0.5
  }

  karyotype_space <- max_y_rectangle * karyotype_rel_size
  lower_limit_karyotype <- upper_limit_karyotype - (max_y_rectangle * karyotype_rel_size)
  p <- p + geom_rect(data = karyotype_data_now,
                     mapping = aes(xmin = start, xmax = end,
                                   ymin = lower_limit_karyotype,
                                   ymax = upper_limit_karyotype,
                                   group = chr),
                     fill = karyotype_data_now$color, color = "black", linewidth = 0.2)

  ## y-axis (secondary, read support) scaling:
  coeff <- ceiling(max(sv$VF / 1000)) * 1000 / sec_axis_adj

  if (highlight_amp) {
    p <- p +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0,
                    ymax = max.cn * ymax_highlight_ratio, group = chr),
                data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[] & cnv.plot$copyNumber > 3 * cnv.plot$ploidy, ],
                fill = "#d92a05", alpha = 0.1)
  }
  if (highlight_hom_del) {
    p <- p +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0,
                    ymax = max.cn * ymax_highlight_ratio, group = chr),
                data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[] & cnv.plot$copyNumber < 0.5, ],
                fill = color_homdel, alpha = 0.05)
  }

  p <- p +
    geom_segment(aes(x = start, xend = end,
                     y = minorAlleleCopyNumber - minorAllele_offset,
                     yend = minorAlleleCopyNumber - minorAllele_offset, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr, ], colour = color_minor_cn, linewidth = cn_size) +
    geom_segment(aes(x = start, y = majorAlleleCopyNumber, xend = end,
                     yend = majorAlleleCopyNumber, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr, ], colour = color_major_cn, linewidth = cn_size) +
    geom_segment(aes(x = start, y = -loh_position_ratio * max.cn, xend = end,
                     yend = -loh_position_ratio * max.cn, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr, ] %>% dplyr::filter(minorAlleleCopyNumber < 0.5),
                 colour = color_loh, linewidth = cn_size) +
    geom_segment(aes(x = start, y = -loh_position_ratio * max.cn, xend = end,
                     yend = -loh_position_ratio * max.cn, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr, ] %>% dplyr::filter(copyNumber < 0.5),
                 colour = color_homdel, linewidth = cn_size) +
    facet_grid(. ~ factor(chr, levels = unique(chr_selection$chr)),
               scales = "free_x", space = "free_x", switch = "x",
               labeller = as_labeller(chromosome_labeller)) +
    ggtitle(paste0(this_sample, " (", wgd_status_this_sample, ")")) +
    theme(
      text = element_text(size = size_text, colour = "black"),
      axis.text.x = element_text(size = size_text, colour = "black"),
      axis.title.x = element_text(size = size_text + 3, colour = "black"),
      axis.title.y.left = element_text(size = size_text + 3, colour = "black", hjust = 0.45, vjust = 0.5),
      axis.title.y.right = element_text(size = size_text + 3, colour = "black", hjust = 0.4, vjust = 0.5),
      plot.title = element_text(size = size_title + 1, colour = "black", margin = unit(c(0, 0, 0.5, 0), "cm"), hjust = 0.5),
      axis.text.y.left = element_text(size = size_text + 2, colour = "black"),
      axis.text.y.right = element_text(size = size_text + 2, colour = "black"),
      panel.background = element_blank(),
      strip.background = element_blank(),
      axis.ticks.length.y = unit(0.2, "cm"),
      axis.ticks.length.x = unit(0.15, "cm"),
      strip.placement = "outside",
      strip.clip = "off",
      plot.margin = unit(c(.25, .1, -0.1, .1), "cm"),
      plot.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.spacing.x = unit(npc_now, "npc"),
      strip.text.x = element_text(size = size_chr_labels),
      strip.switch.pad.grid = unit(0, "cm"),
      panel.grid.minor.y = element_line(colour = "transparent"),
      panel.grid.major.y = element_line(colour = "transparent"),
      panel.border = element_blank()
    ) +
    labs(x = "", y = "Allele specific\ncopy number") +
    coord_cartesian(clip = "off", expand = 0) +
    scale_x_continuous(expand = c(0, 0),
                       labels = scales::number_format(scale = 1 / xscale),
                       breaks = scaler(scale_ticks))

  ## Intrachromosomal SVs:
  if (!is.null(intraSV)) {
    if (nrow(intraSV) >= 1) {
      for (i in 1:nrow(intraSV)) {
        max_coord <- max(cnv.plot[cnv.plot$chr == intraSV$chr1[i], "end"])
        min_coord <- min(cnv.plot[cnv.plot$chr == intraSV$chr1[i], "start"])

        intraSV$curve[which(intraSV$strands %in% c("DEL", "h2hINV", "+-", "--"))] <-
          abs(intraSV$curve[which(intraSV$strands %in% c("DEL", "h2hINV", "+-", "--"))])

        if (intraSV$pos1[i] >= min_coord & intraSV$pos2[i] <= max_coord) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos1[i], xend = intraSV$pos2[i],
                       y = intraSV$VF[i] / coeff, yend = intraSV$VF[i] / coeff,
                       curvature = intraSV$curve[i], linewidth = size_sv_line, colour = intraSV$colour[i])
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos1[i], xend = intraSV$pos1[i],
                       y = 0, yend = intraSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_sv_line, colour = intraSV$colour[i]) +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos2[i], xend = intraSV$pos2[i],
                       y = 0, yend = intraSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_sv_line, colour = intraSV$colour[i])
        }

        if (intraSV$pos1[i] >= min_coord & intraSV$pos1[i] < max_coord & intraSV$pos2[i] > max_coord) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos1[i], xend = intraSV$pos1[i],
                       y = 0, yend = intraSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_sv_line, colour = intraSV$colour[i])
          x_range <- (chr_selection$end[which(chr_selection$chr == intraSV$chr1[i])] -
                        chr_selection$start[which(chr_selection$chr == intraSV$chr1[i])]) * 0.01
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos1[i], xend = intraSV$pos1[i] - x_range,
                       y = 0, yend = -yend_outside_range,
                       angle = 45, curvature = 0, linewidth = size_sv_line, colour = intraSV$colour[i])
        }

        if (intraSV$pos1[i] < min_coord & intraSV$pos2[i] <= max_coord & intraSV$pos2[i] > min_coord) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr1[i]),
                       x = intraSV$pos2[i], xend = intraSV$pos2[i],
                       y = 0, yend = intraSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_sv_line, colour = intraSV$colour[i])
          x_range <- (chr_selection$end[which(chr_selection$chr == intraSV$chr1[i])] -
                        chr_selection$start[which(chr_selection$chr == intraSV$chr1[i])]) * 0.01
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = intraSV$chr2[i]),
                       x = intraSV$pos2[i], xend = intraSV$pos2[i] - x_range,
                       y = 0, yend = -yend_outside_range,
                       angle = 45, curvature = 0, linewidth = size_sv_line, colour = intraSV$colour[i])
        }
      }
    }
  }

  ## Interchromosomal SVs among chromosomes of interest:
  if (interFlag) {
    idx <- which(interSV$chr1 %in% chr_selection$chr & interSV$chr2 %in% chr_selection$chr)
    if (length(idx) > 0) {
      interSV <- interSV[idx, ]
      info_chrs <- data.frame(chrs = chr_selection$chr)

      total_chr_size <- 0
      size_in_plot <- c()
      for (ii in chr_selection$chr) {
        size_now <- max(cnv.plot[cnv.plot$chr == ii, "end"]) - min(cnv.plot[cnv.plot$chr == ii, "start"])
        size_in_plot <- c(size_in_plot, size_now)
        total_chr_size <- total_chr_size + size_now
      }
      info_chrs$size_in_plot <- size_in_plot
      gap <- (total_chr_size * npc_now)

      for (i in 1:nrow(interSV)) {
        position_chr_1 <- which(chr_selection$chr == interSV$chr1[i])
        position_chr_2 <- which(chr_selection$chr == interSV$chr2[i])

        if (position_chr_1 > position_chr_2) {
          chr1_now <- interSV$chr1[i]; pos1_now <- interSV$pos1[i]
          chr2_now <- interSV$chr2[i]; pos2_now <- interSV$pos2[i]
          interSV$chr1[i] <- chr2_now; interSV$chr2[i] <- chr1_now
          interSV$pos1[i] <- pos2_now; interSV$pos2[i] <- pos1_now
          chrs_now <- interSV[i, c("chr1", "chr2")]
        } else {
          chrs_now <- interSV[i, c("chr1", "chr2")]
        }

        min_range <- min(cnv.plot[cnv.plot$chr == chrs_now[1, 2], "start"])
        size_leftmost_chr <- max(cnv.plot[cnv.plot$chr == chrs_now[1, 1], "end"])

        offset <- interSV$pos1[i] +
          (size_leftmost_chr - interSV$pos1[i]) +
          (interSV$pos2[i] - min_range) +
          gap * 1.6

        indexes <- seq(1, length(chr_selection$chr))
        indexes <- which(indexes > position_chr_1 & indexes < position_chr_2)
        if (length(indexes) > 0) {
          offset <- offset + (gap * length(indexes)) + sum(info_chrs$size_in_plot[indexes])
        }

        min_pos_chr1 <- min(cnv.plot[cnv.plot$chr == chrs_now[1, 1], "start"])
        max_pos_chr1 <- max(cnv.plot[cnv.plot$chr == chrs_now[1, 1], "end"])
        min_pos_chr2 <- min(cnv.plot[cnv.plot$chr == chrs_now[1, 2], "start"])
        max_pos_chr2 <- max(cnv.plot[cnv.plot$chr == chrs_now[1, 2], "end"])
        in_range_chr1 <- (interSV$pos1[i] > min_pos_chr1 & interSV$pos1[i] < max_pos_chr1)
        in_range_chr2 <- (interSV$pos2[i] > min_pos_chr2 & interSV$pos2[i] < max_pos_chr2)

        if (in_range_chr1 & in_range_chr2) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x = interSV$pos1[i], xend = interSV$pos1[i],
                       y = 0, yend = interSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_interchr_line, colour = interSV$colour[i]) +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x = interSV$pos2[i], xend = interSV$pos2[i],
                       y = 0, yend = interSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_interchr_line, colour = interSV$colour[i])
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = chrs_now[1, 1]),
                       x = interSV$pos1[i], xend = offset,
                       y = interSV$VF[i] / coeff, yend = interSV$VF[i] / coeff,
                       curvature = interSV$curve[i],
                       linewidth = size_interchr_line, colour = interSV$colour[i])
        }

        if (!in_range_chr1) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x = interSV$pos2[i], xend = interSV$pos2[i],
                       y = 0, yend = interSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_interchr_line, colour = interSV$colour[i])
          x_range <- (chr_selection$end[which(chr_selection$chr == interSV$chr1[i])] -
                        chr_selection$start[which(chr_selection$chr == interSV$chr1[i])]) * 0.01
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                       x = interSV$pos2[i], xend = interSV$pos2[i] - x_range,
                       y = 0, yend = -yend_outside_range,
                       angle = 45, curvature = 0, linewidth = size_interchr_line, colour = interSV$colour[i])
          if (label_interchr_SV) {
            p <- p + geom_text(data = data.frame(cov = 1, chr = interSV$chr2[i]),
                               x = interSV$pos2[i] - x_range, y = interSV$VF[i] / coeff,
                               size = size_text / .pt, colour = interSV$colour[i], label = interSV$chr1[i])
          }
        }

        if (!in_range_chr2) {
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x = interSV$pos1[i], xend = interSV$pos1[i],
                       y = 0, yend = interSV$VF[i] / coeff, curvature = 0,
                       linewidth = size_interchr_line, colour = interSV$colour[i])
          x_range <- (chr_selection$end[which(chr_selection$chr == interSV$chr1[i])] -
                        chr_selection$start[which(chr_selection$chr == interSV$chr1[i])]) * 0.01
          p <- p +
            geom_curve(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                       x = interSV$pos1[i], xend = interSV$pos1[i] - x_range,
                       y = 0, yend = -yend_outside_range,
                       angle = 45, curvature = 0, linewidth = size_interchr_line, colour = interSV$colour[i])
          if (label_interchr_SV) {
            p <- p + geom_text(data = data.frame(cov = 1, chr = interSV$chr1[i]),
                               x = interSV$pos1[i] - x_range, y = interSV$VF[i] / coeff,
                               size = size_text / .pt, colour = interSV$colour[i], label = interSV$chr2[i])
          }
        }
      }
    }
  }

  ## Interchromosomal SVs to other chromosomes:
  if (interFlag) {
    if (exists("interSV_other_chrs")) {
      for (i in 1:nrow(interSV_other_chrs)) {
        chr1_in <- (interSV_other_chrs$chr1[i] %in% chr_selection$chr)
        if (!chr1_in) {
          min_pos_chr2 <- min(cnv.plot[cnv.plot$chr == interSV_other_chrs$chr2[i], "start"])
          max_pos_chr2 <- max(cnv.plot[cnv.plot$chr == interSV_other_chrs$chr2[i], "end"])
          in_range_chr2 <- (interSV_other_chrs$pos2[i] > min_pos_chr2 & interSV_other_chrs$pos2[i] < max_pos_chr2)
          if (in_range_chr2) {
            p <- p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr2[i]),
                         x = interSV_other_chrs$pos2[i], xend = interSV_other_chrs$pos2[i],
                         y = 0, yend = interSV_other_chrs$VF[i] / coeff, curvature = 0,
                         linewidth = size_interchr_line, colour = interSV_other_chrs$colour[i])
            x_range <- (chr_selection$end[which(chr_selection$chr == interSV_other_chrs$chr2[i])] -
                          chr_selection$start[which(chr_selection$chr == interSV_other_chrs$chr2[i])]) * 0.01
            p <- p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr2[i]),
                         x = interSV_other_chrs$pos2[i], xend = interSV_other_chrs$pos2[i] - x_range,
                         y = 0, yend = -yend_outside_range,
                         angle = 45, curvature = 0, linewidth = size_interchr_line, colour = interSV_other_chrs$colour[i])
          }
        }

        chr2_in <- (interSV_other_chrs$chr2[i] %in% chr_selection$chr)
        if (!chr2_in) {
          min_pos_chr1 <- min(cnv.plot[cnv.plot$chr == interSV_other_chrs$chr1[i], "start"])
          max_pos_chr1 <- max(cnv.plot[cnv.plot$chr == interSV_other_chrs$chr1[i], "end"])
          in_range_chr1 <- (interSV_other_chrs$pos1[i] > min_pos_chr1 & interSV_other_chrs$pos1[i] < max_pos_chr1)
          if (in_range_chr1) {
            p <- p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr1[i]),
                         x = interSV_other_chrs$pos1[i], xend = interSV_other_chrs$pos1[i],
                         y = 0, yend = interSV_other_chrs$VF[i] / coeff, curvature = 0,
                         linewidth = size_interchr_line, colour = interSV_other_chrs$colour[i])
            x_range <- (chr_selection$end[which(chr_selection$chr == interSV_other_chrs$chr1[i])] -
                          chr_selection$start[which(chr_selection$chr == interSV_other_chrs$chr1[i])]) * 0.01
            p <- p +
              geom_curve(data = data.frame(cov = 1, chr = interSV_other_chrs$chr1[i]),
                         x = interSV_other_chrs$pos1[i], xend = interSV_other_chrs$pos1[i] - x_range,
                         y = 0, yend = -yend_outside_range,
                         angle = 45, curvature = 0, linewidth = size_interchr_line, colour = interSV_other_chrs$colour[i])
          }
        }
      }
    }
  }

  ## Gene labels: collect all in-window genes into one data frame so labels can
  ## be de-collided together (ggrepel) instead of overprinting.
  if (!is.null(genes)) {
    gene_label_df <- data.frame()
    for (gene in genes) {
      gene_coord_now <- gene_coord[gene_coord$gene %in% gene, ]
      if (nrow(gene_coord_now) == 0) next
      if (sum(gene_coord_now$chr %in% chr_selection$chr) > 0) {
        if (sum(nrow(gene_coord_now) > 0 & gene_coord_now$chr %in% chr_selection$chr &
                gene_coord_now$start > chr_selection$start[which(chr_selection$chr %in% gene_coord_now$chr)] &
                gene_coord_now$start < chr_selection$end[which(chr_selection$chr %in% gene_coord_now$chr)]) > 0) {

          if (displayExon) {
            if (is.null(cds_gr)) {
              stop("displayExon = TRUE requires 'cds_gr' (a GRanges of CDS/exon ranges with a 'gene_name' column).")
            }
            gene_dat <- cds_gr %>% gUtils::gr2dt() %>% dplyr::filter(gene_name == gene) %>%
              dplyr::mutate(chr = paste0("chr", seqnames))
            gene_dat$pos_label <- min(gene_dat$start) + abs((min(gene_dat$start) - max(gene_dat$end))) / 2
            gene_dat$y <- max.cn * offset_gene * 1.01
            p <- p + geom_text(data = gene_dat %>% dplyr::distinct(gene_name, .keep_all = TRUE),
                               mapping = aes(x = pos_label, y = y * 1.05, label = gene_name),
                               size = size_gene_label, fontface = "italic")
            p <- p + geom_rect(data = gene_dat, aes(xmin = start, xmax = end, ymin = y * 0.93, ymax = y * 0.96), color = "grey")
            p <- p + geom_rect(data = gene_dat %>% dplyr::mutate(min_coord = min(start), max_coord = max(end)) %>%
                                 dplyr::distinct(gene_name, .keep_all = TRUE),
                               aes(xmin = min_coord, xmax = max_coord,
                                   ymin = (y * 0.93 + y * 0.96) / 2, ymax = (y * 0.93 + y * 0.96) / 2), color = "grey")
          } else {
            gene_label_df <- rbind(gene_label_df, data.frame(
              label = gene,
              chr = factor(gene_coord_now$chr, levels = chr_selection$chr),
              pos = (gene_coord_now$start + gene_coord_now$end) / 2,
              y = max.cn * offset_gene))
          }
        }
      }
    }

    if (nrow(gene_label_df) > 0) {
      ## Anchor points at the locus, then labels above.
      p <- p + geom_point(data = gene_label_df, mapping = aes(x = pos, y = y * 0.92),
                          shape = 16, size = 1, colour = "black")
      use_repel <- repel_labels && requireNamespace("ggrepel", quietly = TRUE)
      if (use_repel) {
        p <- p + ggrepel::geom_text_repel(
          data = gene_label_df, mapping = aes(x = pos, y = y, label = label),
          size = size_gene_label, fontface = "italic",
          direction = "x", segment.size = 0.2, segment.colour = "grey60",
          min.segment.length = 0, box.padding = 0.15, max.overlaps = Inf,
          seed = 1L)
      } else {
        p <- p + geom_text(data = gene_label_df, mapping = aes(x = pos, y = y, label = label),
                           size = size_gene_label, fontface = "italic")
      }
    }
  }

  ## Stress amplified segments:
  p <- p +
    geom_segment(aes(x = start, y = majorAlleleCopyNumber, xend = end, yend = majorAlleleCopyNumber, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[] & cnv.plot$copyNumber > 3 * cnv.plot$ploidy, ],
                 colour = color_major_cn, linewidth = cn_size - 0.3) +
    geom_point(aes(x = start, y = majorAlleleCopyNumber, group = chr),
               data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[] & cnv.plot$copyNumber > 3 * cnv.plot$ploidy, ],
               colour = color_major_cn, size = cn_size - 1) +
    geom_point(aes(x = end, y = majorAlleleCopyNumber, group = chr),
               data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[] & cnv.plot$copyNumber > 3 * cnv.plot$ploidy, ],
               colour = color_major_cn, size = cn_size - 1)

  ##### Tidy plot: ####
  if (verbose) {
    message("Max majorAlleleCopyNumber: ", ceiling(max(cnv.plot$majorAlleleCopyNumber)))
    if (!is.null(intraSV) && nrow(intraSV) > 0) message("Max intraSV VF: ", max(intraSV$VF))
  }

  ## Copy-number (left) axis: breaks span the actual CN range so labels do not
  ## bunch at the baseline for highly amplified samples. `yend_left` (when set)
  ## overrides the auto axis top.
  cn_axis_max <- if (!is.null(yend_left)) yend_left else nice_ceiling(max.cn)
  cn_breaks   <- unique(c(0, cn_axis_max / 2, cn_axis_max))

  ## Read-support (secondary, right) axis. Arcs sit at y = VF / coeff, so the
  ## SVs occupy the CN axis up to `arc_top`. Place the read-support ticks at
  ## real VF values (positioned via / coeff). When a focal amplicon dwarfs the
  ## SV arcs, the secondary axis is uninformative and would collide, so omit it.
  max_vf  <- if (nrow(sv) > 0) max(sv$VF, na.rm = TRUE) else 0
  arc_top <- if (coeff > 0) max_vf / coeff else 0
  show_read_support <- arc_top >= 0.12 * cn_axis_max

  p <- p +
    geom_segment(aes(x = start, xend = start, y = 0, yend = cn_axis_max, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[1], ] %>% dplyr::filter(start == min(start)),
                 colour = "black", linewidth = 0.4) +
    geom_segment(aes(x = end, xend = end, y = 0, yend = yend_right_line, group = chr),
                 data = cnv.plot[cnv.plot$chr %in% chr_selection$chr[nrow(chr_selection)], ] %>% dplyr::filter(end == max(end)),
                 colour = "black", linewidth = 0.4)

  if (show_read_support) {
    vf_step   <- nice_step(max_vf, 2)
    vf_values <- seq(0, floor(max_vf / vf_step) * vf_step, by = vf_step)
    p <- p + scale_y_continuous(
      breaks = cn_breaks,
      sec.axis = sec_axis(trans = ~ . * coeff, breaks = vf_values / coeff,
                          labels = vf_values, name = "Read support"))
  } else {
    if (verbose) message("Omitting read-support axis (SV arcs are small relative ",
                         "to the copy-number range at this locus).")
    p <- p + scale_y_continuous(breaks = cn_breaks)
  }

  ##### Save plot: ####
  outfile <- NULL
  if (!is.null(outdir) && save) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    ## Encode the plotted region in the filename so a zoomed view does not
    ## overwrite the full-chromosome PDF.
    region_tag <- if (is.null(chromosome_range)) "" else
      paste0("_", paste(sprintf("%s-%s",
                                format(round(chr_selection$start), scientific = FALSE),
                                format(round(chr_selection$end), scientific = FALSE)),
                        collapse = "_"))
    outfile <- file.path(outdir, paste0(this_sample, "_",
                                        paste(chr_selection$chr, collapse = "_"),
                                        region_tag, "_linear_plot.pdf"))
    grDevices::pdf(outfile, width = plot_width, height = plot_height, useDingbats = FALSE)
    print(p)
    grDevices::dev.off()
    if (verbose) message("Wrote ", outfile)
  }

  attr(p, "path") <- outfile
  if (!is.null(outfile)) invisible(p) else p
}
