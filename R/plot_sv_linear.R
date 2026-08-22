#' Linear copy-number and structural rearrangement plot
#'
#' Draws allele-specific copy number, structural variant arcs, karyotype
#' ideograms and gene labels for one or more loci laid out side-by-side on a
#' single concatenated x-axis. Because all loci share one coordinate system, the
#' structural variants that interconnect separate amplicons (e.g. the junctions
#' of a multi-fragment / hub ecDNA) are drawn as arcs spanning the loci.
#'
#' Optionally, supplying `snv_data` adds a second panel of the sample's small
#' mutations directly beneath the copy-number / SV panel, sharing the same
#' concatenated genomic x-axis so each mutation lines up with the copy number
#' and rearrangements it sits within. Its y-axis shows the intermutation
#' distance (a rainfall plot, the default), the variant-allele frequency, or the
#' SNV (mutation) copy number (see `snv_y`), and the points can be coloured by
#' their timing relative to the focal amplification (`snv_timing`). The two
#' panels are returned and saved as one stacked figure (requires \pkg{patchwork});
#' see the `snv_*` parameters for the full set of controls.
#'
#' The loci to display are resolved in this order:
#' \enumerate{
#'   \item `loci` if supplied (explicit windows);
#'   \item `chromosome` + `chromosome_range` (explicit windows, one per
#'     chromosome);
#'   \item `chromosome` alone -- the amplified region on each named chromosome is
#'     auto-detected and padded by `flank_pct`; a chromosome with no amplification
#'     falls back to its whole length;
#'   \item nothing -- every amplified locus in the sample is detected
#'     automatically.
#' }
#' A locus is "amplified" where `copyNumber > min_cn_ratio * ploidy`; its
#' extent is the min-max of amplified segments on the chromosome, padded by
#' `flank_pct`, with single-segment artefacts below `min_amp_width` dropped.
#'
#' @param sample Character scalar; sample identifier (matched in `cnv_data`,
#'   `sv_data`, `wgd_data`).
#' @param cnv_data A data.frame of copy-number segments with columns `sample`,
#'   `seqnames`, `start`, `end`, `copyNumber`, `ploidy`, `majorAlleleCopyNumber`,
#'   `minorAlleleCopyNumber`.
#' @param sv_data A data.frame of SVs with columns `chrom1`, `start1`, `chrom2`,
#'   `start2`, `strand1`, `strand2`, `svclass`, `VF`, `JCN`, `sample`.
#' @param wgd_data Optional data.frame with a sample-identifier column (see
#'   `wgd_sample_col`) and a `Polyploidy` column (`"No"` = diploid, otherwise
#'   WGD). If supplied, the sample's WGD status is annotated in the plot title;
#'   if `NULL` (default) the title shows just the sample name.
#' @param genome Genome build for the bundled references: `"hg38"` (default),
#'   `"hg19"` or `"mm10"`. Selects the karyotype ideogram and oncogene panel used
#'   when `karyotype` / `gene_coord` are not given. Ignored for either reference
#'   that is supplied explicitly.
#' @param karyotype Ideogram bands (UCSC `cytoBand` / `chr_info` style), as a
#'   data.frame (first column = chromosome; needs `start`, `end`, `gieStain`) or
#'   an `.rds` path. `NULL` (default) uses the bundled reference for `genome`.
#' @param gene_coord Gene coordinates as a data.frame with columns
#'   `chr`,`start`,`end`,`strand`,`gene`, or a path to a headerless tab-separated
#'   BED-like file with those columns. `NULL` (default) uses the bundled oncogene
#'   panel for `genome`; supply your own to label other genes.
#' @param chromosome Optional character vector of chromosomes to display, e.g.
#'   `c("chr7", "chr12")`.
#' @param chromosome_range Optional two-column matrix/data.frame of `start`,`end`
#'   window limits, one row per entry in `chromosome`.
#' @param loci Optional explicit loci: a `data.frame` with columns
#'   `chr`,`start`,`end`, or a character vector of `"chr:start-end"` strings.
#'   Takes precedence over `chromosome`/`chromosome_range`.
#' @param events Character vector of copy-number event types to target when
#'   auto-detecting loci (used only when neither `loci` nor `chromosome_range` is
#'   given). Any of: `"amp"` (amplification, `copyNumber > min_cn_ratio * ploidy`),
#'   `"gain"` (`> gain_ratio * ploidy` but not amplified), `"loh"`
#'   (`minorAlleleCopyNumber < loh_thresh`), `"homdel"` (homozygous deletion,
#'   `copyNumber < homdel_thresh`). Default `"amp"`. When explicit loci are
#'   supplied the whole region is plotted regardless of event type -- the
#'   function is a general CN/SV viewer, not amplification-only.
#' @param cluster_gap Numeric; when auto-detecting, consecutive event segments
#'   more than this many bp apart start a new locus, so scattered focal events
#'   become separate panels instead of one whole-chromosome span (default `5e6`).
#' @param flank_pct Numeric percentage by which each auto-detected (CN-status)
#'   region is extended on both sides -- i.e. the flanking window shown around an
#'   amplicon / deletion, as a percent of its width (default `10` = +/-10%). Only
#'   applies to auto-detected loci, not to explicit `chromosome`/`loci`.
#' @param min_cn_ratio Numeric; amplification threshold as a multiple of ploidy
#'   (`copyNumber > min_cn_ratio * ploidy`, default `3`).
#' @param gain_ratio Numeric; gain threshold as a multiple of ploidy (default
#'   `1.4`); a "gain" is `> gain_ratio * ploidy` but not amplified.
#' @param loh_thresh Numeric; LOH minor-allele threshold
#'   (`minorAlleleCopyNumber < loh_thresh`, default `0.5`).
#' @param homdel_thresh Numeric; homozygous-deletion copy-number threshold
#'   (`copyNumber < homdel_thresh`, default `0.5`).
#' @param min_amp_width Numeric; drop auto-detected loci whose total event span
#'   is below this many bp (default `1e5`).
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
#'   `NULL` it is rounded up from the data to two significant figures, with a
#'   floor of 2 so near-diploid regions still use a full 0-2 axis.
#' @param displayExon Logical; if `TRUE`, draw exon models (from `cds_gr`) for
#'   in-window genes instead of a point + label (default `FALSE`).
#' @param cds_gr Optional [GenomicRanges::GRanges] of CDS/exon ranges (metadata
#'   column `gene_name`); required when `displayExon = TRUE`.
#' @param offset_gene,ymax_highlight_ratio,karyotype_rel_size,loh_position_ratio
#'   Layout tuning parameters.
#' @param highlight_amp,highlight_hom_del Logical; shade amplified / homozygously
#'   deleted segments.
#' @param amplicons Optional distinct-amplicon overlay: a
#'   [GenomicRanges::GRanges] or data.frame of amplicon regions (columns
#'   `seqnames`/`chr`, `start`, `end`, and an optional `ID`/`label`). Each
#'   amplicon that falls in the plotted window is drawn as a short horizontal bar
#'   across the TOP of the plot spanning its extent, in its own colour and with
#'   its label above, so several distinct amplicons (e.g. a focal episome sitting
#'   inside a larger, separately detected amplified span) are visually separable
#'   without washing over the data. Genomic coordinates are mapped through the
#'   same per-locus transform as the rest of the plot, so the bars align with the
#'   copy-number track.
#' @param parallel_breakpoints Optional data.frame of *adjacent parallel
#'   breakpoint* pairs to highlight (the breakage-replication/fusion hallmark):
#'   columns `chr`/`seqnames`, `pos1`, `pos2` (the two same-orientation breakends
#'   of a pair) and an optional `strand`. Each pair is drawn near the baseline as
#'   a bracket joining the two breakends, with a caret marking each, so BRF
#'   breakpoints stand out among the other junctions.
#' @param wgd_sample_col Optional name of the sample column in `wgd_data`
#'   (default: `sample`, falling back to `WGS_ID`).
#' @param snv_data Optional SNV/SSM table as a `data.frame` or
#'   [GenomicRanges::GRanges] with (at least) `seqnames`/`start` position columns,
#'   a sample-identifier column (see `snv_sample_col`) and, for `snv_y = "vaf"`, a
#'   variant-allele-frequency column (see `vaf_col`). When supplied, a second SNV
#'   panel (intermutation-distance rainfall by default, or VAF; see `snv_y`) is
#'   drawn directly beneath the copy-number / SV panel, sharing the same
#'   concatenated genomic x-axis, and the two are returned/saved as one stacked
#'   figure (requires the \pkg{patchwork} package). When `NULL` (default) only the
#'   CN/SV panel is drawn and behaviour is unchanged.
#' @param snv_sample_col Optional name of the sample column in `snv_data`
#'   (default: `sampleID`, falling back to `sample`).
#' @param snv_y What the SNV panel's y-axis shows: `"imd"` (default) plots the
#'   intermutation distance -- the bp distance to the previous SNV on the same
#'   chromosome, a rainfall plot, on a log10 axis; `"vaf"` plots the
#'   variant-allele frequency; `"cn"` plots the SNV copy number (mutation copy
#'   number, see `snv_cn_col`) -- in an amplicon this times each SNV against the
#'   amplification. Only single-nucleotide variants are used in every case (indels
#'   / MNVs are excluded, see `snv_type_col`); intermutation distances are computed
#'   across all of the sample's SNVs per chromosome before restricting to the
#'   plotted loci, so window-edge mutations keep their true neighbour distance.
#' @param snv_type_col Optional name of a mutation-type column in `snv_data` used
#'   to keep SNVs only (rows where the value is `"SNV"`). Defaults to `type` when
#'   present; if no such column exists, SNVs are inferred from single-base
#'   `ref`/`mut` columns.
#' @param vaf_col Name of the VAF column in `snv_data` (default `allelic_freq`),
#'   used when `snv_y = "vaf"`. Values outside `[0, vaf_max]` are treated as
#'   artefacts and dropped.
#' @param snv_cn_col Name of the SNV copy-number column in `snv_data` (default
#'   `variant_cn`), used when `snv_y = "cn"` and for `snv_timing`.
#' @param snv_timing Logical; if `TRUE`, classify each SNV by its timing relative
#'   to the focal amplification and colour the points accordingly, with a legend
#'   collected to the right of the figure. Uses the mutation copy number
#'   (`snv_cn_col`) against the amplified-allele copy number (`major_cn`) at the
#'   site (also needs `minor_cn`): a site is amplified when
#'   `major_cn + minor_cn > min_cn_ratio * ploidy`; within it, an SNV is
#'   "Pre-amplification" when its copy number is `>= snv_timing_pre_frac * major_cn`
#'   (and `>= 2`), "Post-amplification" when it is `<= snv_timing_post_mcn`, and
#'   "Unknown" otherwise or when the site is not amplified. Works with any `snv_y`.
#' @param snv_timing_pre_frac Numeric; fraction of the amplified-allele copy number
#'   an SNV's copy number must reach to be called pre-amplification (default `0.5`).
#' @param snv_timing_post_mcn Numeric; SNV copy-number at or below which a variant
#'   in an amplified site is called post-amplification (default `1.5`).
#' @param snv_timing_colours Named character vector of point colours for the
#'   `"Pre-amplification"`, `"Post-amplification"` and `"Unknown"` classes.
#' @param snv_rel_height Numeric; height of the SNV panel relative to the CN/SV
#'   panel (default `0.4`).
#' @param snv_point_size,snv_alpha,snv_colour Point size, alpha and colour for the
#'   SNV VAF scatter (defaults `0.7`, `0.6`, `"#1d3557"`).
#' @param vaf_max Numeric top of the VAF axis; also the upper bound of the plausible
#'   VAF window used to drop artefactual values (default `1`).
#' @param outdir Optional directory in which to write the plot. If `NULL` no file
#'   is written and only the plot object is returned.
#' @param save Logical; if `TRUE` (default) and `outdir` is supplied, write the
#'   plot as a PDF.
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
#' @seealso [call_simple_excision()]
#' @export
#' @import ggplot2
#' @importFrom grid unit
#' @importFrom data.table rbindlist
#' @importFrom grDevices pdf dev.off
#' @importFrom utils read.table head
plot_sv_linear <- function(sample,
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
                           genes_to_highlight = NULL,
                           gene_label_angle = NULL,
                           repel_labels = TRUE,
                           cn_max = NULL,
                           displayExon = FALSE,
                           cds_gr = NULL,
                           offset_gene = 1.15,
                           ymax_highlight_ratio = 1.08,
                           karyotype_rel_size = 0.048,
                           loh_position_ratio = 0.5,
                           highlight_amp = TRUE,
                           highlight_hom_del = TRUE,
                           amplicons = NULL,
                           parallel_breakpoints = NULL,
                           wgd_sample_col = NULL,
                           snv_data = NULL,
                           snv_sample_col = NULL,
                           snv_y = c("imd", "vaf", "cn"),
                           snv_type_col = NULL,
                           vaf_col = "allelic_freq",
                           snv_cn_col = "variant_cn",
                           snv_timing = FALSE,
                           snv_timing_pre_frac = 0.5,
                           snv_timing_post_mcn = 1.5,
                           snv_timing_colours = c("Pre-amplification"  = "#d1495b",
                                                  "Post-amplification" = "#1d3557",
                                                  "Unknown"            = "grey70"),
                           snv_rel_height = 0.4,
                           snv_point_size = 0.7,
                           snv_alpha = 0.6,
                           snv_colour = "#1d3557",
                           vaf_max = 1,
                           outdir = NULL,
                           save = TRUE,
                           plot_width_custom = NULL,
                           plot_height_custom = NULL,
                           verbose = FALSE) {

  this_sample <- sample

  ## Whether a second (SNV VAF) panel is drawn beneath the CN/SV panel. When it
  ## is, the karyotype ideogram is dropped from the CN panel (the inter-panel gap
  ## stays clean) and the shared genomic axis (Mb ticks + chromosome names) is
  ## moved to the bottom of the SNV panel, labelling the figure once at the base.
  has_snv <- !is.null(snv_data)

  ## ---- Reference data -----------------------------------------------------
  ## `genome` selects the bundled karyotype + oncogene panel (hg38 / hg19 /
  ## mm10); an explicit `karyotype` / `gene_coord` (data.frame or path) overrides.
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
  ## default to the whole loaded panel, so labels adapt to `genome` (the mm10
  ## panel uses mouse symbols, which a hardcoded human list would never match).
  genes <- if (is.null(genes_to_highlight)) unique(gene_coord$gene) else genes_to_highlight

  ## ---- Sample data --------------------------------------------------------
  cnv <- cnv_data[cnv_data$sample == this_sample, , drop = FALSE]
  cnv$chr <- as.character(cnv$seqnames)
  if (nrow(cnv) == 0) stop("No copy-number data for sample ", this_sample, ".")

  ## WGD annotation is optional: only shown in the title when `wgd_data` is given.
  wgd_status <- NULL
  if (!is.null(wgd_data)) {
    wgd_data <- as.data.frame(wgd_data)
    wgd_col <- if (!is.null(wgd_sample_col)) wgd_sample_col else
      if ("sample" %in% names(wgd_data)) "sample" else "WGS_ID"
    wrow <- wgd_data[wgd_data[[wgd_col]] == this_sample, , drop = FALSE]
    wgd_status <- if (nrow(wrow) >= 1 && wrow$Polyploidy[1] == "No") "Diploid" else "WGD"
  }

  ## ---- Resolve loci -------------------------------------------------------
  events <- match.arg(events, c("amp", "gain", "loh", "homdel"), several.ok = TRUE)
  ## Copy-number event predicate over a set of segments:
  ##   amp    : copyNumber > min_cn_ratio * ploidy
  ##   gain   : copyNumber > gain_ratio  * ploidy, but not amplified
  ##   loh    : minorAlleleCopyNumber < loh_thresh
  ##   homdel : copyNumber < homdel_thresh
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
    ## Cluster event segments per chromosome: start a new locus wherever
    ## consecutive events are separated by more than `cluster_gap`, so scattered
    ## focal events (e.g. homozygous deletions) become separate loci rather than
    ## one whole-chromosome span.
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
      if (nrow(agg) == 0) stop("No amplified (", paste(events, collapse = "/"), ") loci found for ",
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
  ## `nice_ceiling()` lives in R/utils-axis.R (shared with plot_sv_reconstruction).
  ## Left (copy-number) axis top: rounded up from the data, but never below 2 so a
  ## near-diploid region still uses a full 0-2 axis instead of a squashed 0-1.1 one.
  cn_axis_max <- if (!is.null(cn_max)) cn_max else max(2, nice_ceiling(max.cn))

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
  ## Right (read-support / VF) axis. When the loci contain no SVs there is no VF to
  ## scale to, but the axis is still drawn (for a consistent two-axis layout) by
  ## mirroring the copy-number axis (coeff = 1), so the top panel always shows both.
  rs_axis_max <- if (max_vf > 0) nice_ceiling(max_vf) else cn_axis_max
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

  ## ---- Distinct-amplicon bands (background) -------------------------------
  ## Draw each supplied amplicon as its own full-height translucent band so
  ## multiple amplicons in one window (e.g. a focal episome enclosed by a larger
  ## detected span) are visually separable. Mapped through the same per-locus
  ## `offset + (pos - locus_start)` transform as every other layer.
  if (!is.null(amplicons)) {
    amp_df <- if (methods::is(amplicons, "GRanges")) {
      data.frame(chr = as.character(GenomicRanges::seqnames(amplicons)),
                 start = GenomicRanges::start(amplicons),
                 end = GenomicRanges::end(amplicons),
                 label = if (!is.null(amplicons$ID)) as.character(amplicons$ID)
                         else as.character(seq_along(amplicons)),
                 stringsAsFactors = FALSE)
    } else {
      dd <- as.data.frame(amplicons)
      chrc <- if ("chr" %in% names(dd)) dd$chr else dd$seqnames
      lab  <- if ("ID" %in% names(dd)) dd$ID else if ("label" %in% names(dd)) dd$label
              else seq_len(nrow(dd))
      data.frame(chr = as.character(chrc), start = dd$start, end = dd$end,
                 label = as.character(lab), stringsAsFactors = FALSE)
    }
    amp_df$chr <- ifelse(grepl("^chr", amp_df$chr), amp_df$chr, paste0("chr", amp_df$chr))
    amp_bands <- do.call(rbind, lapply(seq_len(nrow(amp_df)), function(k) {
      do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
        if (amp_df$chr[k] != loci$chr[i]) return(NULL)
        s <- max(amp_df$start[k], loci$start[i]); e <- min(amp_df$end[k], loci$end[i])
        if (e < s) return(NULL)
        data.frame(gx_start = loci$offset[i] + (s - loci$start[i]),
                   gx_end   = loci$offset[i] + (e - loci$start[i]),
                   label = amp_df$label[k], ki = k, stringsAsFactors = FALSE)
      }))
    }))
    if (!is.null(amp_bands) && nrow(amp_bands) > 0) {
      amp_pal <- c("#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756",
                   "#72B7B2", "#EECA3B", "#9D755D")
      amp_bands$fill <- amp_pal[((amp_bands$ki - 1L) %% length(amp_pal)) + 1L]
      ## Draw each amplicon as a short horizontal bar across the TOP of the plot
      ## (spanning its extent), with its label above -- a compact extent marker
      ## rather than a full-height wash over the data. Placed ABOVE the gene-label
      ## row (which sits at cn_axis_max * offset_gene) so gene annotations stay
      ## clear; coord_cartesian(clip = "off") keeps it visible above the panel.
      amp_y <- cn_axis_max * (offset_gene + 0.12)
      p <- p +
        geom_segment(data = amp_bands, aes(x = gx_start, xend = gx_end,
                  y = amp_y, yend = amp_y),
                  colour = amp_bands$fill, linewidth = 1.8, lineend = "butt", alpha = 0.55) +
        geom_text(data = amp_bands, aes(x = (gx_start + gx_end) / 2,
                  y = amp_y, label = label),
                  vjust = -0.5, size = 2.3, fontface = "bold", colour = amp_bands$fill, alpha = 0.7)
    }
  }

  ## ---- Adjacent parallel breakpoints (BRF hallmark) ----------------------
  ## Highlight each same-orientation breakend pair with a bracket near the top of
  ## the copy-number track and a caret at each breakend.
  if (!is.null(parallel_breakpoints) && nrow(as.data.frame(parallel_breakpoints)) > 0) {
    pbp <- as.data.frame(parallel_breakpoints)
    pbp$chr <- as.character(if ("chr" %in% names(pbp)) pbp$chr else pbp$seqnames)
    pbp$chr <- ifelse(grepl("^chr", pbp$chr), pbp$chr, paste0("chr", pbp$chr))
    gxof <- function(ch, pos) {
      for (i in seq_len(nrow(loci)))
        if (ch == loci$chr[i] && pos >= loci$start[i] && pos <= loci$end[i])
          return(loci$offset[i] + (pos - loci$start[i]))
      NA_real_
    }
    pb_rows <- do.call(rbind, lapply(seq_len(nrow(pbp)), function(k) {
      g1 <- gxof(pbp$chr[k], pbp$pos1[k]); g2 <- gxof(pbp$chr[k], pbp$pos2[k])
      if (is.na(g1) || is.na(g2)) return(NULL)
      data.frame(g1 = g1, g2 = g2,
                 strand = if ("strand" %in% names(pbp)) as.character(pbp$strand[k]) else "",
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(pb_rows) && nrow(pb_rows) > 0) {
      ## Mark each parallel breakpoint with a black down-arrowhead across the TOP
      ## of the plot, just above the amplicon-bar row so the two overlays coexist.
      pb_y <- cn_axis_max * (offset_gene + 0.20)
      p <- p +
        geom_point(data = data.frame(gx = c(pb_rows$g1, pb_rows$g2)),
                   aes(x = gx, y = pb_y), shape = 25, size = 1.8,
                   fill = "black", colour = "black")
    }
  }

  ideo <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    d <- karyotype[karyotype$chr == loci$chr[i] & karyotype$end >= loci$start[i] & karyotype$start <= loci$end[i], , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d$start <- pmax(d$start, loci$start[i]); d$end <- pmin(d$end, loci$end[i])
    d$gx_start <- loci$offset[i] + (d$start - loci$start[i])
    d$gx_end   <- loci$offset[i] + (d$end   - loci$start[i])
    d
  }))
  if (!has_snv && !is.null(ideo) && nrow(ideo) > 0)
    p <- p + geom_rect(data = ideo, aes(xmin = gx_start, xmax = gx_end,
             ymin = lower_limit_karyotype, ymax = upper_limit_karyotype),
             fill = ideo$color, colour = "black", linewidth = 0.2)

  ## Telomere markers: label a window edge "tel" when it reaches a chromosome
  ## terminus -- the p-telomere at coordinate ~0 or the q-telomere at the
  ## chromosome length -- so a reader can tell whether an amplicon's end sits at
  ## an actual telomere.
  tel_tol <- 5e5
  tel_rows <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
    L <- suppressWarnings(max(karyotype$end[karyotype$chr == loci$chr[i]]))
    rr <- NULL
    if (loci$start[i] <= tel_tol)
      rr <- rbind(rr, data.frame(gx = loci$offset[i]))                     # p-telomere
    if (is.finite(L) && (L - loci$end[i]) <= tel_tol)
      rr <- rbind(rr, data.frame(gx = loci$offset[i] + loci$width[i]))     # q-telomere
    rr
  }))
  if (!has_snv && !is.null(tel_rows) && nrow(tel_rows) > 0)
    p <- p + geom_text(data = tel_rows, aes(x = gx, y = lower_limit_karyotype),
                       label = "tel", vjust = 1.5, size = 2.4, fontface = "bold",
                       colour = "#0b6b6b")

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
      } else if (gx1 <= gx2) -abs(curv_inter) else abs(curv_inter)  # TRA: always bow upward regardless of chrom1/chrom2 order
      ## draw the connecting arc only when the two breakends map to distinct x
      ## positions (geom_curve errors on identical end points); the breakpoint
      ## verticals are drawn regardless.
      if (abs(gx2 - gx1) >= 1)
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

  ## Gene labels (in-window). `displayExon` draws exon models; otherwise a
  ## point + de-collided ggrepel label. The gene track is anchored to the axis top
  ## (`cn_axis_max`), not the local data max, so genes sit at the top of the panel
  ## regardless of the underlying allelic CN (e.g. over a near-diploid region).
  gene_y_ref <- cn_axis_max
  if (displayExon) {
    if (is.null(cds_gr)) stop("displayExon = TRUE requires 'cds_gr' (a GRanges of ",
                              "CDS/exon ranges with a 'gene_name' column).")
    cds <- as.data.frame(cds_gr)
    cds$chr <- if (all(grepl("^chr", cds$seqnames))) as.character(cds$seqnames)
               else paste0("chr", cds$seqnames)
    yb <- gene_y_ref * offset_gene
    ## Minimum rendered exon width so short exons remain visible at genomic scale:
    min_ex_w <- total_w * 0.0015
    ex_df <- data.frame(); body_df <- data.frame(); lab_df <- data.frame()
    for (i in seq_len(nrow(loci))) {
      for (gn in intersect(genes, cds$gene_name[cds$chr == loci$chr[i]])) {
        ex <- cds[cds$gene_name == gn & cds$chr == loci$chr[i] &
                  cds$end >= loci$start[i] & cds$start <= loci$end[i], , drop = FALSE]
        if (nrow(ex) == 0) next
        gxs <- loci$offset[i] + (pmax(ex$start, loci$start[i]) - loci$start[i])
        gxe <- loci$offset[i] + (pmin(ex$end,   loci$end[i])   - loci$start[i])
        ## widen exons narrower than min_ex_w symmetrically about their centre:
        short <- (gxe - gxs) < min_ex_w
        mid <- (gxs + gxe) / 2
        gxs[short] <- mid[short] - min_ex_w / 2
        gxe[short] <- mid[short] + min_ex_w / 2
        ex_df   <- rbind(ex_df, data.frame(gxs = gxs, gxe = gxe))
        body_df <- rbind(body_df, data.frame(gxs = min(gxs), gxe = max(gxe)))
        lab_df  <- rbind(lab_df, data.frame(gx = (min(gxs) + max(gxe)) / 2, gene = gn))
      }
    }
    if (nrow(ex_df) > 0) {
      p <- p +
        geom_segment(data = body_df, aes(x = gxs, xend = gxe, y = yb * 0.945, yend = yb * 0.945),
                     colour = "grey40", linewidth = 0.4) +
        geom_rect(data = ex_df, aes(xmin = gxs, xmax = gxe, ymin = yb * 0.90, ymax = yb * 0.99),
                  fill = "grey25", colour = "grey15", linewidth = 0.1) +
        geom_text(data = lab_df, aes(x = gx, y = yb * 1.08, label = gene),
                  size = size_gene_label, fontface = "italic")
    }
  } else {
    gl <- do.call(rbind, lapply(seq_len(nrow(loci)), function(i) {
      g <- gene_coord[gene_coord$chr == loci$chr[i] &
                      gene_coord$start >= loci$start[i] & gene_coord$end <= loci$end[i] &
                      gene_coord$gene %in% genes, , drop = FALSE]
      if (nrow(g) == 0) return(NULL)
      g$gx <- loci$offset[i] + ((g$start + g$end) / 2 - loci$start[i])
      g
    }))
    if (!is.null(gl) && nrow(gl) > 0) {
      gl$y <- gene_y_ref * offset_gene
      ## "Crowded" = at least two gene labels closer than 5% of the plotted width.
      ## When crowded, labels are angled and ggrepel draws a leader line back to
      ## each point; when the labels are nicely separated, the leader lines are
      ## suppressed (min.segment.length = Inf) as they only add clutter.
      gp <- sort(gl$gx)
      crowded <- length(gp) > 1 && min(diff(gp)) < 0.05 * total_w
      if (is.null(gene_label_angle)) gene_label_angle <- if (crowded) 45 else 0
      p <- p + geom_point(data = gl, aes(x = gx, y = y * 0.90), shape = 16, size = 1, colour = "black")
      if (repel_labels && requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(data = gl, aes(x = gx, y = y, label = gene),
          size = size_gene_label, fontface = "italic", angle = gene_label_angle, hjust = 0,
          direction = "both", nudge_y = gene_y_ref * 0.02, ylim = c(gene_y_ref * 0.99, NA),
          segment.size = 0.2, segment.colour = "grey60",
          min.segment.length = if (crowded) 0 else Inf,
          box.padding = 0.25, max.overlaps = Inf, seed = 1L)
      } else {
        p <- p + geom_text(data = gl, aes(x = gx, y = y, label = gene),
                           size = size_gene_label, fontface = "italic", angle = gene_label_angle)
      }
    }
  }

  ## ---- Bottom labelling (Mb ticks + Mb labels + chromosome names) ---------
  unit_y      <- cn_axis_max * 0.07
  tick_y0     <- lower_limit_karyotype
  tick_y1     <- lower_limit_karyotype - 0.7 * unit_y   # longer x tick marks
  mb_label_y  <- lower_limit_karyotype - 1.6 * unit_y
  chr_label_y <- lower_limit_karyotype - 3.0 * unit_y

  ## Minimum spacing (in gx units) between Mb tick labels so they never overlap
  ## on narrow panels; derived from the output width so labels stay legible.
  plot_width_est <- if (!is.null(plot_width_custom)) plot_width_custom else max(5, 2.2 * nrow(loci) + 1)
  min_tick_gap <- 0.62 * total_w / plot_width_est
  tick_df <- data.frame(); mb_df <- data.frame()
  for (i in seq_len(nrow(loci))) {
    ticks <- pretty(c(loci$start[i], loci$end[i]), n = 3)
    ticks <- ticks[ticks >= loci$start[i] & ticks <= loci$end[i]]
    ticks <- ticks[ticks != 0]                          # drop a "0.0" at a locus start
    if (length(ticks) == 0) next
    gxt <- loci$offset[i] + (ticks - loci$start[i])
    ## greedily thin ticks that are closer than min_tick_gap:
    keep <- rep(FALSE, length(gxt)); last <- -Inf
    for (j in seq_along(gxt)) if (gxt[j] - last >= min_tick_gap) { keep[j] <- TRUE; last <- gxt[j] }
    gxt <- gxt[keep]; ticks <- ticks[keep]
    tick_df <- rbind(tick_df, data.frame(x = gxt))
    mb_df <- rbind(mb_df, data.frame(x = gxt, lab = formatC(ticks / 1e6, format = "f", digits = 1)))
  }
  chr_lab <- data.frame(gx = loci$gx_mid, chr = gsub("chr", "Chr ", loci$chr))
  cn_breaks <- unique(c(0, cn_axis_max / 2, cn_axis_max))

  ## Reusable genomic-axis decoration: draws the Mb tick marks, Mb labels and
  ## chromosome names below `y_base` (in that panel's own y units, `y_unit`).
  draw_genomic_axis <- function(gg, y_base, y_unit) {
    gg +
      geom_segment(data = tick_df, aes(x = x, xend = x, y = y_base, yend = y_base - 0.7 * y_unit),
                   colour = "black", linewidth = 0.4) +
      geom_text(data = mb_df, aes(x = x, y = y_base - 1.6 * y_unit, label = lab),
                size = size_text / 2.1, vjust = 1) +
      geom_text(data = chr_lab, aes(x = gx, y = y_base - 3.0 * y_unit, label = chr),
                size = size_text / 1.78, vjust = 1, fontface = "bold")
  }

  p <- p +
    geom_segment(data = data.frame(x = loci$gx_start[1]),
                 aes(x = x, xend = x, y = 0, yend = cn_axis_max), colour = "black", linewidth = 0.4) +
    geom_segment(data = data.frame(x = loci$gx_end[nrow(loci)]),
                 aes(x = x, xend = x, y = 0, yend = cn_axis_max), colour = "black", linewidth = 0.4)
  ## The CN panel carries the genomic axis only when there is no SNV panel below.
  if (!has_snv) p <- draw_genomic_axis(p, y_base = lower_limit_karyotype, y_unit = unit_y)
  p <- p +
    ggtitle(if (is.null(wgd_status)) this_sample else paste0(this_sample, " (", wgd_status, ")")) +
    labs(x = NULL, y = "Allele specific\ncopy number") +
    coord_cartesian(clip = "off", expand = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    theme(
      text = element_text(size = size_text, colour = "black"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_text(size = size_text + 3, colour = "black"),
      axis.text.y.right = element_text(size = size_text + 3, colour = "black"),
      axis.title.y = element_text(size = size_text + 4, colour = "black"),
      axis.title.y.right = element_text(size = size_text + 4, colour = "black"),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.y.right = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length.y = unit(0.2, "cm"),           # longer y tick marks
      panel.background = element_blank(), plot.background = element_blank(), panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = size_text + 4),
      ## Drop the tall bottom margin when the SNV panel (which now carries the
      ## genomic axis labels) sits below.
      plot.margin = unit(c(.3, .5, if (has_snv) .1 else 1.3, .2), "cm")
    )

  ## Always carry the read-support secondary axis (see `coeff` above for the
  ## no-SV case), so the panel's layout is identical with or without SVs.
  rs_breaks <- unique(c(0, rs_axis_max / 2, rs_axis_max))
  p <- p + scale_y_continuous(breaks = cn_breaks,
    sec.axis = sec_axis(trans = ~ . * coeff, breaks = rs_breaks, name = "Read support"))

  ## ---- SNV panel ----------------------------------------------------------
  ## A second track of the sample's small mutations, sharing the CN/SV panel's
  ## concatenated genomic x-axis (same `loci` offsets + `locus_idx`), so mutations
  ## line up beneath the copy-number and rearrangement they sit within. The y-axis
  ## is either the intermutation distance ("imd", a rainfall plot -- the default)
  ## or the variant-allele frequency ("vaf").
  plot_obj <- p
  if (has_snv) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Drawing the SNV panel requires the 'patchwork' package; install it or set snv_data = NULL.")
    snv_y <- match.arg(snv_y, c("imd", "vaf", "cn"))
    snv <- as.data.frame(snv_data)
    scol <- if (!is.null(snv_sample_col)) snv_sample_col else
      if ("sampleID" %in% names(snv)) "sampleID" else "sample"
    if (!scol %in% names(snv)) stop("SNV sample column '", scol, "' not found in `snv_data`.")
    snv <- snv[snv[[scol]] == this_sample, , drop = FALSE]
    ## Keep single-nucleotide variants only (exclude indels / MNVs): use the
    ## mutation-type column when present, else infer from single-base ref/alt.
    tcol <- if (!is.null(snv_type_col)) snv_type_col else if ("type" %in% names(snv)) "type" else NA
    n_all <- nrow(snv)
    if (!is.na(tcol) && tcol %in% names(snv)) {
      snv <- snv[toupper(as.character(snv[[tcol]])) == "SNV", , drop = FALSE]
    } else if (all(c("ref", "mut") %in% names(snv))) {
      snv <- snv[nchar(as.character(snv$ref)) == 1 & nchar(as.character(snv$mut)) == 1, , drop = FALSE]
    }
    snv$chr <- as.character(snv$seqnames)
    snv$pos <- as.integer(snv$start)

    ## Optional amplification-timing classification of each SNV, from its mutation
    ## copy number `m` (snv_cn_col) relative to the amplified-allele copy number
    ## (`major_cn`) at that site:
    ##   pre-amplification  -- on (nearly) all amplified copies: m >= pre_frac*major
    ##                         and m >= 2 (so it predates the amplification);
    ##   post-amplification -- on a single copy: m <= post_mcn (arose afterwards);
    ##   unknown            -- site not amplified (major+minor <= min_cn_ratio*ploidy)
    ##                         or an intermediate m (ambiguous stepwise timing).
    ## Computed here, before any y-mode row filtering, so the label rides along.
    if (snv_timing) {
      need <- c(snv_cn_col, "major_cn", "minor_cn")
      miss <- need[!need %in% names(snv)]
      if (length(miss))
        stop("snv_timing = TRUE needs column(s) not in `snv_data`: ", paste(miss, collapse = ", "))
      m_cn  <- suppressWarnings(as.numeric(snv[[snv_cn_col]]))
      maj   <- suppressWarnings(as.numeric(snv$major_cn))
      minr  <- suppressWarnings(as.numeric(snv$minor_cn))
      amp   <- is.finite(maj) & is.finite(minr) & (maj + minr) > min_cn_ratio * ploidy
      pre   <- amp & is.finite(m_cn) & is.finite(maj) & m_cn >= 2 & m_cn >= snv_timing_pre_frac * maj
      post  <- amp & !pre & is.finite(m_cn) & m_cn <= snv_timing_post_mcn
      cls <- rep("Unknown", nrow(snv))
      cls[pre]  <- "Pre-amplification"
      cls[post] <- "Post-amplification"
      ## Level order (and hence legend order) follows the colour vector.
      snv$timing <- factor(cls, levels = names(snv_timing_colours))
    }

    if (snv_y == "imd") {
      ## Intermutation distance: distance (bp) to the previous SNV on the same
      ## chromosome, computed across ALL of this sample's SNVs (not just in-window
      ## ones), so window-edge mutations get their true neighbour distance. Plotted
      ## as log10 on a linear axis, which keeps the below-panel ideogram/axis (drawn
      ## at negative y) valid -- a genuine log scale cannot render y <= 0.
      snv <- snv[order(snv$chr, snv$pos), , drop = FALSE]
      snv$imd <- stats::ave(snv$pos, snv$chr, FUN = function(z) c(NA, diff(z)))
      snv$yv  <- log10(snv$imd)
      snv <- snv[is.finite(snv$yv), , drop = FALSE]   # drop first-per-chr + imd <= 0
      y_title <- "Intermutation\ndistance (bp)"
    } else if (snv_y == "cn") {
      ## SNV copy number: the number of tumour genome copies carrying the mutation
      ## (mutation copy number). In an amplicon this times each SNV against the
      ## amplification -- pre-amplification mutations ride up to high copy number,
      ## post-amplification ones stay near 1.
      if (!snv_cn_col %in% names(snv)) stop("SNV copy-number column '", snv_cn_col,
                                            "' not found in `snv_data`.")
      snv$yv <- suppressWarnings(as.numeric(snv[[snv_cn_col]]))
      snv <- snv[is.finite(snv$yv) & snv$yv >= 0, , drop = FALSE]
      y_title <- "SNV copy\nnumber"
    } else {
      if (!vaf_col %in% names(snv)) stop("VAF column '", vaf_col, "' not found in `snv_data`.")
      snv$yv <- suppressWarnings(as.numeric(snv[[vaf_col]]))
      snv <- snv[snv$yv >= 0 & snv$yv <= vaf_max, , drop = FALSE]   # drop VAF artefacts
      y_title <- "SNV VAF"
    }

    ## Restrict to the plotted loci and map to the shared x coordinate.
    snv$li <- locus_idx(snv$chr, snv$pos)
    snv <- snv[!is.na(snv$li), , drop = FALSE]
    snv$gx <- loci$offset[snv$li] + (snv$pos - loci$start[snv$li])
    if (verbose) {
      message(sprintf("SNVs in loci: %d of %d small mutations for %s (y = %s)",
                      nrow(snv), n_all, this_sample, snv_y))
      if (snv_timing && nrow(snv))
        message("  timing: ", paste(sprintf("%s=%d", names(table(snv$timing)),
                                            as.integer(table(snv$timing))), collapse = ", "))
    }

    ## y-axis geometry (`y_top`, breaks, labels) per mode.
    if (snv_y == "imd") {
      y_top <- max(8, if (nrow(snv)) ceiling(max(snv$yv)) else 8)   # log10 bp; >= 100 Mb headroom
      y_breaks <- seq(0, y_top, by = 2)
      y_labels <- vapply(10^y_breaks, function(x)
        if (x < 1e3) formatC(x, format = "d")
        else if (x < 1e6) paste0(x / 1e3, "kb")
        else if (x < 1e9) paste0(x / 1e6, "Mb")
        else paste0(x / 1e9, "Gb"), character(1))
    } else if (snv_y == "cn") {
      y_top <- nice_ceiling(if (nrow(snv)) max(snv$yv) else 1)
      y_breaks <- unique(c(0, y_top / 2, y_top))
      y_labels <- ggplot2::waiver()
    } else {
      y_top <- vaf_max
      y_breaks <- unique(c(0, vaf_max / 2, vaf_max))
      y_labels <- ggplot2::waiver()
    }

    y_unit <- y_top * 0.07
    ## Karyotype ideogram strip beneath the SNV panel: as the bottom-most track it
    ## carries the shared genomic axis (Mb ticks + chromosome names) directly below
    ## it, mirroring how the ideogram sat at the CN panel's base in the no-SNV case.
    ## A clear gap (`0.9 * y_unit`) separates the panel's x-axis (y = 0) from the
    ## ideogram so it does not sit flush against the bottom of the plot.
    ideo_top    <- -0.9 * y_unit
    ideo_bottom <- ideo_top - y_top * karyotype_rel_size
    p_snv <- ggplot()
    ## Faint amplified-region shading, matching the top panel, for visual anchor.
    amp_shade <- cn_plot[cn_plot$copyNumber > 3 * cn_plot$ploidy, , drop = FALSE]
    if (nrow(amp_shade) > 0)
      p_snv <- p_snv + geom_rect(data = amp_shade, aes(xmin = gx_start, xmax = gx_end,
        ymin = 0, ymax = y_top), fill = "#d92a05", alpha = 0.1)
    if (nrow(snv) > 0) {
      if (snv_timing) {
        ## Colour by amplification-timing class, with a legend (collected to the
        ## right of the stacked figure below, so panel widths stay aligned). A
        ## zero-size, fully transparent seed point per class guarantees every class
        ## appears in the legend even when a class has no SNVs in view.
        seed <- data.frame(gx = loci$gx_start[1], yv = 0,
                           timing = factor(names(snv_timing_colours),
                                           levels = names(snv_timing_colours)))
        p_snv <- p_snv +
          geom_point(data = seed, aes(x = gx, y = yv, colour = timing),
                     size = 0, alpha = 0, na.rm = TRUE) +
          geom_point(data = snv, aes(x = gx, y = yv, colour = timing),
                     size = snv_point_size, alpha = snv_alpha, stroke = 0) +
          scale_colour_manual(name = "SNV timing", values = snv_timing_colours, drop = FALSE,
                              guide = guide_legend(override.aes = list(size = 2, alpha = 1)))
      } else {
        p_snv <- p_snv + geom_point(data = snv, aes(x = gx, y = yv),
          colour = snv_colour, size = snv_point_size, alpha = snv_alpha, stroke = 0)
      }
    }
    p_snv <- p_snv +
      geom_segment(data = data.frame(x = loci$gx_start[1]),
                   aes(x = x, xend = x, y = 0, yend = y_top), colour = "black", linewidth = 0.4) +
      geom_segment(data = data.frame(x = loci$gx_end[nrow(loci)]),
                   aes(x = x, xend = x, y = 0, yend = y_top), colour = "black", linewidth = 0.4)
    if (!is.null(ideo) && nrow(ideo) > 0)
      p_snv <- p_snv + geom_rect(data = ideo, aes(xmin = gx_start, xmax = gx_end,
        ymin = ideo_bottom, ymax = ideo_top), fill = ideo$color, colour = "black", linewidth = 0.2)
    p_snv <- draw_genomic_axis(p_snv, y_base = ideo_bottom, y_unit = y_unit)
    p_snv <- p_snv +
      labs(x = NULL, y = y_title) +
      ## Fix the y-range with coord (not scale limits) so the genomic-axis labels
      ## drawn just below 0 are kept and rendered into the bottom margin
      ## (clip = "off"); scale `limits` would drop them before clipping applies.
      coord_cartesian(ylim = c(0, y_top), clip = "off", expand = FALSE) +
      scale_x_continuous(expand = expansion(mult = 0.01)) +
      scale_y_continuous(breaks = y_breaks, labels = y_labels) +
      theme(
        text = element_text(size = size_text, colour = "black"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        axis.text.y = element_text(size = size_text + 3, colour = "black"),
        axis.title.y = element_text(size = size_text + 4, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
        axis.ticks.length.y = unit(0.2, "cm"),
        panel.background = element_blank(), plot.background = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(.1, .5, 1.8, .2), "cm")
      )
    plot_obj <- patchwork::wrap_plots(p, p_snv, ncol = 1, heights = c(1, snv_rel_height),
                                      guides = if (snv_timing) "collect" else "keep")
    if (snv_timing)
      ## Collect the timing legend to the right of the whole figure (not inside the
      ## SNV panel), so both panels keep the same plotting width and stay aligned.
      plot_obj <- plot_obj & theme(legend.position = "right",
                                   legend.title = element_text(size = size_text + 2, colour = "black"),
                                   legend.text  = element_text(size = size_text + 1, colour = "black"),
                                   legend.key = element_blank())
  }

  ## ---- Size + save --------------------------------------------------------
  plot_width  <- if (!is.null(plot_width_custom)) plot_width_custom else max(5, 2.2 * nrow(loci) + 1)
  plot_height <- if (!is.null(plot_height_custom)) plot_height_custom else 3.0
  ## Grow the canvas to accommodate the stacked SNV panel (unless overridden).
  if (has_snv && is.null(plot_height_custom))
    plot_height <- plot_height * (1 + snv_rel_height) + 0.3

  outfile <- NULL
  if (!is.null(outdir) && save) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    region_tag <- if (nrow(loci) == 1)
      sprintf("_%.0f-%.0f", round(loci$start[1]), round(loci$end[1])) else ""
    outfile <- file.path(outdir, paste0(this_sample, "_", paste(loci$chr, collapse = "_"),
                                        region_tag, "_linear_plot.pdf"))
    grDevices::pdf(outfile, width = plot_width, height = plot_height, useDingbats = FALSE)
    print(plot_obj); grDevices::dev.off()
    if (verbose) message("Wrote ", outfile)
  }
  attr(plot_obj, "path") <- outfile
  if (!is.null(outfile)) invisible(plot_obj) else plot_obj
}
