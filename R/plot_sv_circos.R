## ---------------------------------------------------------------------------
## Whole-genome circos plot for one sample: ideogram + copy-number-status ring,
## intra- and inter-chromosomal breakpoint-density rings, and class-coloured SV
## links. Self-contained (hg38 by default): the ideogram cytoband and genome
## tiling are built from EpiTracer's bundled chr_info, so no network or
## GitHub-only packages (regioneR / gUtils) are needed -- only circlize.
## ---------------------------------------------------------------------------

## Circlize-ready cytoband (chr/start/end/name/gieStain) from bundled chr_info.
.load_cytoband <- function(genome) {
  f <- system.file("extdata", paste0("chr_info_", genome, ".rds"), package = "EpiTracer")
  if (!nzchar(f)) stop("No bundled cytoband for genome '", genome, "'.", call. = FALSE)
  ci <- as.data.frame(readRDS(f))
  data.frame(chr = as.character(ci$seqnames), start = as.integer(ci$start),
             end = as.integer(ci$end), name = as.character(ci$name),
             gieStain = as.character(ci$gieStain), stringsAsFactors = FALSE)
}

#' Whole-genome SV circos plot for one sample
#'
#' Draws a per-sample circos: chromosome-label + banded ideogram, a copy-number
#' status ring (amplification / gain / LOH / homozygous deletion / neutral,
#' relative to ploidy), separate intra- and inter-chromosomal breakpoint-density
#' rings, and structural-variant links coloured by class (DUP / TRA / h2hINV /
#' t2tINV / other). Titled with the sample and its whole-genome-doubling status.
#'
#' A port of a bespoke circos plotter to EpiTracer conventions: hg38 by default,
#' the ideogram and genome tiling built from the bundled `chr_info`, and the
#' gUtils/regioneR helpers replaced by the package's own [to_granges()] / `\%$\%`.
#' `chrY` and `chrM` are dropped. Requires the \pkg{circlize} package.
#'
#' @param sample Character scalar; the sample to plot (matched against the
#'   `sample` column of `sv_data` / `cnv_data`).
#' @param sv_data A `data.frame`/`data.table` of structural variants in BEDPE
#'   layout: `sample`, `chrom1`, `start1`, `chrom2`, `start2`, `svclass` (and
#'   optionally `end1`/`end2`; when absent the start positions are used).
#' @param cnv_data A `data.frame`/`data.table` of copy-number segments: `sample`,
#'   `seqnames`, `start`, `end`, `copyNumber`, `ploidy`, `minorAlleleCopyNumber`.
#' @param wgd_data Optional `data.frame` with `sample` and `Polyploidy`
#'   (`"wgd"` vs anything else) for the title; `NULL` labels the sample only.
#' @param genome Reference build; one of `"hg38"` (default), `"hg19"`, `"mm10"`.
#' @param bin_size Integer; genome tile width (bp) for the density and
#'   copy-number rings (default `1e6`).
#' @param highlight_events Optional character vector of SV identifiers to draw
#'   bold on top of the faint full-genome links (e.g. the junctions of one
#'   chromoplexy cycle). Matched against `highlight_id_col`.
#' @param highlight_id_col Optional name of the column in `sv_data` holding the
#'   identifiers matched by `highlight_events`; `NULL` (default) auto-detects
#'   `name` then `event`.
#' @param highlight_colour Colour for highlighted links (default `"#d95f0e"`).
#' @param dim_unhighlighted Logical; when `TRUE`, non-highlighted links are greyed
#'   so the highlighted set stands out (default `FALSE`).
#' @param outdir Optional directory; when supplied the plot is written to
#'   `<outdir>/<sample>_circos.pdf`. When `NULL` the plot is drawn on the current
#'   graphics device (open one yourself, e.g. with [grDevices::pdf()]).
#' @param overwrite Logical; when `FALSE` and the output PDF already exists, the
#'   sample is skipped (default `TRUE`).
#' @return Invisibly `NULL`; called for its plotting side effect.
#' @seealso [plot_sv_linear()]
#' @export
plot_sv_circos <- function(sample, sv_data, cnv_data, wgd_data = NULL,
                           genome = c("hg38", "hg19", "mm10"), bin_size = 1e6,
                           highlight_events = NULL, highlight_id_col = NULL,
                           highlight_colour = "#d95f0e", dim_unhighlighted = FALSE,
                           outdir = NULL, overwrite = TRUE) {
  genome <- match.arg(genome)
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("plot_sv_circos() requires the 'circlize' package.", call. = FALSE)
  .s <- sample                                       # avoid the sample/column clash

  outfile <- if (!is.null(outdir)) file.path(outdir, paste0(.s, "_circos.pdf")) else NULL
  if (!is.null(outfile) && !overwrite && file.exists(outfile)) return(invisible(NULL))

  pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  drop_chr <- c("chrY", "chrM")

  cytoband <- .load_cytoband(genome)
  cytoband <- cytoband[!cytoband$chr %in% drop_chr, ]
  seqlens <- load_chrom_lengths(genome)
  seqlens <- seqlens[!names(seqlens) %in% drop_chr]

  ## ---- SV data ------------------------------------------------------------
  sv <- data.table::as.data.table(sv_data)[sample == .s]
  if (!nrow(sv)) stop("No structural variants for sample '", .s, "'.", call. = FALSE)
  sv[, `:=`(chrom1 = pfx(as.character(chrom1)), chrom2 = pfx(as.character(chrom2)))]
  sv <- sv[!chrom1 %in% drop_chr & !chrom2 %in% drop_chr]
  if (!"end1" %in% names(sv)) sv[, end1 := start1]
  if (!"end2" %in% names(sv)) sv[, end2 := start2]
  sv[, cols := data.table::fcase(
    svclass == "DUP", "#c1447e", svclass == "TRA", "#fac881",
    svclass == "h2hINV", "#a5a6ae", svclass == "t2tINV", "#384351",
    default = "#bde0fe")]
  ## optional event highlight (shared contract with the linear plotters): matched
  ## junctions are re-drawn bold on top; dim_unhighlighted greys the rest first.
  hl <- .resolve_highlight(sv, highlight_events, highlight_id_col)
  if (isTRUE(dim_unhighlighted) && any(hl)) sv[!hl, cols := "grey85"]
  nuc1 <- as.data.frame(sv[, .(chrom1, start1, end1)])
  nuc2 <- as.data.frame(sv[, .(chrom2, start2, end2)])

  ## ---- genome tiling + breakpoint-density rings ---------------------------
  ## cut.last.tile.in.chrom = TRUE returns a GRanges (not a GRangesList)
  bins <- GenomicRanges::tileGenome(seqlens, tilewidth = bin_size,
                                    cut.last.tile.in.chrom = TRUE)
  bp_gr <- function(d) to_granges(data.frame(
    seqnames = c(d$chrom1, d$chrom2), start = c(d$start1, d$start2),
    end = c(d$end1, d$end2)))
  density_ring <- function(d) {
    dens <- if (nrow(d)) GenomicRanges::countOverlaps(bins, bp_gr(d)) else integer(length(bins))
    data.frame(chr = as.character(GenomeInfoDb::seqnames(bins)),
               start = GenomicRanges::start(bins), end = GenomicRanges::end(bins),
               density = dens)
  }
  intra_ring <- density_ring(sv[svclass != "TRA"])
  inter_ring <- density_ring(sv[svclass == "TRA"])

  ## ---- copy-number status ring --------------------------------------------
  cnv <- data.table::as.data.table(cnv_data)[sample == .s]
  cnv[, seqnames := pfx(as.character(seqnames))]
  cnv <- cnv[!seqnames %in% drop_chr]
  binned <- gr2dt(bins %$% to_granges(as.data.frame(cnv)))
  binned[, status := data.table::fcase(
    copyNumber > 3 * ploidy, 5,          # amplification
    copyNumber > 1.4 * ploidy, 4,        # gain
    copyNumber < 0.5, 1,                 # homozygous deletion
    minorAlleleCopyNumber < 0.5, 2,      # LOH
    default = 3)]                        # neutral
  binned[is.na(status), status := 3]
  cn_ring <- as.data.frame(binned[, .(chr = seqnames, start, end, status)])

  ## ---- WGD status for the title -------------------------------------------
  ## Use the supplied wgd_data if given; otherwise derive WGD from the copy
  ## number itself with PURPLE's rule (see call_wgd()).
  wgd_status <- "NA"
  if (!is.null(wgd_data)) {
    w <- data.table::as.data.table(wgd_data)[sample == .s]
    if (nrow(w)) wgd_status <- if (identical(as.character(w$Polyploidy[1]), "wgd")) "WGD" else "Diploid"
  } else {
    ww <- .wgd_from_dt(cnv)
    if (nrow(ww)) wgd_status <- if (isTRUE(ww$wgd[1])) "WGD" else "Diploid"
  }

  ## ---- draw ---------------------------------------------------------------
  if (!is.null(outfile)) { grDevices::pdf(outfile); on.exit(grDevices::dev.off(), add = TRUE) }
  circlize::circos.clear()
  circlize::circos.par(gap.degree = 1, start.degree = 90, cell.padding = c(0, 0, 0, 0))
  circlize::circos.initializeWithIdeogram(cytoband = cytoband, plotType = "labels")
  circlize::circos.genomicIdeogram(cytoband = cytoband, track.height = 0.03)

  circlize::circos.genomicHeatmap(
    bed = cn_ring, numeric.column = "status", heatmap_height = 0.04, side = "inside",
    connection_height = NULL,
    col = circlize::colorRamp2(c(1, 2, 3, 4, 5),
      c("#0077b6", "#bde0fe", "#e0e1dd", "#d1aac2", "#A30000")))
  circlize::circos.genomicHeatmap(
    bed = intra_ring, numeric.column = "density", heatmap_height = 0.04, side = "inside",
    connection_height = NULL,
    col = circlize::colorRamp2(c(0, 1, 5, 10, 50, 100),
      c("#e0e1dd", "#a4ac86", "#656d4a", "#414833", "#333d29", "#051923")))
  circlize::circos.genomicHeatmap(
    bed = inter_ring, numeric.column = "density", heatmap_height = 0.04, side = "inside",
    connection_height = NULL,
    col = circlize::colorRamp2(c(0, 1, 5, 10, 50, 100),
      c("#e0e1dd", "#e0b1cb", "#be95c4", "#9f86c0", "#5e548e", "#231942")))

  circlize::circos.genomicLink(nuc1, nuc2, col = sv$cols, h.ratio = 0.65)
  if (any(hl))                                    # highlighted junctions, bold, on top
    circlize::circos.genomicLink(nuc1[hl, , drop = FALSE], nuc2[hl, , drop = FALSE],
                                 col = highlight_colour, h.ratio = 0.65, lwd = 2.5)
  graphics::title(paste0(.s, " (", wgd_status, ")"))
  invisible(NULL)
}
