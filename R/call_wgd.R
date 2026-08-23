## ---------------------------------------------------------------------------
## Whole-genome doubling (WGD) from allele-specific copy number, per PURPLE's
## definition: a WGD is called when the MAJOR-allele copy number exceeds 1.5
## across at least half of the bases on at least 11 of the 22 autosomes.
## ---------------------------------------------------------------------------

## Per-sample WGD from a copy-number data.table (seqnames/start/end/sample/
## majorAlleleCopyNumber). Returns sample, n_wgd_autosomes, wgd.
.wgd_from_dt <- function(cn) {
  pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  cn <- data.table::as.data.table(cn)
  if (!"majorAlleleCopyNumber" %in% names(cn))
    stop("call_wgd() needs a 'majorAlleleCopyNumber' column.", call. = FALSE)
  if (!"sample" %in% names(cn)) cn[, sample := "sample"]   # single-sample fallback
  cn[, seqnames := pfx(as.character(seqnames))]
  cn <- cn[seqnames %in% paste0("chr", 1:22)]                 # autosomes only
  if (!nrow(cn))
    return(data.table::data.table(sample = character(), n_wgd_autosomes = integer(), wgd = logical()))
  cn[, w := as.numeric(end) - as.numeric(start) + 1]
  ## fraction of each autosome's bases with major-allele CN > 1.5
  per_chr <- cn[, .(frac = sum(w[majorAlleleCopyNumber > 1.5], na.rm = TRUE) /
                            sum(w, na.rm = TRUE)),
                by = .(sample, seqnames)]
  res <- per_chr[, .(n_wgd_autosomes = sum(frac >= 0.5, na.rm = TRUE)), by = sample]
  res[, wgd := n_wgd_autosomes >= 11L]
  res[]
}

#' Call whole-genome doubling (WGD) from allele-specific copy number
#'
#' Applies PURPLE's whole-genome-doubling rule: a sample is WGD when its
#' **major-allele copy number exceeds 1.5 across at least half of the bases** on
#' **at least 11 of the 22 autosomes**. Operates directly on the allele-specific
#' copy-number segments EpiTracer already uses, so WGD status needs no separate
#' input file.
#'
#' @param cnv_gr A [GenomicRanges::GRanges], `data.frame` or
#'   [data.table::data.table] of copy-number segments carrying
#'   `majorAlleleCopyNumber` (and `seqnames`/`start`/`end`, taken from the ranges
#'   for a GRanges). An optional `sample` column splits multi-sample input; when
#'   absent the whole table is treated as one sample. Chromosome names may be with
#'   or without the `chr` prefix.
#' @return A [data.table::data.table] with one row per `sample`:
#'   `n_wgd_autosomes` (how many autosomes meet the >1.5-over-half-bases test) and
#'   `wgd` (logical; `n_wgd_autosomes >= 11`).
#' @seealso [plot_sv_circos()]
#' @export
call_wgd <- function(cnv_gr) {
  cn <- if (methods::is(cnv_gr, "GRanges")) gr2dt(cnv_gr) else data.table::as.data.table(cnv_gr)
  .wgd_from_dt(cn)
}
