#' Detect focal-amplicon seeds from copy number alone
#'
#' Derives candidate focal-amplicon regions ("seeds") directly from
#' allele-specific copy-number segments, so `call_episomal_ecdna()` can run
#' **without** an AmpliconArchitect amplicon catalogue. Per sample, segments with
#' `copyNumber > min_cn_ratio * ploidy` are merged across gaps up to `gap` bp and
#' regions narrower than `min_width` are dropped; each surviving region is
#' labelled with an `ID` and `WGS_ID`, matching the `ecdna_gr` contract.
#'
#' @param cnv_gr A [GenomicRanges::GRanges] of allele-specific copy-number
#'   segments with metadata columns `sample`, `copyNumber`, `ploidy`.
#' @param min_cn_ratio Numeric; a segment is amplified where
#'   `copyNumber > min_cn_ratio * ploidy` (default `3`).
#' @param gap Integer; merge amplified segments separated by at most this many bp
#'   into one seed (default `1e6`).
#' @param min_width Integer; drop seeds narrower than this (default `1e5`).
#'
#' @return A [GenomicRanges::GRanges] of amplicon seeds with metadata columns
#'   `ID` and `WGS_ID` (empty if no amplification is present).
#' @seealso [call_episomal_ecdna()]
#' @examples
#' seeds <- detect_amplicon_seeds(ex_caller_inputs$cnv_gr)
#' seeds
#' @export
#' @importFrom GenomicRanges reduce width
detect_amplicon_seeds <- function(cnv_gr, min_cn_ratio = 3,
                                  gap = 1e6, min_width = 1e5) {
  stopifnot(methods::is(cnv_gr, "GRanges"),
            !is.null(cnv_gr$sample), !is.null(cnv_gr$copyNumber),
            !is.null(cnv_gr$ploidy))

  empty <- GenomicRanges::GRanges()
  S4Vectors::mcols(empty)$ID     <- character(0)
  S4Vectors::mcols(empty)$WGS_ID <- character(0)

  amp <- cnv_gr[cnv_gr$copyNumber > min_cn_ratio * cnv_gr$ploidy]
  if (length(amp) == 0L) return(empty)

  seeds <- lapply(unique(amp$sample), function(s) {
    r <- GenomicRanges::reduce(amp[amp$sample == s], min.gapwidth = gap)
    r <- r[GenomicRanges::width(r) >= min_width]
    if (length(r) == 0L) return(NULL)
    S4Vectors::mcols(r)$WGS_ID <- s
    S4Vectors::mcols(r)$ID     <- paste0(s, "_amp", seq_along(r))
    r
  })
  seeds <- Filter(Negate(is.null), seeds)
  if (length(seeds) == 0L) return(empty)
  do.call(c, unname(seeds))
}
