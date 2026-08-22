#' Detect focal-amplicon seeds from copy number alone
#'
#' Derives candidate focal-amplicon regions ("seeds") directly from
#' allele-specific copy-number segments, so `call_simple_excision()` can run
#' **without** an AmpliconArchitect amplicon catalogue. Per sample, segments with
#' `copyNumber > min_cn_ratio * ploidy` are merged across gaps up to `gap` bp;
#' adjacent regions of SIMILAR copy number are then merged across larger gaps (up
#' to `merge_gap`), so a single amplicon split by internal deletions is not
#' mistaken for several separate amplicons (which would let a duplication that is
#' internal to the true amplicon pass as a boundary DUP). Regions narrower than
#' `min_width` are dropped; each surviving region is labelled with an `ID` and
#' `WGS_ID`, matching the `ecdna_gr` contract.
#'
#' @param cnv_gr A [GenomicRanges::GRanges] of allele-specific copy-number
#'   segments with metadata columns `sample`, `copyNumber`, `ploidy`.
#' @param min_cn_ratio Numeric; a segment is amplified where
#'   `copyNumber > min_cn_ratio * ploidy` (default `3`).
#' @param gap Integer; merge amplified segments separated by at most this many bp
#'   into one seed (default `1e6`).
#' @param min_width Integer; drop seeds narrower than this (default `1e5`).
#' @param merge_gap Integer; after the initial `gap` merge, also merge two
#'   adjacent same-chromosome seeds separated by at most this many bp when their
#'   representative copy numbers are similar (see `cn_ratio`). Models one amplicon
#'   broken only by internal deletions (default `3e6`). Set to `gap` to disable
#'   copy-number-aware merging.
#' @param cn_ratio Numeric in `(0, 1]`; two adjacent seeds count as "similar copy
#'   number" (and are eligible to merge across up to `merge_gap`) when
#'   `min(cn1, cn2) / max(cn1, cn2) >= cn_ratio` (default `0.5`, i.e. within 2x).
#' @param breakpoints Optional per-breakend table (a [GenomicRanges::GRanges] or
#'   data.frame with `seqnames`/`chr`, `start`/`pos`, `WGS_ID`/`sample`, `event`,
#'   `PURPLE_CN`). When supplied the seeds become **SV-aware**: two seeds on the
#'   same chromosome joined by a structural variant whose *both* breakends sit in
#'   amplified copy number (`PURPLE_CN > min_cn_ratio * ploidy`, the same bar a
#'   segment must clear to be a seed) are linked into one amplicon, even across a
#'   large gap (e.g. a centromere) -- so a single rearranged amplicon that
#'   traverses a centromere is not split. A junction that dips into non-amplified
#'   sequence on either side does not link. Each returned
#'   seed also carries `cn_only` (`TRUE` when no breakend falls inside it, i.e. it
#'   is called by copy number alone with no supporting junction).
#' @param link_tol Integer; tolerance (bp) when mapping a breakend to a seed for
#'   SV-aware linking and the `cn_only` test (default `1e4`). A boundary-defining
#'   breakend typically sits *at* the amplified-segment edge, so the reduced seed
#'   may begin a base or two after it; strict containment would then miss exactly
#'   the junction that anchors the boundary. Matching within `link_tol` (the same
#'   window used to read a breakend's `PURPLE_CN`) lets an edge-anchored breakend
#'   map to its seed.
#'
#' @return A [GenomicRanges::GRanges] of amplicon seeds with metadata columns
#'   `ID`, `WGS_ID`, and (when `breakpoints` is supplied) `cn_only` (empty if no
#'   amplification is present).
#' @seealso [call_simple_excision()]
#' @examples
#' seeds <- detect_amplicon_seeds(ex_caller_inputs$cnv_gr)
#' seeds
#' @export
#' @importFrom GenomicRanges reduce width start end seqnames findOverlaps
detect_amplicon_seeds <- function(cnv_gr, min_cn_ratio = 3,
                                  gap = 2e6, min_width = 1e5,
                                  merge_gap = 3e6, cn_ratio = 0.5,
                                  breakpoints = NULL, link_tol = 1e4) {
  stopifnot(methods::is(cnv_gr, "GRanges"),
            !is.null(cnv_gr$sample), !is.null(cnv_gr$copyNumber),
            !is.null(cnv_gr$ploidy))

  empty <- GenomicRanges::GRanges()
  S4Vectors::mcols(empty)$ID     <- character(0)
  S4Vectors::mcols(empty)$WGS_ID <- character(0)

  ## Normalise the optional breakpoint table to chr / pos / sample / event / cn.
  bp <- NULL
  if (!is.null(breakpoints)) {
    bp <- if (methods::is(breakpoints, "GRanges")) {
      data.frame(chr = as.character(GenomicRanges::seqnames(breakpoints)),
                 pos = GenomicRanges::start(breakpoints),
                 sample = breakpoints$WGS_ID, event = breakpoints$event,
                 cn = as.numeric(breakpoints$PURPLE_CN), stringsAsFactors = FALSE)
    } else {
      d <- as.data.frame(breakpoints)
      data.frame(chr = as.character(d[[if ("seqnames" %in% names(d)) "seqnames" else "chr"]]),
                 pos = d[[if ("start" %in% names(d)) "start" else "pos"]],
                 sample = d[[if ("WGS_ID" %in% names(d)) "WGS_ID" else "sample"]],
                 event = d$event, cn = as.numeric(d$PURPLE_CN), stringsAsFactors = FALSE)
    }
  }

  amp <- cnv_gr[cnv_gr$copyNumber > min_cn_ratio * cnv_gr$ploidy]
  if (length(amp) == 0L) return(empty)

  seeds <- lapply(unique(amp$sample), function(s) {
    as <- amp[amp$sample == s]
    r <- GenomicRanges::reduce(as, min.gapwidth = gap)
    r <- r[GenomicRanges::width(r) >= min_width]
    if (length(r) == 0L) return(NULL)

    ## Representative copy number of each region: the width-weighted mean copy
    ## number of the amplified segments it covers.
    ov <- GenomicRanges::findOverlaps(r, as)
    qh <- S4Vectors::queryHits(ov); sh <- S4Vectors::subjectHits(ov)
    w  <- pmin(GenomicRanges::end(as)[sh], GenomicRanges::end(r)[qh]) -
          pmax(GenomicRanges::start(as)[sh], GenomicRanges::start(r)[qh]) + 1
    num <- tapply(as$copyNumber[sh] * w, qh, sum)
    den <- tapply(w, qh, sum)
    region_cn <- rep(NA_real_, length(r))
    region_cn[as.integer(names(num))] <- num / den

    reg <- data.frame(
      chr   = as.character(GenomicRanges::seqnames(r)),
      start = GenomicRanges::start(r),
      end   = GenomicRanges::end(r),
      cn    = region_cn,
      stringsAsFactors = FALSE
    )
    reg <- reg[order(reg$chr, reg$start), , drop = FALSE]

    ## Copy-number-aware merge: fuse adjacent same-chromosome regions separated by
    ## <= merge_gap whose copy numbers are within cn_ratio. Iterate to convergence
    ## (region counts per sample are small).
    if (merge_gap > gap && nrow(reg) > 1) {
      repeat {
        merged <- FALSE
        for (i in seq_len(nrow(reg) - 1L)) {
          if (reg$chr[i] != reg$chr[i + 1L]) next
          g  <- reg$start[i + 1L] - reg$end[i]
          rr <- min(reg$cn[i], reg$cn[i + 1L]) / max(reg$cn[i], reg$cn[i + 1L])
          if (!is.na(rr) && g <= merge_gap && rr >= cn_ratio) {
            wi <- reg$end[i] - reg$start[i] + 1; wj <- reg$end[i + 1L] - reg$start[i + 1L] + 1
            reg$cn[i]  <- (reg$cn[i] * wi + reg$cn[i + 1L] * wj) / (wi + wj)
            reg$end[i] <- max(reg$end[i], reg$end[i + 1L])
            reg        <- reg[-(i + 1L), , drop = FALSE]
            merged     <- TRUE
            break
          }
        }
        if (!merged) break
      }
    }

    ## SV-aware linking: fuse two seeds on the same chromosome joined by a
    ## structural variant whose BOTH breakends sit in AMPLIFIED copy number
    ## (PURPLE_CN > min_cn_ratio * ploidy -- the same bar a segment must clear to
    ## be a seed). Such an amplicon-internal junction spans the two regions (e.g.
    ## across a centromere) as one amplicon; a junction that dips into
    ## non-amplified sequence on either side does not link -- extending across a
    ## breakend that is not itself amplified would stitch together regions that
    ## are not part of one contiguous amplicon.
    cn_only <- rep(NA, nrow(reg))
    if (!is.null(bp) && nrow(reg) >= 1) {
      pl <- stats::median(as$ploidy)
      bp_s <- bp[bp$sample == s & !is.na(bp$pos), , drop = FALSE]
      ## Map a breakend to a seed within link_tol bp (a boundary-anchoring
      ## breakend sits at the seed edge, so allow a small window rather than
      ## strict containment); if it falls in more than one, take the nearest.
      seed_of <- function(ch, p) {
        w <- which(reg$chr == ch & reg$start - link_tol <= p & reg$end + link_tol >= p)
        if (!length(w)) return(NA_integer_)
        if (length(w) == 1L) return(w)
        d <- pmax(reg$start[w] - p, p - reg$end[w], 0)
        w[which.min(d)]
      }
      if (nrow(reg) >= 2 && nrow(bp_s)) {
        parent <- seq_len(nrow(reg))
        root <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
        for (e in split(bp_s, bp_s$event)) {
          if (nrow(e) < 2 || length(unique(e$chr)) != 1) next
          if (any(is.na(e$cn)) || !all(e$cn > min_cn_ratio * pl)) next  # both breakends amplified
          sids <- unique(stats::na.omit(mapply(seed_of, e$chr, e$pos)))
          if (length(sids) >= 2) for (k in 2:length(sids)) parent[root(sids[k])] <- root(sids[1])
        }
        comp <- vapply(seq_len(nrow(reg)), root, integer(1))
        reg <- do.call(rbind, lapply(unique(comp), function(cc) {
          rr <- reg[comp == cc, , drop = FALSE]
          data.frame(chr = rr$chr[1], start = min(rr$start), end = max(rr$end),
                     cn = stats::weighted.mean(rr$cn, rr$end - rr$start + 1), stringsAsFactors = FALSE)
        }))
        reg <- reg[order(reg$chr, reg$start), , drop = FALSE]
      }
      ## cn_only: no breakend falls inside the (final) seed (within link_tol, so
      ## an edge-anchored junction still counts as supporting the seed)
      cn_only <- vapply(seq_len(nrow(reg)), function(i)
        !any(bp_s$chr == reg$chr[i] & bp_s$pos >= reg$start[i] - link_tol &
             bp_s$pos <= reg$end[i] + link_tol),
        logical(1))
    }

    out <- GenomicRanges::GRanges(reg$chr, IRanges::IRanges(reg$start, reg$end))
    S4Vectors::mcols(out)$WGS_ID  <- s
    S4Vectors::mcols(out)$ID      <- paste0(s, "_amp", seq_along(out))
    S4Vectors::mcols(out)$cn_only <- cn_only
    out
  })
  seeds <- Filter(Negate(is.null), seeds)
  if (length(seeds) == 0L) return(empty)
  do.call(c, unname(seeds))
}
