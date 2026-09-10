## ---------------------------------------------------------------------------
## Episomal ecDNA from AmpliconSuite (AmpliconArchitect + AmpliconClassifier).
##
## This function uses AmpliconSuite output as distributed by AmpliconRepository to call episomal ## ecDNA:
##   * an AmpliconArchitect `*_graph.txt` (sequence edges = CN segments;
##     discordant edges = junctions with predicted copy number and read-pair
##     support), and
##   * the AmpliconClassifier ecDNA interval BED (`*_ecDNA_<k>_intervals.bed`),
##     which defines the amplicon extent -- so the amplicon is AA's own call, not
##     re-derived by a copy-number threshold here.
##
## An amplicon is EPISOMAL if a BOUNDARY DUP exists: an intra-chromosomal
## duplication junction whose two ends bracket the outer boundary of the ecDNA
## interval on its chromosome. The boundary DUP is the highest read-pair-support
## spanning DUP (a selector, not a gate). A multi-chromosome ecDNA is a CHIMERIC
## episome when every chromosome is bracketed by its own boundary DUP and the
## chromosomes are joined by >= 1 translocation. Internal junctions are ignored:
## an episome may be internally rearranged (chromothripsis within a micronucleus).
## ---------------------------------------------------------------------------

## Coordinate-invariant SV class from a junction's two breakends. AmpliconArchitect
## writes each discordant edge in an arbitrary direction, so strands must be read
## AFTER sorting the two breakends by coordinate (lo first): lo`-`/hi`+` is a
## head-to-tail duplication (DUP, the circulariser), `+`/`-` is a DEL, matching
## strands an inversion (INV), and different chromosomes a translocation (TRA).
.aa_svclass <- function(cA, pA, sA, cB, pB, sB) {
  lo_s <- ifelse(pA <= pB, sA, sB)
  hi_s <- ifelse(pA <= pB, sB, sA)
  data.table::fifelse(cA != cB, "TRA",
    data.table::fifelse(lo_s == hi_s, "INV",
      data.table::fifelse(lo_s == "-" & hi_s == "+", "DUP", "DEL")))
}

#' Read an AmpliconArchitect graph file
#'
#' Parses an AmpliconArchitect `*_graph.txt` into copy-number segments (from the
#' `sequence` edges) and structural-variant junctions (from the `discordant`
#' edges), the latter with a **coordinate-invariant** SV class. Chromosome names
#' are accepted with or without a `chr` prefix.
#'
#' @param graph_file Path to an AmpliconArchitect `*_graph.txt`.
#'
#' @return A list with three [data.table::data.table]s:
#' \describe{
#'   \item{`segments`}{`chr`, `start`, `end`, `cn` (predicted copy number).}
#'   \item{`junctions`}{`cA`,`pA`,`sA` and `cB`,`pB`,`sB` (the two breakends'
#'     chromosome/position/strand), `jcn` (junction predicted copy number),
#'     `pe_support` (number of supporting read pairs), `homlen` (breakpoint
#'     homology length, `NA` if absent), `svclass` (`DUP`/`DEL`/`INV`/`TRA`) and
#'     `lo`/`hi` (the lower/upper breakend coordinate).}
#'   \item{`source`}{`chr`, `pos`, `cn` for each AA *source* edge -- an open end
#'     where copy-number flow enters/leaves the amplicon via the rest of the
#'     genome instead of looping back through a junction (a self-contained circle
#'     has none). Used by the structural QC of
#'     [call_episomal_from_ampliconsuite()].}
#' }
#' Any table may have zero rows.
#' @seealso [call_episomal_from_ampliconsuite()]
#' @export
#' @importFrom data.table data.table tstrsplit fifelse
read_ampliconsuite <- function(graph_file) {
  seg <- data.table::data.table(chr = character(), start = numeric(),
                                end = numeric(), cn = numeric())
  jun <- data.table::data.table(cA = character(), pA = numeric(), sA = character(),
                                cB = character(), pB = numeric(), sB = character(),
                                jcn = numeric(), pe_support = numeric(),
                                homlen = numeric(), svclass = character(),
                                lo = numeric(), hi = numeric())
  src <- data.table::data.table(chr = character(), pos = numeric(), cn = numeric())
  if (is.null(graph_file) || is.na(graph_file) || !file.exists(graph_file))
    return(list(segments = seg, junctions = jun, source = src))
  L <- readLines(graph_file, warn = FALSE)

  ## sequence edges -> CN segments
  se <- L[startsWith(L, "sequence")]
  if (length(se)) {
    f  <- data.table::tstrsplit(se, "\t", fixed = TRUE)
    a  <- regmatches(f[[2]], regexec("^([^:]+):([0-9]+)", f[[2]]))
    b  <- regmatches(f[[3]], regexec(":([0-9]+)", f[[3]]))
    seg <- data.table::data.table(
      chr   = vapply(a, function(x) if (length(x) == 3) x[2] else NA_character_, ""),
      start = as.numeric(vapply(a, function(x) if (length(x) == 3) x[3] else NA_character_, "")),
      end   = as.numeric(vapply(b, function(x) if (length(x) == 2) x[2] else NA_character_, "")),
      cn    = suppressWarnings(as.numeric(f[[4]])))
    seg <- seg[!is.na(chr) & !is.na(cn)]
  }

  ## discordant edges -> junctions
  d <- L[startsWith(L, "discordant")]
  if (length(d)) {
    g  <- data.table::tstrsplit(d, "\t", fixed = TRUE)
    sp <- data.table::tstrsplit(g[[2]], "->", fixed = TRUE)
    pl <- regmatches(sp[[1]], regexec("^([^:]+):(-?[0-9]+)([+-])$", sp[[1]]))
    pr <- regmatches(sp[[2]], regexec("^([^:]+):(-?[0-9]+)([+-])$", sp[[2]]))
    gg <- function(m, i) vapply(m, function(x) if (length(x) == 4) x[i] else NA_character_, "")
    hom <- if (length(g) >= 5) suppressWarnings(as.numeric(g[[5]])) else NA_real_
    jun <- data.table::data.table(
      cA = gg(pl, 2), pA = as.numeric(gg(pl, 3)), sA = gg(pl, 4),
      cB = gg(pr, 2), pB = as.numeric(gg(pr, 3)), sB = gg(pr, 4),
      jcn = suppressWarnings(as.numeric(g[[3]])),
      pe_support = suppressWarnings(as.numeric(g[[4]])),
      homlen = hom)
    jun <- jun[!is.na(cA) & !is.na(cB)]
    if (nrow(jun)) {
      jun[is.na(pe_support), pe_support := 0]
      jun[, svclass := .aa_svclass(cA, pA, sA, cB, pB, sB)]
      jun[, lo := pmin(pA, pB)][, hi := pmax(pA, pB)]
    }
  }

  ## source edges -> open ends. AA writes 'source  chr:-1-->chr:pos(strand)  CN ...';
  ## a source edge is copy-number flow that enters/leaves the amplicon via the rest
  ## of the genome (an OPEN end) rather than looping back through a junction. The
  ## real breakend is the last "->"-separated field; `cn` (col 3) is the open flow.
  so <- L[startsWith(L, "source")]
  if (length(so)) {
    s  <- data.table::tstrsplit(so, "\t", fixed = TRUE)
    sp <- data.table::tstrsplit(s[[2]], "->", fixed = TRUE)
    be <- sp[[length(sp)]]                              # 2nd endpoint = actual breakend
    m  <- regmatches(be, regexec("^([^:]+):(-?[0-9]+)([+-])$", be))
    src <- data.table::data.table(
      chr = vapply(m, function(x) if (length(x) == 4) x[2] else NA_character_, ""),
      pos = as.numeric(vapply(m, function(x) if (length(x) == 4) x[3] else NA_character_, "")),
      cn  = suppressWarnings(as.numeric(s[[3]])))
    src <- src[!is.na(chr) & !is.na(pos)]
  }
  list(segments = seg, junctions = jun, source = src)
}

## Estimate tumour ploidy from a genome-wide copy-number BED (e.g. the AA-input
## `*_CNV_CALLS.bed`, CNVkit format `chr start end method CN`). Ploidy is the
## width-weighted median copy number across the whole genome; the CN column is
## taken as the last numeric column (so both 4- and 5-column BEDs work). Returns
## NA when the file is missing/unreadable.
.cnv_ploidy <- function(cnv_file) {
  if (is.null(cnv_file) || length(cnv_file) != 1L || is.na(cnv_file) ||
      !file.exists(cnv_file)) return(NA_real_)
  x <- tryCatch(data.table::fread(cnv_file, header = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) < 1L || ncol(x) < 4L) return(NA_real_)
  cn_col <- suppressWarnings(max(which(vapply(x, function(c) is.numeric(c) || !any(is.na(suppressWarnings(as.numeric(as.character(c))))), logical(1)))))
  if (!is.finite(cn_col) || cn_col < 4L) return(NA_real_)
  v <- suppressWarnings(as.numeric(x[[cn_col]]))
  w <- suppressWarnings(as.numeric(x[[3]]) - as.numeric(x[[2]]))
  k <- is.finite(v) & is.finite(w) & w > 0
  if (!any(k)) return(NA_real_)
  o <- order(v[k]); as.numeric(v[k][o][which(cumsum(w[k][o]) / sum(w[k]) >= 0.5)[1]])
}

## Read an AmpliconClassifier ecDNA interval BED into inclusive-coordinate
## intervals. AC writes 1-based starts and an exclusive end (end = last base + 1),
## so the inclusive end is `end - 1`; this makes the boundary-DUP bracket test
## exact against graph/junction coordinates with no fuzzy tolerance.
.read_ac_intervals <- function(intervals) {
  if (is.character(intervals) && length(intervals) == 1L) {
    b <- data.table::fread(intervals, header = FALSE)
    b <- b[, 1:3]; data.table::setnames(b, c("chr", "start", "end"))
  } else {
    b <- data.table::as.data.table(intervals)
    nm <- names(b); pick <- function(o) o[o %in% nm][1]
    b <- b[, c(pick(c("chr", "seqnames", "chrom")),
               pick(c("start", "chromStart")),
               pick(c("end", "chromEnd"))), with = FALSE]
    data.table::setnames(b, c("chr", "start", "end"))
  }
  b[, `:=`(chr = as.character(chr), start = as.numeric(start), end = as.numeric(end) - 1)]
  b[!is.na(chr) & !is.na(start) & !is.na(end)]
}

#' Derive ecDNA intervals from an AmpliconArchitect cycles file
#'
#' Generalises [call_episomal_from_ampliconsuite()] to **any** AmpliconArchitect
#' output, including bare AA runs without AmpliconClassifier. The ecDNA extent is
#' taken from the reconstructed cycle: the segments of the dominant (highest
#' `Copy_count`) cycle, or of all cycles. The AA "source" vertex (segment id `0`)
#' is dropped. Returned coordinates use the exclusive-end convention of the
#' AmpliconClassifier BED, so the result plugs straight into the `intervals`
#' argument of [call_episomal_from_ampliconsuite()].
#'
#' @param cycles_file Path to an AmpliconArchitect `*_cycles.txt`.
#' @param cycle Which cycle(s) to use: `"dominant"` (default, highest copy count)
#'   or `"all"`.
#' @return A [data.table::data.table] of `chr`, `start`, `end` (exclusive end), or
#'   `NULL` if no usable cycle/segments are found.
#' @seealso [call_episomal_from_ampliconsuite()], [read_ampliconsuite()]
#' @export
#' @importFrom data.table data.table tstrsplit
aa_cycle_intervals <- function(cycles_file, cycle = c("dominant", "all")) {
  cycle <- match.arg(cycle)
  if (is.null(cycles_file) || length(cycles_file) != 1L || is.na(cycles_file) ||
      !file.exists(cycles_file)) return(NULL)
  L <- readLines(cycles_file, warn = FALSE)
  segl <- L[startsWith(L, "Segment")]
  cyl  <- L[startsWith(L, "Cycle=")]
  if (!length(segl) || !length(cyl)) return(NULL)
  sf  <- data.table::tstrsplit(segl, "\t", fixed = TRUE)
  seg <- data.table::data.table(id = as.character(sf[[2]]), chr = as.character(sf[[3]]),
                                start = as.numeric(sf[[4]]), end = as.numeric(sf[[5]]))
  cc     <- suppressWarnings(as.numeric(sub(".*Copy_count=([0-9.eE+-]+).*", "\\1", cyl)))
  segstr <- sub(".*Segments=", "", cyl)
  keep   <- if (cycle == "dominant") which.max(cc) else seq_along(cyl)
  ids    <- unique(unlist(lapply(segstr[keep],
                    function(s) gsub("[+-]", "", strsplit(s, ",", fixed = TRUE)[[1]]))))
  ids    <- setdiff(ids, "0")                       # segment 0 = AA source vertex
  iv <- seg[id %in% ids & !is.na(start) & !is.na(end), .(chr, start, end = end + 1)]
  if (!nrow(iv)) return(NULL)
  unique(iv)
}

#' Call episomal ecDNA from AmpliconSuite output
#'
#' Classifies a single AmpliconClassifier **ecDNA feature** as episomal
#' (simple-excision circle) or not, from its AmpliconArchitect graph and its
#' ecDNA interval(s). The amplicon extent is taken from the AmpliconClassifier
#' ecDNA interval BED (AA's own amplicon call) -- no copy-number threshold is
#' applied here. The feature is **episomal** when a BOUNDARY DUP exists: an
#' intra-chromosomal DUP whose ends bracket the outer boundary of the ecDNA
#' interval on its chromosome (`lo <= interval start` and `hi >= interval end`).
#' The boundary DUP reported is the one with the highest read-pair support among
#' spanning DUPs. A multi-chromosome ecDNA is **chimeric** episomal when every
#' chromosome is bracketed by a boundary DUP and joined by >= 1 translocation.
#'
#' @param graph_file Path to the amplicon's AmpliconArchitect `*_graph.txt`.
#' @param intervals The ecDNA feature's intervals: a path to an AmpliconClassifier
#'   `*_ecDNA_<k>_intervals.bed`, or a `data.frame`/`data.table`/GRanges-like
#'   object with `chr`/`start`/`end`. AmpliconClassifier's exclusive-end
#'   convention is handled internally.
#' @param tol Integer; coordinate tolerance (bp) for the boundary-DUP bracket
#'   test (default `0`; graph and interval coordinates come from the same graph
#'   so no tolerance is normally needed).
#' @param min_occupancy Numeric in `(0, 1]`; minimum fraction of the boundary-DUP
#'   bracketed span that must be covered by the ecDNA intervals for an episomal
#'   call (default `0.05`, i.e. reject when >95% of the bracketed span is
#'   non-amplified). A simple-excision circle is near-contiguous, so a DUP
#'   spanning a largely non-amplified region (e.g. a lone DUP across a diploid
#'   chromosome arm) is rejected as not a self-contained circle. The default was
#'   set by a PCAWG-vs-PURPLE concordance sweep (optimal at 2-5%; 5% additionally
#'   removes scattered low-occupancy secondary features).
#' @param trim_to_amplified Logical; trim the AmpliconClassifier ecDNA interval to
#'   its amplified extent before the boundary-DUP test (default `TRUE`).
#'   AmpliconClassifier occasionally pads an interval with flanking diploid
#'   sequence, pushing the outer extent past the boundary DUP and mis-calling a
#'   clean episome as complex; trimming uses the AA graph's own segment copy
#'   numbers to drop leading/trailing near-diploid segments (edges only -- internal
#'   low-CN segments are kept). It can only shrink the extent, so it may rescue a
#'   false-complex call but never breaks a true-episomal one.
#' @param min_seg_cn_frac Numeric in `(0, 1]`; a segment counts as amplified for
#'   `trim_to_amplified` when its copy number is at least this fraction of the
#'   amplicon's peak segment copy number (default `0.10`). Relative to the
#'   amplicon's own amplification, not a global ploidy cutoff.
#' @param structural_qc Logical (default `TRUE`); apply structural quality control
#'   to a bracketed amplicon before calling it episomal, using the AA graph's own
#'   segment copy numbers, source (open-end) edges and junction copy numbers.
#'   Rejects (→ `status = "complex"`) an amplicon that is not amplified at both
#'   boundary-DUP breakends, is OPEN at the boundary (a large source edge at a
#'   breakend), is closed by a boundary DUP whose copy number is far below the
#'   circle CN, or has an INV/TRA at a boundary out-copying the DUP. Ported from
#'   the PURPLE-native [call_simple_excision()]. Set `FALSE` for the older
#'   bracket-and-occupancy behaviour.
#' @param qc_amp_frac,qc_src_frac,qc_jcn_frac Numeric in `(0, 1]`; the structural-QC
#'   thresholds, each relative to the width-weighted median copy number across the
#'   boundary-DUP span (the "circle-body" CN). A boundary breakend must be amplified
#'   to `>= qc_amp_frac` × circle CN (default `0.5`); a boundary source edge
#'   `>= qc_src_frac` × circle CN marks an open boundary (default `0.5`); the
#'   boundary DUP's junction copy number must be `>= qc_jcn_frac` × circle CN
#'   (default `0.1`; a low floor that rejects only a boundary DUP explaining a
#'   tiny fraction of the amplification, while tolerating internally re-amplified
#'   episomes whose resident locus is boosted above the circularising junction's
#'   copy number — e.g. DO11854 EGFR).
#' @param qc_edge_tol Integer; window (bp) around a boundary-DUP breakend within
#'   which a source edge or INV/TRA breakend counts as "at the boundary"
#'   (default `1e4`).
#' @param ploidy,cnv_file Tumour ploidy for the amplification / tandem-duplication
#'   gate. Supply `ploidy` directly, or `cnv_file` -- a path to a genome-wide
#'   copy-number BED (e.g. the AA-input `*_CNV_CALLS.bed`) whose width-weighted
#'   median CN is used. When ploidy is known and structural QC is on, a bracketed
#'   amplicon whose circle-body CN does not exceed `min_cn_ratio` x ploidy is
#'   downgraded (`qc_fail = "below_ploidy"`): at low copy number a boundary DUP +
#'   copy-number step is indistinguishable from an in-place tandem duplication by
#'   unphased short-read data, so a genuine ecDNA must be amplified well above
#'   ploidy. Both `NULL` (the default) skips the gate.
#' @param min_cn_ratio Numeric; the circle-body CN must exceed this multiple of
#'   ploidy to pass the amplification gate (default `3`, matching the
#'   PURPLE-native [call_simple_excision()]).
#'
#' A **bridged/ambiguous downgrade** applies even when a bracketing boundary DUP
#' exists: if the amplicon's highest read-pair-support junction (a breakend inside
#' the ecDNA extent) has its other breakend on the **same chromosome**, > `bridge_gap`
#' beyond the extent, and in an amplified segment, the dominant rearrangement
#' co-mobilises a co-amplified region, so `status` is set to `"bridged"` and
#' `episomal` to `FALSE`. Cross-chromosome TRAs are the chimeric-episome signature
#' and are NOT downgraded.
#'
#' @return A one-row [data.table::data.table] with `episomal` (logical; `TRUE`
#'   only for a clean episome), `status` (`"episomal"`/`"bridged"`/`"complex"`),
#'   `chimeric` (logical), `n_chr` (chromosomes spanned by the ecDNA), the
#'   boundary-DUP record `bd_chr`, `bd_lo`, `bd_hi`, `bd_pe_support`, `bd_jcn`,
#'   `bd_homlen` (`NA` when no boundary DUP), `bd_occupancy` (ecDNA-interval
#'   coverage of the boundary-DUP span; the occupancy-guard metric), the
#'   structural-QC record `bd_circle_cn` (width-weighted median CN of the
#'   boundary-DUP span), `bd_jcn_frac` (boundary-DUP jcn ÷ circle CN), `bd_src_frac`
#'   (boundary source-edge CN ÷ circle CN), `ploidy` (resolved tumour ploidy or
#'   `NA`), `cn_over_ploidy` (circle CN ÷ ploidy), `qc_pass` (logical) and
#'   `qc_fail` (`;`-separated failed checks incl. `below_ploidy`, `NA` if passed
#'   or QC off), the
#'   bridging-junction record
#'   `bridge_pe_support`, `bridge_partner` (`chr:pos` of the out-of-amplicon
#'   breakend; `NA` unless downgraded), and `n_junctions`.
#' @seealso [read_ampliconsuite()], [call_simple_excision()]
#' @examples
#' \dontrun{
#' call_episomal_from_ampliconsuite(
#'   graph_file = "SAMPLE_amplicon1_graph.txt",
#'   intervals  = "SAMPLE_amplicon1_ecDNA_1_intervals.bed")
#' }
#' @export
#' @importFrom data.table data.table rbindlist
call_episomal_from_ampliconsuite <- function(graph_file, intervals, tol = 0L,
                                             bridge_gap = 1e6, min_occupancy = 0.05,
                                             trim_to_amplified = TRUE, min_seg_cn_frac = 0.10,
                                             structural_qc = TRUE, qc_amp_frac = 0.5,
                                             qc_src_frac = 0.5, qc_jcn_frac = 0.1,
                                             qc_edge_tol = 1e4, ploidy = NULL,
                                             cnv_file = NULL, min_cn_ratio = 3) {
  none <- data.table::data.table(
    episomal = NA, status = NA_character_, chimeric = NA, n_chr = NA_integer_,
    bd_chr = NA_character_, bd_lo = NA_real_, bd_hi = NA_real_, bd_pe_support = NA_real_,
    bd_jcn = NA_real_, bd_homlen = NA_real_, bd_occupancy = NA_real_,
    bd_circle_cn = NA_real_, bd_jcn_frac = NA_real_, bd_src_frac = NA_real_,
    ploidy = NA_real_, cn_over_ploidy = NA_real_,
    qc_pass = NA, qc_fail = NA_character_,
    bridge_pe_support = NA_real_, bridge_partner = NA_character_,
    n_junctions = NA_integer_)

  ## accept a bare AA cycles file as the interval source (generalises to any AA output)
  if (is.character(intervals) && length(intervals) == 1L && !is.na(intervals) &&
      grepl("_cycles\\.txt$", intervals))
    intervals <- aa_cycle_intervals(intervals)
  if (is.null(intervals) ||
      (is.character(intervals) &&
       (length(intervals) != 1L || is.na(intervals) || !file.exists(intervals))))
    return(none)
  iv <- .read_ac_intervals(intervals)
  if (!nrow(iv)) return(none)
  gr <- read_ampliconsuite(graph_file)
  j  <- gr$junctions

  ## Trim the AC interval to its AMPLIFIED extent (default on). AmpliconClassifier
  ## sometimes pads an ecDNA interval with flanking DIPLOID sequence -- e.g. GBM39EC:
  ## ~100 kb of CN~4 either side of a CN~91 EGFR circle -- which pushes the outer
  ## extent past the boundary DUP and mis-calls a clean episome as complex. Using the
  ## AA graph's own segment copy numbers, drop leading/trailing segments whose CN is
  ## < `min_seg_cn_frac` of the amplicon's peak CN, per chromosome. Only the EDGES are
  ## trimmed: internal low-CN segments (excision scars, chromothriptic deletions) sit
  ## between amplified segments and are retained. Shrinking the extent only makes the
  ## boundary-DUP bracket easier, so this can rescue a false-complex but never break a
  ## true-episomal call.
  if (isTRUE(trim_to_amplified) && !is.null(gr$segments) && nrow(gr$segments)) {
    seg <- gr$segments
    bnd <- iv[, {
      s   <- seg[seg$chr == .BY$chr & seg$end >= min(start) & seg$start <= max(end), ]
      ## Robust peak CN: ignore sub-2kb segments so a tiny high-CN spike cannot set
      ## the reference and trim away a genuine amplicon body.
      ref  <- s[is.finite(s$cn) & (s$end - s$start) >= 2000, ]
      peak <- if (nrow(ref)) max(ref$cn) else if (nrow(s)) max(s$cn, na.rm = TRUE) else NA_real_
      amp  <- if (nrow(s) && is.finite(peak)) s[is.finite(s$cn) & s$cn >= min_seg_cn_frac * peak, ] else s
      ## clamp the amplified span to WITHIN the original interval so trimming can
      ## only shrink the extent (a segment straddling the edge must not expand it)
      if (nrow(amp)) .(lo2 = max(min(start), min(amp$start)), hi2 = min(max(end), max(amp$end)))
      else           .(lo2 = min(start),     hi2 = max(end))
    }, by = chr]
    iv2 <- merge(iv, bnd, by = "chr", all.x = TRUE)
    iv2 <- iv2[end >= lo2 & start <= hi2]
    if (nrow(iv2)) {                              # keep original interval if trimming would empty it
      iv2[, start := pmax(start, lo2)][, end := pmin(end, hi2)][, c("lo2", "hi2") := NULL]
      iv <- iv2
    }
  }

  ## per-chromosome outer extent of the ecDNA feature
  ext <- iv[, .(lo = min(start), hi = max(end)), by = chr]
  chrs <- ext$chr

  ## spanning (boundary) DUPs: intra-chromosomal DUP bracketing the extent on chr
  span <- data.table::rbindlist(lapply(chrs, function(cc) {
    e <- ext[chr == cc]
    if (is.null(j) || !nrow(j)) return(NULL)
    j[svclass == "DUP" & cA == cc & cB == cc & lo <= e$lo + tol & hi >= e$hi - tol]
  }), fill = TRUE)
  has_span <- vapply(chrs, function(cc)
    !is.null(span) && nrow(span) > 0 && nrow(span[cA == cc]) > 0, logical(1))

  if (length(chrs) == 1L) {
    episomal <- isTRUE(has_span[[1]]); chimeric <- FALSE
  } else {
    ntra <- if (is.null(j) || !nrow(j)) 0L else nrow(j[svclass == "TRA"])
    episomal <- all(has_span) && ntra >= 1L
    chimeric <- episomal
  }

  ## Occupancy of the selected boundary DUP on each bracketed chromosome:
  ## (summed ecDNA-interval width on the chr) / (boundary-DUP-bracketed span).
  occ_of <- function(cc) {
    sp <- if (is.null(span) || !nrow(span)) NULL else span[cA == cc]
    if (is.null(sp) || !nrow(sp)) return(NA_real_)
    d <- sp[which.max(pe_support)]; dspan <- d$hi - d$lo
    if (is.na(dspan) || dspan <= 0) return(NA_real_)
    sum(iv[chr == cc, end - start]) / dspan
  }
  occ <- vapply(chrs, occ_of, numeric(1))
  bd_occupancy <- suppressWarnings(min(occ, na.rm = TRUE))
  if (!is.finite(bd_occupancy)) bd_occupancy <- NA_real_

  ## OCCUPANCY GUARD. A simple-excision circle is near-contiguous. Reject an
  ## episomal call whose boundary DUP brackets a span that is mostly NON-amplified
  ## -- a lone DUP spanning a largely diploid region (e.g. a whole chromosome arm)
  ## is not a self-contained circle. Every bracketed chromosome's ecDNA-interval
  ## coverage of its boundary-DUP span must reach `min_occupancy`.
  if (isTRUE(episomal) && any(!is.na(occ) & occ < min_occupancy)) {
    episomal <- FALSE; chimeric <- FALSE
  }

  bd <- if (episomal && !is.null(span) && nrow(span)) span[which.max(pe_support)] else NULL

  ## Bridged/ambiguous downgrade. Even with a valid bracketing boundary DUP, the
  ## amplicon is not a self-contained circle if its BEST-SUPPORTED junction reaches
  ## OUT of the ecDNA extent into another amplified region: the dominant
  ## rearrangement co-mobilises sequence beyond the amplicon. Rank amplicon
  ## junctions (>= 1 breakend inside the extent of one of the feature's
  ## chromosomes) by read-pair support (pe_support, NOT jcn); if the top one has
  ## its OTHER breakend outside every feature interval and that breakend sits in an
  ## amplified segment, downgrade `episomal` to the "bridged" status.
  status <- if (isTRUE(episomal)) "episomal" else "complex"

  ## STRUCTURAL QC (default on). A bracketing boundary DUP is necessary but not
  ## sufficient for a clean simple-excision circle. Using the AA graph's own
  ## segment copy numbers, source (open-end) edges and junction copy numbers,
  ## reject a bracketed amplicon that is:
  ##   (a) not amplified at BOTH boundary-DUP breakends (CN < qc_amp_frac x the
  ##       circle-body CN) -- a DUP reaching into non-amplified sequence;
  ##   (b) OPEN at the boundary -- a source edge >= qc_src_frac x circle CN within
  ##       qc_edge_tol of a breakend (copy number leaks into the genome there, so
  ##       the boundary is not a closed circularisation). An INTERNAL open
  ##       co-amplification (away from the boundary) does NOT disqualify;
  ##   (c) closed by a boundary DUP whose junction copy number is << the circle CN
  ##       (jcn < qc_jcn_frac x circle CN) -- the DUP explains only a fraction of
  ##       the amplification, i.e. it is a passenger/partial tandem duplication;
  ##   (d) has an INV/TRA at a boundary breakend that out-copies the DUP -- the
  ##       edge is a fold-back / translocation junction, not a circularisation.
  ## The circle-body CN is the width-weighted median segment CN across the
  ## boundary-DUP span (robust to internal high-CN spikes). Ported from the
  ## PURPLE-native call_simple_excision(); calibrated on the GBM39 gold
  ## (jcn/circle 0.98, src/circle 0.11) and validated to retain 13/13 primary-GBM
  ## EGFR episomes while rejecting passenger, open and chromothriptic amplicons.
  ## Resolve tumour ploidy for the amplification (tandem-dup) gate: a direct
  ## `ploidy` value takes precedence, else the genome-wide median of `cnv_file`.
  ploidy_val <- if (is.numeric(ploidy) && length(ploidy) == 1L && is.finite(ploidy)) ploidy
                else .cnv_ploidy(cnv_file)

  bd_circle_cn <- NA_real_; bd_jcn_frac <- NA_real_; bd_src_frac <- NA_real_
  cn_over_ploidy <- NA_real_
  qc_pass <- NA; qc_fail <- NA_character_
  if (isTRUE(structural_qc) && !is.null(bd)) {
    seg <- gr$segments; bc <- bd$cA
    scc <- if (!is.null(seg) && nrow(seg)) seg[chr == bc & end >= bd$lo & start <= bd$hi] else seg
    ccn <- if (!is.null(scc) && nrow(scc)) {
      w <- pmax(1, pmin(scc$end, bd$hi) - pmax(scc$start, bd$lo)); o <- order(scc$cn)
      scc$cn[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]]
    } else NA_real_
    cn_at <- function(p) {
      s <- if (!is.null(seg) && nrow(seg)) seg[chr == bc & end >= p - 300 & start <= p + 300] else NULL
      if (!is.null(s) && nrow(s)) max(s$cn, na.rm = TRUE) else NA_real_
    }
    cn_lo <- cn_at(bd$lo); cn_hi <- cn_at(bd$hi)
    sr <- gr$source
    src_bnd <- if (!is.null(sr) && nrow(sr))
      suppressWarnings(max(c(0, sr[chr == bc &
        (abs(pos - bd$lo) <= qc_edge_tol | abs(pos - bd$hi) <= qc_edge_tol)]$cn), na.rm = TRUE)) else 0
    ej <- if (!is.null(j) && nrow(j)) j[cA == bc & cB == bc &
        svclass %in% c("h2hINV", "t2tINV", "INV", "TRA") &
        (abs(lo - bd$lo) <= qc_edge_tol | abs(hi - bd$hi) <= qc_edge_tol |
         abs(lo - bd$hi) <= qc_edge_tol | abs(hi - bd$lo) <= qc_edge_tol)] else NULL
    cn_ok  <- is.finite(ccn) && ccn > 0
    a_ok <- cn_ok && is.finite(cn_lo) && is.finite(cn_hi) &&
            cn_lo >= qc_amp_frac * ccn && cn_hi >= qc_amp_frac * ccn
    b_ok <- !cn_ok || (src_bnd / ccn) < qc_src_frac
    c_ok <- is.finite(bd$jcn) && cn_ok && (bd$jcn / ccn) >= qc_jcn_frac
    d_ok <- is.null(ej) || !nrow(ej) || all(ej$jcn <= (if (is.finite(bd$jcn)) bd$jcn else 0))
    ## (e) Amplification / tandem-duplication gate: the circle-body CN must exceed
    ## min_cn_ratio x tumour ploidy. At low copy number a boundary DUP + CN step is
    ## indistinguishable from an in-place tandem duplication by unphased short reads,
    ## so require the amplification to be well above ploidy. Skipped (always passes)
    ## when ploidy is unknown -- supply `ploidy` or `cnv_file` to enable it.
    e_ok <- !(is.finite(ploidy_val) && ploidy_val > 0 && cn_ok) ||
            (ccn > min_cn_ratio * ploidy_val)
    bd_circle_cn <- ccn
    bd_jcn_frac  <- if (is.finite(bd$jcn) && cn_ok) bd$jcn / ccn else NA_real_
    bd_src_frac  <- if (cn_ok) src_bnd / ccn else NA_real_
    cn_over_ploidy <- if (is.finite(ploidy_val) && ploidy_val > 0 && cn_ok) ccn / ploidy_val else NA_real_
    qc_pass <- isTRUE(a_ok && b_ok && c_ok && d_ok && e_ok)
    qc_fail <- paste(c("both_ends_amplified", "open_boundary", "jcn_below_cn",
                       "invtra_boundary", "below_ploidy")[
      c(!a_ok, !b_ok, !c_ok, !d_ok, !e_ok)], collapse = ";")
    if (!qc_pass) { episomal <- FALSE; chimeric <- FALSE; status <- "complex" }
  }

  bridge_pe <- NA_real_; bridge_partner <- NA_character_
  if (isTRUE(episomal) && !is.null(j) && nrow(j)) {
    in_ext <- function(ch, pos)
      ch %in% ext$chr && any(ext$chr == ch & pos >= ext$lo - tol & pos <= ext$hi + tol)
    amp_j <- j[mapply(in_ext, cA, pA) | mapply(in_ext, cB, pB)]
    if (nrow(amp_j)) {
      top <- amp_j[which.max(pe_support)]
      a_in <- in_ext(top$cA, top$pA); b_in <- in_ext(top$cB, top$pB)
      if (xor(a_in, b_in)) {                                   # exactly one end inside
        in_c  <- if (a_in) top$cA else top$cB
        far_c <- if (a_in) top$cB else top$cA
        far_p <- if (a_in) top$pB else top$pA
        e <- ext[chr == in_c]
        ## SAME-CHROMOSOME bridge only: far end on the amplicon's own chromosome but
        ## > bridge_gap beyond the extent. Cross-chromosome TRAs are the
        ## chimeric-episome signature and are deliberately NOT downgraded.
        far_region <- far_c == in_c && nrow(e) > 0 &&
          (far_p < e$lo - bridge_gap || far_p > e$hi + bridge_gap)
        seg <- gr$segments
        far_amp <- if (is.null(seg) || !nrow(seg)) TRUE else {
          cov <- seg[chr == far_c & start <= far_p & end >= far_p]
          nrow(cov) == 0 || any(cov$cn >= 0.5 * max(seg$cn, na.rm = TRUE))
        }
        if (far_region && far_amp) {
          status <- "bridged"; episomal <- FALSE
          bridge_pe <- top$pe_support
          bridge_partner <- paste0(far_c, ":", format(far_p, scientific = FALSE, trim = TRUE))
        }
      }
    }
  }

  data.table::data.table(
    episomal = episomal, status = status, chimeric = chimeric, n_chr = length(chrs),
    bd_chr = if (!is.null(bd)) bd$cA else NA_character_,
    bd_lo  = if (!is.null(bd)) bd$lo else NA_real_,
    bd_hi  = if (!is.null(bd)) bd$hi else NA_real_,
    bd_pe_support = if (!is.null(bd)) bd$pe_support else NA_real_,
    bd_jcn = if (!is.null(bd)) bd$jcn else NA_real_,
    bd_homlen = if (!is.null(bd)) bd$homlen else NA_real_,
    bd_occupancy = bd_occupancy,
    bd_circle_cn = bd_circle_cn, bd_jcn_frac = bd_jcn_frac, bd_src_frac = bd_src_frac,
    ploidy = ploidy_val, cn_over_ploidy = cn_over_ploidy,
    qc_pass = qc_pass,
    qc_fail = if (!is.na(qc_fail) && nzchar(qc_fail)) qc_fail else NA_character_,
    bridge_pe_support = bridge_pe, bridge_partner = bridge_partner,
    n_junctions = if (is.null(j)) 0L else nrow(j))
}
