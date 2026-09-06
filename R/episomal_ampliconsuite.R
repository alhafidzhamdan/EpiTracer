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
#' @return A list with two [data.table::data.table]s:
#' \describe{
#'   \item{`segments`}{`chr`, `start`, `end`, `cn` (predicted copy number).}
#'   \item{`junctions`}{`cA`,`pA`,`sA` and `cB`,`pB`,`sB` (the two breakends'
#'     chromosome/position/strand), `jcn` (junction predicted copy number),
#'     `pe_support` (number of supporting read pairs), `homlen` (breakpoint
#'     homology length, `NA` if absent), `svclass` (`DUP`/`DEL`/`INV`/`TRA`) and
#'     `lo`/`hi` (the lower/upper breakend coordinate).}
#' }
#' Either table may have zero rows.
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
  if (is.null(graph_file) || is.na(graph_file) || !file.exists(graph_file))
    return(list(segments = seg, junctions = jun))
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
  list(segments = seg, junctions = jun)
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
                                             bridge_gap = 1e6, min_occupancy = 0.05) {
  none <- data.table::data.table(
    episomal = NA, chimeric = NA, n_chr = NA_integer_, bd_chr = NA_character_,
    bd_lo = NA_real_, bd_hi = NA_real_, bd_pe_support = NA_real_,
    bd_jcn = NA_real_, bd_homlen = NA_real_, bd_occupancy = NA_real_,
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
    bridge_pe_support = bridge_pe, bridge_partner = bridge_partner,
    n_junctions = if (is.null(j)) 0L else nrow(j))
}
