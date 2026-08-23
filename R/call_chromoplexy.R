## ---------------------------------------------------------------------------
## Chromoplexy: closed chains of balanced, low-copy rearrangements whose
## breakends meet in "deletion bridges" (Baca et al., Cell 2013; ChainFinder).
##
## Unlike the amplicon-formation callers (call_simple_excision / call_brf /
## call_bfb / call_micronucleation / call_translocation_bridge_amp), chromoplexy
## is NOT amplicon-centric: it is a sample-level event over the whole genome, and
## it is defined by LOW copy number (a balanced rearrangement, not an
## amplification). It therefore runs on its own sample-level engine rather than
## the per-amplicon annotate_amplicon() engine.
##
## Model (closed cycles only): reconstruct each SV as a junction between two
## breakends; keep balanced (low-copy) junctions; collapse breakends that sit
## within `max_dist` of each other into "loci" (a locus holding ends from two
## different junctions is a deletion bridge); build a graph whose nodes are loci
## and whose edges are the long-range junctions; and call a chromoplexy for every
## connected component that is a CLEAN SIMPLE CYCLE -- every locus of degree two,
## with at least `min_pairs` (default 3, per Baca et al.) bridges closing the
## ring (chr_a -> chr_b -> chr_c -> chr_a). Open chains and higher-degree tangles
## are reported (`topology`) but not called.
## ---------------------------------------------------------------------------

## Single-linkage cluster of breakend positions into loci: walking positions in
## (chromosome, coordinate) order, a gap greater than `max_dist` (or a change of
## chromosome) starts a new locus. Returns an integer locus id per input element.
.assign_loci <- function(chr, pos, max_dist) {
  n <- length(pos)
  if (!n) return(integer(0))
  o <- order(chr, pos)
  locus <- integer(n)
  cur <- 0L; prev_chr <- NA_character_; prev_pos <- NA_real_
  for (i in o) {
    if (is.na(prev_chr) || chr[i] != prev_chr || (pos[i] - prev_pos) > max_dist)
      cur <- cur + 1L
    locus[i] <- cur
    prev_chr <- chr[i]; prev_pos <- pos[i]
  }
  locus
}

## Connected components by union-find over edges (node ids must be 1..N). Returns
## the component root of every node. No igraph dependency.
.uf_components <- function(a, b, N) {
  parent <- seq_len(N)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (k in seq_along(a)) {
    ra <- find(a[k]); rb <- find(b[k])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(N), find, integer(1))
}

## Reconstruct one junction (row) per `event` from the per-breakend table, keeping
## only events with exactly two breakends. `cn`/`ploidy` are taken at the higher
## copy-number end (a junction counts as amplified if EITHER end is amplified);
## `span` is the intrachromosomal distance, or Inf for an interchromosomal join.
.junctions_from_breakpoints <- function(bp) {
  data.table::setorder(bp, event)
  jl <- bp[, {
    if (.N != 2L) NULL else {
      cnv <- PURPLE_CN; cnv[is.na(cnv)] <- 0
      hi <- which.max(cnv)
      list(chr1 = seqnames[1], pos1 = start[1],
           chr2 = seqnames[2], pos2 = start[2],
           svclass = svclass[1],
           cn = cnv[hi], ploidy = ploidy[hi],
           interchrom = seqnames[1] != seqnames[2],
           span = if (seqnames[1] == seqnames[2]) as.numeric(abs(start[2] - start[1])) else Inf)
    }
  }, by = event]
  jl
}

## Attach the local per-chromosome ploidy baseline to each breakend by overlap
## with the sample copy-number segments (default 2 where no segment overlaps).
.attach_breakend_ploidy <- function(bp, cnv_dt) {
  bp$ploidy <- 2
  if (!nrow(cnv_dt) || !"ploidy" %in% names(cnv_dt)) return(bp)
  bgr <- GenomicRanges::GRanges(bp$seqnames, IRanges::IRanges(bp$start, width = 1L))
  cgr <- GenomicRanges::GRanges(cnv_dt$seqnames, IRanges::IRanges(cnv_dt$start, cnv_dt$end))
  ov <- GenomicRanges::findOverlaps(bgr, cgr, select = "first")
  hit <- !is.na(ov)
  bp$ploidy[hit] <- as.numeric(cnv_dt$ploidy[ov[hit]])
  bp
}

## Per-sample chromoplexy detector: everything from a breakend table + CN table
## to the event rows for one sample. Returns a data.table (possibly 0 rows).
.detect_chromoplexy <- function(bp, cnv_dt, sample_id, min_pairs, min_span,
                                max_dist, max_cn, footprint_width, max_small,
                                min_chromosomes) {
  empty <- .chromoplexy_schema()
  if (!nrow(bp)) return(empty)
  bp <- .attach_breakend_ploidy(bp, cnv_dt)
  jl <- .junctions_from_breakpoints(bp)
  if (!nrow(jl)) return(empty)

  ## (1) drop small local dup/del noise
  small <- jl$svclass %in% c("DEL", "DUP") & jl$span <= max_small
  jl <- jl[!small]
  if (!nrow(jl)) return(empty)

  ## (2) keep only BALANCED (low-copy) junctions: ploidy-scaled, like every other
  ## EpiTracer caller (robust on polysomic chromosomes, e.g. GBM chr7).
  jl <- jl[cn <= max_cn * ploidy / 2]
  if (nrow(jl) < min_pairs) return(empty)
  jl[, jid := .I]

  ## (3) collapse the 2*nrow(jl) breakends into loci (deletion-bridge candidates)
  ends <- data.table::data.table(
    jid = rep(jl$jid, 2L),
    end = rep(c(1L, 2L), each = nrow(jl)),
    chr = c(jl$chr1, jl$chr2),
    pos = c(jl$pos1, jl$pos2))
  ends[, locus := .assign_loci(chr, pos, max_dist)]

  ## edge = one junction connecting the loci of its two ends; drop self-loops
  el <- data.table::dcast(ends, jid ~ end, value.var = "locus")
  data.table::setnames(el, c("1", "2"), c("l1", "l2"))
  el <- merge(el, jl, by = "jid")
  el <- el[l1 != l2]
  ## (4) chain edges must be long-range: interchromosomal, or a long intrachrom span
  el <- el[interchrom | span >= min_span]
  if (nrow(el) < min_pairs) return(empty)

  ## relabel participating loci to 1..N and compute degree over the kept edges
  used <- sort(unique(c(el$l1, el$l2)))
  el[, `:=`(n1 = match(l1, used), n2 = match(l2, used))]
  N <- length(used)
  deg <- tabulate(c(el$n1, el$n2), nbins = N)
  comp <- .uf_components(el$n1, el$n2, N)
  el[, comp := comp[n1]]

  ## (5) a CLEAN SIMPLE CYCLE = connected component whose every locus has degree 2
  ## (a connected 2-regular graph is a single cycle) and that has >= min_pairs
  ## bridges. Report the shape of every component; call only the clean cycles.
  ends_gr_all <- GenomicRanges::GRanges(ends$chr, IRanges::IRanges(ends$pos, width = 1L),
                                        jid = ends$jid)
  rows <- vector("list", 0L); cid <- 0L
  for (cc in unique(el$comp)) {
    ce <- el[comp == cc]
    nodes <- unique(c(ce$n1, ce$n2))
    clean_cycle <- all(deg[nodes] == 2L) && nrow(ce) == length(nodes)
    topology <- if (clean_cycle) "cycle" else if (any(deg[nodes] > 2L)) "complex" else "open"
    if (!(clean_cycle && nrow(ce) >= min_pairs)) next
    ## chromoplexy must genuinely bridge chromosomes, not close a ring on one:
    ## require the cycle to span at least `min_chromosomes` distinct chromosomes.
    chrs <- sort(unique(c(ce$chr1, ce$chr2)))
    if (length(chrs) < min_chromosomes) next
    cid <- cid + 1L
    ## loci footprints for this cycle (min-max breakend position within each locus)
    cyc_loci <- used[nodes]
    lf <- ends[locus %in% cyc_loci, .(lo = min(pos), hi = max(pos)), by = .(chr, locus)]
    foot <- paste(sprintf("%s:%d-%d", lf$chr, lf$lo, lf$hi), collapse = ";")
    ## intruding junctions: breakends of OTHER junctions landing in the cycle loci
    loci_gr <- GenomicRanges::GRanges(lf$chr, IRanges::IRanges(pmax(1, lf$lo - max_dist), lf$hi + max_dist))
    other <- ends_gr_all[!(ends_gr_all$jid %in% ce$jid)]
    num_other <- length(S4Vectors::queryHits(GenomicRanges::findOverlaps(other, loci_gr)))
    rows[[cid]] <- data.table::data.table(
      WGS_ID = sample_id, chromoplexy_id = cid, topology = topology,
      n_junctions = nrow(ce), n_bridges = length(nodes),
      n_chromosomes = length(chrs), chromosomes = paste(chrs, collapse = ","),
      max_cn = max(ce$cn), min_cn = min(ce$cn),
      min_span = min(ce$span), n_interchrom = sum(ce$interchrom),
      num_other = num_other,
      frac_cp = nrow(ce) / (nrow(ce) + num_other),
      footprint = foot, events = paste(ce$event, collapse = ","))
  }
  if (!length(rows)) return(empty)
  data.table::rbindlist(rows)
}

## Typed 0-row result (the callers' "always return a typed table" contract).
.chromoplexy_schema <- function() {
  data.table::data.table(
    WGS_ID = character(), chromoplexy_id = integer(), topology = character(),
    n_junctions = integer(), n_bridges = integer(), n_chromosomes = integer(),
    chromosomes = character(), max_cn = numeric(), min_cn = numeric(),
    min_span = numeric(), n_interchrom = integer(), num_other = integer(),
    frac_cp = numeric(), footprint = character(), events = character())
}

#' Call chromoplexy (closed balanced rearrangement chains)
#'
#' Sample-level caller for chromoplexy, the punctuated formation of a closed chain
#' of balanced, low-copy rearrangements whose breakends meet in *deletion bridges*
#' (Baca et al., *Cell* 2013, doi:10.1016/j.cell.2013.03.021; graph criteria after
#' the gGnome `chromoplexy()` caller, Hadi et al., *Cell* 2020,
#' doi:10.1016/j.cell.2020.08.006). Each structural variant is reconstructed as a
#' junction between two breakends; balanced (ploidy-scaled low-copy) junctions are
#' kept; breakends within `max_dist` collapse into loci (a locus joining two
#' distinct junctions is a deletion bridge); and a chromoplexy is called for every
#' connected component of the locus/junction graph that forms a **clean simple
#' cycle** -- every locus of degree two, with at least `min_pairs` bridges closing
#' the ring. Runs independently of the amplicon-formation callers.
#'
#' @param breakpoints_gr A [GenomicRanges::GRanges] of SV breakends (two per
#'   `event`), as used by [call_simple_excision()]: metadata `WGS_ID`, `event`,
#'   `svclass`, `PURPLE_CN`.
#' @param cnv_gr A [GenomicRanges::GRanges] of allele-specific copy-number
#'   segments with metadata `sample`, `copyNumber`, `ploidy` (the per-chromosome
#'   baseline used to ploidy-scale the low-copy filter).
#' @param min_pairs Integer; minimum number of deletion bridges (equivalently the
#'   cycle length) required to call a chromoplexy (default `3L`, per Baca et al.).
#' @param min_chromosomes Integer; minimum number of DISTINCT chromosomes the
#'   cycle must span (default `3L`). Chromoplexy is characteristically
#'   interchromosomal (Baca et al.); this rejects closed rings that merely link
#'   loci on one or two chromosomes. Set to `2L` (or `1L`) to relax.
#' @param min_span Numeric; minimum span (bp) for an intrachromosomal junction to
#'   count as a long-range chain edge (default `1e7`); interchromosomal junctions
#'   always qualify.
#' @param max_dist Integer; how close (bp) two breakends must sit to collapse into
#'   the same locus / deletion bridge (default `1e4`).
#' @param max_cn Numeric; a junction is balanced (kept) when its copy number is
#'   at most `max_cn * ploidy / 2` (default `3`, i.e. <=3 copies at diploid
#'   baseline, scaled up on polysomic chromosomes).
#' @param footprint_width Numeric; padding (bp) applied when reporting locus
#'   footprints (default `1e6`).
#' @param max_small Numeric; a DUP/DEL junction spanning at most this (bp) is
#'   dropped as local noise before chain-building (default `5e4`).
#' @param mc.cores Integer; samples processed in parallel (default `1`).
#' @return A [data.table::data.table] with one row per called chromoplasy cycle:
#'   `WGS_ID`, `chromoplexy_id`, `topology` (always `"cycle"` for called events),
#'   `n_junctions`, `n_bridges`, `n_chromosomes`, `chromosomes`, `max_cn`,
#'   `min_cn`, `min_span`, `n_interchrom`, `num_other` (intruding foreign
#'   breakends), `frac_cp` (cleanliness `n_junctions / (n_junctions + num_other)`),
#'   `footprint` and `events` (comma-separated `event` ids in the cycle, for
#'   joining back onto amplicons). Empty (0 rows) when no chromoplexy is found.
#' @seealso [call_simple_excision()], [call_translocation_bridge_amp()]
#' @export
call_chromoplexy <- function(breakpoints_gr, cnv_gr, min_pairs = 3L, min_span = 1e7,
                             max_dist = 1e4, max_cn = 3, footprint_width = 1e6,
                             max_small = 5e4, min_chromosomes = 3L, mc.cores = 1) {
  bp_all <- gr2dt(breakpoints_gr)
  cn_all <- gr2dt(cnv_gr)
  if (!nrow(bp_all)) return(.chromoplexy_schema())
  samples <- unique(bp_all$WGS_ID)
  res <- parallel::mclapply(samples, function(s) {
    bp <- bp_all[WGS_ID == s]
    cnv_dt <- cn_all[sample == s]
    .detect_chromoplexy(bp, cnv_dt, s, min_pairs, min_span, max_dist, max_cn,
                        footprint_width, max_small, min_chromosomes)
  }, mc.cores = mc.cores)
  out <- data.table::rbindlist(Filter(function(x) nrow(x) > 0L, res), fill = TRUE)
  if (!nrow(out)) return(.chromoplexy_schema())
  out
}
