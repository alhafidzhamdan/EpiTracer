## ---------------------------------------------------------------------------
## Generative simulation of episomal ecDNA amplicons: the core data structure.
##
## The design follows Rausch/Korbel-style amplicon simulators (and in particular
## the fragment-list formulation of Bernhard et al., "Chromothripsis followed by
## circular recombination drives oncogene amplification in human cancer",
## github.com/seismicon/SeismicAmplification): an amplicon is an ORDERED LIST OF
## FRAGMENTS, each a genomic interval traversed forwards or inverted. Every
## mechanism is a transformation of that list.
##
## The pay-off is that copy number and structural variants are *both derived from
## the same object*: copy number is the coverage of the fragment list, and the
## junctions are read off consecutive fragment pairs. They therefore cannot
## disagree, which is exactly what EpiTracer's callers cross-check -- a simulator
## that emits hand-written CN and SV tables (as validation/simulate_benchmark.R
## does) can never exercise that consistency.
##
## `sim_episome()` gives an amplicon its origin; `sim_micronucleation()` and
## `sim_fuse_episomes()` evolve it; `sim_to_epitracer()` renders it as the
## PURPLE-style inputs the callers take, with a sequencing noise model on top.
## ---------------------------------------------------------------------------

#' Simulated ecDNA amplicon
#'
#' The object returned by [sim_episome()] and transformed by the mechanism
#' operators ([sim_micronucleation()], [sim_fuse_episomes()], [sim_evolve()]).
#' It records the *ground truth* structure of one simulated amplicon, from which
#' [sim_to_epitracer()] renders caller inputs and [sim_to_plot_inputs()] renders
#' plotter inputs.
#'
#' @section Components:
#' \describe{
#'   \item{`sample`}{Character sample identifier.}
#'   \item{`fragments`}{A [GenomicRanges::GRanges] in **amplicon order** -- the
#'     sequence of intervals a polymerase would traverse walking the circle once.
#'     Metadata `inv` (`0` forward, `1` inverted) and `origin` (the label of the
#'     parental episome a fragment came from, which distinguishes the two halves
#'     of a chimera).}
#'   \item{`circular`}{`TRUE` for an episome; the last fragment joins the first.}
#'   \item{`copies`}{Numeric ecDNA copies per tumour cell. Copy number is
#'     `copies * coverage(fragments)` on top of the host baseline, so this sets
#'     the amplicon's copy-number level while `fragments` sets its shape.}
#'   \item{`junctions`}{A [data.table::data.table] ledger of every junction with
#'     `chrom1`, `start1`, `strand1`, `chrom2`, `start2`, `strand2`, `svclass`,
#'     `multiplicity` (traversals per circle), `round`, `mechanism`, `origin`
#'     (`"amplicon"` or `"host"`) and `homology` (the repair pathway).}
#'   \item{`host`}{A [GenomicRanges::GRanges] of host-chromosome background with
#'     `ploidy`, giving the local copy-number baseline the flanks sit at.}
#'   \item{`variants`}{`NULL` for a homogeneous amplicon. Once a single molecule
#'     is altered (see [sim_internal_deletion()]) the cell carries a MIXTURE, and
#'     this holds it as a list of `list(fragments, copies, label)`. `fragments`
#'     and `copies` above then report the most abundant species and the total.
#'     The mixture is what makes the ORDER of events observable: a junction that
#'     arose late is carried by fewer circles, so it is emitted at lower copy
#'     number and lower read support. Each species also carries a `fitness`, its
#'     replication rate relative to the founder, which is what lets a minority
#'     species expand rather than sit at a frozen share.}
#'   \item{`preserved`}{Oncogene loci that mechanisms must not delete.}
#'   \item{`history`}{Character log of the mechanism steps applied, in order.}
#'   \item{`round`}{Integer count of evolutionary rounds applied so far.}
#' }
#' @name ecdna_sim
#' @seealso [sim_episome()], [sim_micronucleation()], [sim_to_epitracer()]
NULL

## Empty junction ledger with the canonical column set / types.
.sim_empty_junctions <- function() {
  data.table::data.table(
    chrom1 = character(), start1 = numeric(), strand1 = character(),
    chrom2 = character(), start2 = numeric(), strand2 = character(),
    svclass = character(), multiplicity = numeric(),
    round = integer(), mechanism = character(),
    origin = character(), homology = character()
  )
}

## Structural-variant class from the two breakend orientations, following
## EpiTracer's convention (see read_purple_sv_vcf): +/- DEL, -/+ DUP, +/+ h2hINV,
## -/- t2tINV, inter-chromosomal TRA. Assumes breakends are already sorted so
## that (chrom1, start1) <= (chrom2, start2), as .sim_junction() guarantees.
.sim_svclass <- function(chrom1, strand1, chrom2, strand2) {
  data.table::fifelse(chrom1 != chrom2, "TRA",
    data.table::fifelse(strand1 == "+" & strand2 == "-", "DEL",
      data.table::fifelse(strand1 == "-" & strand2 == "+", "DUP",
        data.table::fifelse(strand1 == "+" & strand2 == "+", "h2hINV", "t2tINV"))))
}

## The breakend at which a junction LEAVES fragment `i` (its tail in amplicon
## order) and the breakend at which it ENTERS fragment `j` (its head). A forward
## fragment is exited at its reference right edge (orientation "+") and entered at
## its reference left edge ("-"); an inverted fragment is traversed right-to-left,
## so the two are mirrored. This single rule is what makes the emitted svclass
## mixture mechanism-specific: a circularisation gives -/+ (DUP), a fold-back
## gives +/+ (h2hINV), a reference-order rejoin gives +/- (DEL).
.sim_exit <- function(frag) {
  inv <- frag$inv[1]
  list(chrom = as.character(GenomeInfoDb::seqnames(frag))[1],
       pos    = if (inv == 0) GenomicRanges::end(frag)[1] else GenomicRanges::start(frag)[1],
       strand = if (inv == 0) "+" else "-")
}

.sim_entry <- function(frag) {
  inv <- frag$inv[1]
  list(chrom = as.character(GenomeInfoDb::seqnames(frag))[1],
       pos    = if (inv == 0) GenomicRanges::start(frag)[1] else GenomicRanges::end(frag)[1],
       strand = if (inv == 0) "-" else "+")
}

## Build one junction row from an exit breakend and an entry breakend, sorted so
## (chrom1, start1) <= (chrom2, start2) -- the same normalisation read_sv()
## applies, so simulated and real inputs are directly comparable.
.sim_junction <- function(a, b, mechanism, round, homology = "nhej",
                          origin = "amplicon", multiplicity = 1) {
  swap <- (a$chrom > b$chrom) || (a$chrom == b$chrom && a$pos > b$pos)
  if (swap) { tmp <- a; a <- b; b <- tmp }
  data.table::data.table(
    chrom1 = a$chrom, start1 = as.numeric(a$pos), strand1 = a$strand,
    chrom2 = b$chrom, start2 = as.numeric(b$pos), strand2 = b$strand,
    svclass = .sim_svclass(a$chrom, a$strand, b$chrom, b$strand),
    multiplicity = as.numeric(multiplicity),
    round = as.integer(round), mechanism = mechanism,
    origin = origin, homology = homology
  )
}

## TRUE when fragments i and j are contiguous in the reference AND both forward,
## i.e. the "junction" between them merely continues the reference and would not
## be called as an SV. Shattering-and-rejoining can restore such a join by chance;
## it must not be emitted as a breakpoint.
.sim_is_reference_join <- function(fi, fj) {
  fi$inv[1] == 0 && fj$inv[1] == 0 &&
    as.character(GenomeInfoDb::seqnames(fi))[1] == as.character(GenomeInfoDb::seqnames(fj))[1] &&
    GenomicRanges::start(fj)[1] == GenomicRanges::end(fi)[1] + 1L
}

## Read every junction off a fragment list: consecutive pairs, plus the closing
## pair (last -> first) when the amplicon is circular. Junctions traversed more
## than once per circle are collapsed and their `multiplicity` summed, so a
## repeated fragment order raises the junction copy number rather than duplicating
## the row (junction copy number is what PURPLE's JCN reports).
.sim_junctions_from_fragments <- function(frags, circular, mechanism, round,
                                          homology = "nhej") {
  n <- length(frags)
  if (n == 0L) return(.sim_empty_junctions())
  pairs <- if (n >= 2L) cbind(seq_len(n - 1L), 2:n) else matrix(nrow = 0, ncol = 2)
  if (circular && n >= 2L) pairs <- rbind(pairs, c(n, 1L))
  if (circular && n == 1L) pairs <- rbind(pairs, c(1L, 1L))   # self-circularisation
  if (!nrow(pairs)) return(.sim_empty_junctions())

  rows <- lapply(seq_len(nrow(pairs)), function(k) {
    fi <- frags[pairs[k, 1]]; fj <- frags[pairs[k, 2]]
    if (.sim_is_reference_join(fi, fj)) return(NULL)
    .sim_junction(.sim_exit(fi), .sim_entry(fj), mechanism, round, homology)
  })
  out <- data.table::rbindlist(Filter(Negate(is.null), rows))
  if (!nrow(out)) return(.sim_empty_junctions())
  .sim_collapse_junctions(out)
}

## Collapse identical junctions, summing multiplicity. Keeps the FIRST mechanism /
## round label, so a junction created at origin and re-traversed later is still
## attributed to the event that made it.
.sim_collapse_junctions <- function(j) {
  if (!nrow(j)) return(j)
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  j <- data.table::as.data.table(j)
  j[, .(multiplicity = sum(multiplicity), svclass = svclass[1],
        round = round[1], mechanism = mechanism[1],
        origin = origin[1], homology = homology[1]),
    by = key]
}

## Split a fragment list at `n` uniformly random positions, choosing the fragment
## to cut with probability proportional to its width so that breakpoints are
## uniform over SEQUENCE rather than over fragments. Returns the fragment list
## with the cut fragments split in place (amplicon order preserved).
.sim_shatter <- function(frags, n_breaks, min_width = 4L) {
  if (n_breaks <= 0L || length(frags) == 0L) return(frags)
  for (i in seq_len(n_breaks)) {
    w <- GenomicRanges::width(frags)
    idx <- which(w > min_width)
    if (!length(idx)) break
    r <- if (length(idx) == 1L) idx else sample(idx, 1L, prob = w[idx])
    ## draw the cut offset rather than materialising the whole coordinate range
    ## (a fragment can be megabases wide, and this runs once per break)
    bp <- GenomicRanges::start(frags)[r] + sample.int(w[r] - 2L, 1L)
    left <- frags[r]; right <- frags[r]
    GenomicRanges::end(left) <- bp
    GenomicRanges::start(right) <- bp + 1L
    before <- if (r > 1L) frags[seq_len(r - 1L)] else frags[0]
    after  <- if (r < length(frags)) frags[(r + 1L):length(frags)] else frags[0]
    frags <- c(before, left, right, after)
  }
  frags
}

## Fragments overlapping any preserved (oncogene) locus must survive shattering:
## an ecDNA that loses its oncogene confers no selective advantage and would not
## be observed. Mirrors the `preservedRegions` handling in the Seismic simulator.
.sim_is_preserved <- function(frags, preserved) {
  if (is.null(preserved) || !length(preserved) || !length(frags))
    return(rep(FALSE, length(frags)))
  GenomicRanges::countOverlaps(frags, preserved) > 0L
}

## Copy-number profile of a fragment list: coverage() gives the per-base
## multiplicity of each locus within ONE circle; multiplying by `copies` gives the
## amplicon's contribution to the tumour copy number. Returned as segments.
.sim_fragment_coverage <- function(frags) {
  if (!length(frags)) return(GenomicRanges::GRanges())
  cov <- GenomicRanges::coverage(frags)
  out <- lapply(names(cov), function(chr) {
    r <- cov[[chr]]
    if (methods::is(r, "Rle") && length(r) == 0L) return(NULL)
    ends <- cumsum(S4Vectors::runLength(r))
    starts <- c(1L, utils::head(ends, -1L) + 1L)
    v <- S4Vectors::runValue(r)
    keep <- v > 0
    if (!any(keep)) return(NULL)
    GenomicRanges::GRanges(chr, IRanges::IRanges(starts[keep], ends[keep]),
                           n = as.numeric(v[keep]))
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(GenomicRanges::GRanges())
  do.call(c, out)
}

## ---- population of circle variants -----------------------------------------
## An amplicon is usually one structure present at many copies, but once a single
## molecule acquires a change (an internal deletion, say) the cell carries a
## MIXTURE: the intact circle and the altered one, at different abundances. That
## mixture is what makes the ORDER of events observable -- a junction that arose
## late sits at lower copy number, and therefore lower read support, than one
## present from the start.
##
## `variants` holds that mixture as a list of `list(fragments, copies, label)`.
## It is NULL for a homogeneous amplicon, in which case the population is just
## the single (fragments, copies) pair, so everything behaves as before.
.sim_population <- function(x) {
  pop <- if (!is.null(x$variants) && length(x$variants)) x$variants
         else list(list(fragments = x$fragments, copies = x$copies,
                        label = "founder"))
  lapply(pop, function(v) { if (is.null(v$fitness)) v$fitness <- 1; v })
}

## Write a population back onto the object, keeping the representative fields in
## step: `fragments` is the MOST ABUNDANT species' structure (what print(),
## summary() and the mechanism operators act on) and `copies` the total across
## species. Species that have died out are dropped.
.sim_set_population <- function(x, pop) {
  pop <- Filter(function(v) isTRUE(v$copies > 0), pop)
  stopifnot(length(pop) >= 1L)
  cp <- vapply(pop, function(v) v$copies, numeric(1))
  x$variants   <- if (length(pop) > 1L) pop else NULL
  x$fragments  <- pop[[which.max(cp)]]$fragments
  x$copies     <- sum(cp)
  x
}

## Total ecDNA contribution to copy number at each locus: the per-base coverage
## of every species weighted by its abundance, summed. For a homogeneous amplicon
## this is simply coverage(fragments) * copies.
.sim_amplicon_cn <- function(x) {
  pop <- .sim_population(x)
  parts <- Filter(length, lapply(pop, function(v) {
    g <- .sim_fragment_coverage(v$fragments)
    if (length(g)) g$n <- g$n * v$copies
    g
  }))
  if (!length(parts)) return(GenomicRanges::GRanges())
  if (length(parts) == 1L) return(parts[[1]])

  bare <- lapply(parts, function(g) { S4Vectors::mcols(g) <- NULL; g })
  seg <- GenomicRanges::disjoin(do.call(c, bare))
  mid <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(seg),
                                IRanges::IRanges(GenomicRanges::start(seg), width = 1L))
  seg$n <- Reduce(`+`, lapply(parts, function(p) .sim_lookup(mid, p, "n", 0)))
  seg[seg$n > 0]
}

## Junction copy number per junction across the population: a junction is carried
## only by the species whose structure traverses it, so its copy number is the
## sum of `multiplicity * copies` over those species. This is what puts a late
## internal deletion below the boundary duplication it sits inside.
.sim_junction_jcn <- function(x) {
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  parts <- Filter(Negate(is.null), lapply(.sim_population(x), function(v) {
    j <- .sim_junctions_from_fragments(v$fragments, circular = TRUE,
                                       mechanism = "", round = 0L)
    if (!nrow(j)) return(NULL)
    j[, .(jcn = sum(multiplicity) * v$copies), by = key]
  }))
  if (!length(parts)) return(data.table::data.table())
  data.table::rbindlist(parts)[, .(jcn = sum(jcn)), by = key]
}

## Collapse runs in the mechanism history for display: seven replication rounds
## read as "replication x7" rather than seven repeated words.
.sim_history_str <- function(h) {
  if (!length(h)) return("")
  r <- rle(h)
  paste(ifelse(r$lengths > 1L, sprintf("%s x%d", r$values, r$lengths), r$values),
        collapse = " -> ")
}

## Construct the object. Not exported: users enter through sim_episome().
.new_ecdna_sim <- function(sample, fragments, circular, copies, junctions,
                           host, preserved, history, round = 0L, excised = NULL,
                           variants = NULL) {
  structure(list(
    sample = sample, fragments = fragments, circular = circular,
    copies = copies, junctions = junctions, host = host,
    preserved = preserved, history = history, round = as.integer(round),
    excised = excised, variants = variants
  ), class = "ecdna_sim")
}

#' @export
print.ecdna_sim <- function(x, ...) {
  chrs <- unique(as.character(GenomeInfoDb::seqnames(x$fragments)))
  span <- sum(as.numeric(GenomicRanges::width(x$fragments)))
  cat("<ecdna_sim>", x$sample, "\n")
  cat("  structure   :", if (isTRUE(x$circular)) "circular (episome)" else "linear", "\n")
  cat("  fragments   :", length(x$fragments), "on", paste(chrs, collapse = ", "), "\n")
  cat("  circle size :", format(round(span / 1e6, 2), nsmall = 2), "Mb\n")
  cat("  copies/cell :", round(x$copies, 1), "\n")
  pop <- .sim_population(x)
  if (length(pop) > 1L) {
    cat("  population  :", length(pop), "circle species\n")
    for (v in pop)
      cat(sprintf("      %-14s %7.1f copies (%4.1f%%)  fitness %.2f  %3d frag  %6.0f kb\n",
                  v$label, v$copies, 100 * v$copies / x$copies, v$fitness,
                  length(v$fragments),
                  sum(as.numeric(GenomicRanges::width(v$fragments))) / 1e3))
  }
  cat("  junctions   :", nrow(x$junctions),
      if (nrow(x$junctions))
        paste0("(", paste(sprintf("%s=%d", names(table(x$junctions$svclass)),
                                  as.integer(table(x$junctions$svclass))),
                          collapse = ", "), ")") else "", "\n")
  cat("  rounds      :", x$round, "\n")
  cat("  history     :", .sim_history_str(x$history), "\n")
  invisible(x)
}

#' Summarise a simulated amplicon
#'
#' Ground-truth summary statistics of an [ecdna_sim] object -- the quantities the
#' mechanism callers try to recover, so a simulation can be checked against what
#' EpiTracer reports for it.
#'
#' @param object An [ecdna_sim] from [sim_episome()] and friends.
#' @param ... Ignored.
#' @return A one-row [data.table::data.table] with the sample, structure,
#'   fragment and junction counts, per-class junction counts, circle size,
#'   copies per cell, peak copy number, number of chromosomes involved and the
#'   mechanism history.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01")
#' summary(ec)
#' @export
summary.ecdna_sim <- function(object, ...) {
  j <- object$junctions
  cls <- function(k) sum(j$svclass == k)
  cn <- .sim_amplicon_cn(object)
  pop <- .sim_population(object)
  data.table::data.table(
    sample     = object$sample,
    circular   = isTRUE(object$circular),
    n_fragments = length(object$fragments),
    n_chr      = length(unique(as.character(GenomeInfoDb::seqnames(object$fragments)))),
    circle_mb  = sum(as.numeric(GenomicRanges::width(object$fragments))) / 1e6,
    copies     = object$copies,
    n_variants = length(pop),
    max_cn     = if (length(cn)) max(cn$n) else 0,
    n_junctions = nrow(j),
    n_del = cls("DEL"), n_dup = cls("DUP"),
    n_h2h = cls("h2hINV"), n_t2t = cls("t2tINV"), n_tra = cls("TRA"),
    rounds     = object$round,
    history    = .sim_history_str(object$history)
  )
}

#' A seed locus for simulation
#'
#' Convenience wrapper that resolves an oncogene symbol to the amplicon seed and
#' preserved-region pair that [sim_episome()] expects, using the oncogene panel
#' bundled with EpiTracer. The returned region is the gene span padded by `flank`
#' -- a plausible excised neighbourhood -- and carries the gene itself as the
#' locus that must be preserved through shattering.
#'
#' @param gene Gene symbol(s), e.g. `"EGFR"` or `c("CDK4", "MDM2")`.
#' @param genome One of `"hg38"`, `"hg19"`, `"mm10"`.
#' @param flank Bp of neighbouring sequence co-excised with the gene
#'   (default `2.5e5`, giving a sub-megabase episome typical of a simple
#'   excision).
#' @return A [GenomicRanges::GRanges] of the seed region(s), with a `preserved`
#'   attribute holding the gene spans and a `gene` metadata column.
#' @examples
#' seed_locus("EGFR")
#' attr(seed_locus("EGFR"), "preserved")
#' @seealso [sim_episome()], [gene_locus()]
#' @export
seed_locus <- function(gene, genome = c("hg38", "hg19", "mm10"), flank = 2.5e5) {
  genome <- match.arg(genome)
  loc <- gene_locus(gene, genome = genome, flank = 0)
  seed <- GenomicRanges::GRanges(
    loc$chr, IRanges::IRanges(pmax(1, round(loc$start - flank)), round(loc$end + flank)),
    gene = loc$gene
  )
  attr(seed, "preserved") <- GenomicRanges::GRanges(
    loc$chr, IRanges::IRanges(loc$start, loc$end), gene = loc$gene
  )
  seed
}
