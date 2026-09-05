## ---------------------------------------------------------------------------
## Mechanism operators: the evolutionary trajectory of an episomal ecDNA.
##
##   sim_episome()        birth -- excision of a locus and its circularisation
##   sim_micronucleation()one round of micronuclear encapsulation + chromothripsis
##   sim_fuse_episomes()  two circles co-encapsulated -> one chimeric circle
##   sim_evolve()         several rounds of the above
##
## Each takes and returns an `ecdna_sim`, so they compose: an amplicon can be
## born, shattered once, shattered again, and fused with a second episome, and
## the junction ledger accumulates the scars of every step with the round and
## mechanism that produced it. That ledger is the ground truth against which
## call_simple_excision(), call_chimeric_amplicon() and call_chromothripsis() are
## scored.
##
## The biology being modelled (Shoshani et al., Nature 2021; Rosswog et al., Nat
## Genet 2021; Bernhard et al. 2022): an ecDNA lacking a centromere missegregates
## into a micronucleus, whose envelope ruptures; the enclosed DNA is shattered and
## religated at random by non-homologous end joining, and -- because the circle is
## selected on its oncogene dosage -- the products that keep the oncogene sweep.
## Repeating this turns a clean, single-junction episome into the heavily
## rearranged, copy-number-oscillating amplicon seen in patients.
## ---------------------------------------------------------------------------

#' Simulate the birth of an episomal ecDNA
#'
#' Creates a simple-excision episome: one or more genomic regions are excised
#' from the chromosome and circularised. This is the mechanism
#' [call_simple_excision()] is built to detect, and the simulated amplicon
#' carries its three signatures by construction --- a **duplication-orientation
#' junction at the amplicon boundary** (the circularisation junction, `-/+`,
#' emitted as `svclass == "DUP"`), **non-gained flanks** (the host chromosome
#' stays at its baseline ploidy immediately outside the amplicon), and an
#' **excision scar** (a low-copy deletion junction joining the flanks on the
#' chromosome the circle left behind).
#'
#' @param regions A [GenomicRanges::GRanges] of the region(s) to excise. Several
#'   regions produce a circle assembled from several segments, as when an
#'   episome forms from a locus already broken by internal deletions. Use
#'   [seed_locus()] to build one from a gene symbol.
#' @param sample Character sample identifier carried into the emitted inputs.
#' @param copies Numeric ecDNA copies per tumour cell (default `40`). Copy number
#'   is `copies * coverage(fragments)` over the host baseline, so this sets the
#'   amplification level; patient ecDNA amplicons typically sit at 20-100.
#' @param preserved Oncogene loci that later shattering must not delete. Defaults
#'   to the `preserved` attribute of `regions` when [seed_locus()] made it, else
#'   to `regions` themselves.
#' @param host_ploidy Numeric baseline copy number of the host chromosome
#'   (default `2`). Set higher to model an episome on a polysomic chromosome --
#'   e.g. `3` for the chr7 gain that accompanies EGFR ecDNA in glioblastoma.
#' @param host_flank Bp of host chromosome to represent on each side of the
#'   amplicon (default `1.5e7`), giving the flanks the callers test.
#' @param scar Logical; emit the chromosomal excision scar (default `TRUE`). Set
#'   `FALSE` to model an episome whose source chromosome was lost, so the scar is
#'   unobservable.
#' @param retain_source Logical; does the source chromosome keep the excised
#'   region on its other allele (default `TRUE`)? When `TRUE` the excised
#'   interval sits one copy below the flanks on the chromosome, the usual
#'   single-allele excision.
#' @param homology Repair pathway that sealed the circularisation junction, one
#'   of `"nhej"` (blunt), `"mmej"` (2-10 bp microhomology) or `"nahr"`
#'   (>=14 bp homology). Sets `HOMLEN` on the emitted breakends, the
#'   breakpoint-homology signature `call_simple_excision()` reports.
#' @param genome Genome build for chromosome lengths (default `"hg38"`).
#' @param seed Optional integer passed to [set.seed()] for reproducibility.
#' @return An [ecdna_sim] object.
#' @references Shoshani, O. *et al.* Chromothripsis drives the evolution of gene
#'   amplification in cancer. *Nature* **591**, 137-141 (2021).
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01", copies = 50)
#' ec
#' summary(ec)
#' @seealso [sim_micronucleation()], [sim_to_epitracer()], [call_simple_excision()]
#' @export
sim_episome <- function(regions, sample = "SIM01", copies = 40,
                        preserved = NULL, host_ploidy = 2, host_flank = 1.5e7,
                        scar = TRUE, retain_source = TRUE,
                        homology = c("nhej", "mmej", "nahr"),
                        genome = "hg38", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)
  stopifnot(methods::is(regions, "GRanges"), length(regions) >= 1L)

  if (is.null(preserved)) {
    preserved <- attr(regions, "preserved")
    if (is.null(preserved)) preserved <- regions
  }

  regions_plain <- regions
  S4Vectors::mcols(regions_plain) <- NULL
  attr(regions_plain, "preserved") <- NULL

  frags <- regions
  S4Vectors::mcols(frags) <- NULL
  frags$inv <- 0
  frags$origin <- sample

  ## The circle's own junctions: the closing (circularisation) junction plus any
  ## junction implied by joining several excised regions head to tail.
  j <- .sim_junctions_from_fragments(frags, circular = TRUE,
                                     mechanism = "simple_excision", round = 0L,
                                     homology = homology)

  ## The excision scar left on the chromosome: the flanks either side of the
  ## excised interval are joined, a deletion-orientation junction at low copy
  ## number. One scar per excised interval.
  if (scar && retain_source) {
    scars <- data.table::rbindlist(lapply(seq_along(regions), function(i) {
      chr <- as.character(GenomeInfoDb::seqnames(regions))[i]
      lo <- GenomicRanges::start(regions)[i]; hi <- GenomicRanges::end(regions)[i]
      .sim_junction(list(chrom = chr, pos = lo - 1L, strand = "+"),
                    list(chrom = chr, pos = hi + 1L, strand = "-"),
                    mechanism = "excision_scar", round = 0L,
                    homology = homology, origin = "host")
    }))
    j <- rbind(j, scars)
  }

  host <- .sim_host_background(regions, host_ploidy, host_flank, genome)

  .new_ecdna_sim(
    sample = sample, fragments = frags, circular = TRUE, copies = copies,
    junctions = j, host = host, preserved = preserved,
    history = "simple_excision", round = 0L,
    excised = if (retain_source) regions_plain else NULL
  )
}

## Host-chromosome background: one window per involved chromosome, spanning the
## amplicon padded by `flank` and clipped to the chromosome. Establishes the local
## copy-number baseline that the callers' flank tests read.
.sim_host_background <- function(regions, ploidy, flank, genome) {
  L <- tryCatch(load_chrom_lengths(genome), error = function(e) numeric(0))
  chrs <- unique(as.character(GenomeInfoDb::seqnames(regions)))
  gr <- do.call(c, lapply(chrs, function(ch) {
    on <- as.character(GenomeInfoDb::seqnames(regions)) == ch
    lo <- max(1, min(GenomicRanges::start(regions)[on]) - flank)
    hi <- max(GenomicRanges::end(regions)[on]) + flank
    if (ch %in% names(L)) hi <- min(hi, L[[ch]])
    GenomicRanges::GRanges(ch, IRanges::IRanges(lo, hi), ploidy = ploidy)
  }))
  gr
}

#' Simulate one round of micronucleation and chromothripsis
#'
#' Encapsulates the ecDNA in a micronucleus, shatters it, and religates the
#' surviving fragments in random order and orientation back into a circle --- the
#' ecDNA -> micronucleus -> chromothripsis route by which a clean episome becomes
#' a complex amplicon.
#'
#' Three things happen, in order:
#'
#' 1. **Shattering.** The circle is cut at `n_breaks` positions drawn uniformly
#'    over its sequence (so wide fragments are cut more often than narrow ones).
#' 2. **Fragment loss and reduplication.** Each fragment is lost with probability
#'    `del_p`, except fragments carrying a preserved oncogene, which always
#'    survive --- a circle that loses its oncogene is not selected. Each surviving
#'    fragment is incorporated `1 + rbinom(1, max_dup - 1, dup_p)` times. This
#'    reduplication is what generates the **oscillating copy-number profile**
#'    that distinguishes a chromothriptic amplicon from a clean episome (which
#'    would otherwise be a flat two-level profile), and it is the step
#'    [call_chromothripsis()] scores as `cn_oscillations`.
#' 3. **Random rejoining.** Fragments are shuffled and each is inverted with
#'    probability `inv_p`, then religated into a circle. Because order and
#'    orientation are random, the four intrachromosomal junction orientations
#'    (`DEL`, `DUP`, `h2hINV`, `t2tINV`) come out close to equally represented ---
#'    the **random fragment joins** hallmark, which [call_chromothripsis()] tests
#'    with a chi-squared goodness-of-fit and which separates chromothripsis from
#'    orientation-biased mechanisms such as BFB.
#'
#' Applying this repeatedly (see [sim_evolve()]) compounds the rearrangement:
#' junction counts grow roughly linearly in the number of rounds while the circle
#' contracts, reproducing the small, junction-dense amplicons seen in patients.
#'
#' @section Circle populations:
#' A shattering event happens inside ONE micronucleus, to ONE molecule, and the
#' surviving product is what sweeps. So when the amplicon is a mixture of circle
#' species (see [sim_internal_deletion()]), this operator acts on the most
#' abundant species and the population collapses back to one --- the shattered
#' descendant. Minority species present before the event are not carried through.
#'
#' @param x An [ecdna_sim], from [sim_episome()] or a previous round.
#' @param n_breaks Integer number of chromothriptic breakpoints (default `10`).
#'   Patient chromothripsis footprints typically carry 10-50.
#' @param del_p Probability a non-preserved fragment is lost (default `0.2`).
#' @param inv_p Probability a surviving fragment is inverted (default `0.5`).
#' @param dup_p Per-trial probability of an extra incorporation of a fragment
#'   (default `0.25`), giving the copy-number oscillation described above. Set
#'   `0` for strict shatter-and-rejoin with no reduplication.
#' @param max_dup Maximum times one fragment may appear in the religated circle
#'   (default `3`).
#' @param amplify Multiplicative change in copies per cell across the round
#'   (default `1.5`), modelling selection for the fitter, oncogene-retaining
#'   circle. Set `1` for no amplification.
#' @param max_cn Ceiling on the amplicon's peak copy number (default `200`).
#'   Copies per cell are scaled back so `copies * max(coverage)` stays at or
#'   below it, modelling the plateau beyond which further oncogene dosage buys no
#'   selective advantage (the `N_fitness` saturation of the Bernhard *et al.*
#'   model). Without it, compounding `amplify` and fragment reduplication drive
#'   copy number far past anything observed in patients.
#' @param homology Repair pathway sealing the new junctions (default `"nhej"`;
#'   micronuclear religation is canonically blunt-ended NHEJ).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim] with the religated fragment list, the accumulated
#'   junction ledger and `round` incremented.
#' @references Zhang, C.-Z. *et al.* Chromothripsis from DNA damage in
#'   micronuclei. *Nature* **522**, 179-184 (2015).
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01")
#' mn <- sim_micronucleation(ec, n_breaks = 12, seed = 1)
#' summary(mn)
#' @seealso [sim_evolve()], [sim_fuse_episomes()], [call_chromothripsis()]
#' @export
sim_micronucleation <- function(x, n_breaks = 10L, del_p = 0.2, inv_p = 0.5,
                                dup_p = 0.25, max_dup = 3L, amplify = 1.5,
                                max_cn = 200,
                                homology = c("nhej", "mmej", "nahr"),
                                seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"))
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)
  rnd <- x$round + 1L

  frags <- .sim_reassemble(x$fragments, x$preserved, n_breaks = n_breaks,
                           del_p = del_p, inv_p = inv_p,
                           dup_p = dup_p, max_dup = max_dup)

  ## Junctions of the religated circle. Junctions from earlier rounds that the
  ## new order happens to preserve reappear here and are merged with the ledger,
  ## keeping their original round/mechanism attribution.
  new_j <- .sim_junctions_from_fragments(frags, circular = TRUE,
                                         mechanism = "micronucleation",
                                         round = rnd, homology = homology)
  j <- .sim_merge_ledger(x$junctions, new_j)

  .new_ecdna_sim(
    sample = x$sample, fragments = frags, circular = TRUE,
    copies = .sim_cap_copies(x$copies * amplify, frags, max_cn),
    junctions = j, host = x$host,
    preserved = x$preserved, excised = x$excised,
    history = c(x$history, "micronucleation"), round = rnd
  )
}

## Scale copies per cell back so the amplicon's PEAK copy number (copies times
## the highest per-circle fragment multiplicity) respects `max_cn`. Reduplication
## raises copy number for free, so without this the level compounds without bound
## across rounds.
.sim_cap_copies <- function(copies, frags, max_cn) {
  if (!is.finite(max_cn) || max_cn <= 0 || !length(frags)) return(copies)
  cov <- .sim_fragment_coverage(frags)
  peak <- if (length(cov)) max(cov$n) else 1
  min(copies, max_cn / peak)
}

## Shatter -> lose -> reduplicate -> shuffle -> invert. Shared by
## sim_micronucleation() and sim_fuse_episomes().
.sim_reassemble <- function(frags, preserved, n_breaks, del_p, inv_p,
                            dup_p, max_dup) {
  frags <- .sim_shatter(frags, n_breaks)
  if (!length(frags)) return(frags)

  ## Fragment loss, sparing anything carrying a preserved oncogene.
  keep_gene <- .sim_is_preserved(frags, preserved)
  lost <- !keep_gene & stats::runif(length(frags)) < del_p
  if (all(lost)) lost[sample(which(lost), 1L)] <- FALSE   # never lose everything
  frags <- frags[!lost]
  if (!length(frags)) return(frags)

  ## Reduplication: extra incorporations of surviving fragments, the source of
  ## the multi-level, oscillating copy-number profile.
  n_copies <- 1L + stats::rbinom(length(frags), max(0L, as.integer(max_dup) - 1L), dup_p)
  frags <- rep(frags, n_copies)

  ## Random order and random orientation -- the chromothripsis signature.
  frags <- frags[sample(length(frags))]
  flip <- stats::runif(length(frags)) < inv_p
  frags$inv[flip] <- 1 - frags$inv[flip]
  frags
}

## Merge a new junction set into the accumulated ledger. A junction present in
## both keeps its ORIGINAL round/mechanism (the event that created it) but takes
## the new multiplicity (its traversals in the current structure); a junction only
## in the old ledger is retained when it lives on the host chromosome (scars
## persist) and dropped otherwise -- an amplicon junction not traversed by the
## current circle has been recombined away.
.sim_merge_ledger <- function(old, new) {
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  old <- data.table::as.data.table(old); new <- data.table::as.data.table(new)
  host <- old[origin == "host"]
  amp_old <- old[origin != "host"]
  if (nrow(new) && nrow(amp_old)) {
    m <- merge(new, amp_old[, c(key, "round", "mechanism"), with = FALSE],
               by = key, all.x = TRUE, suffixes = c("", ".old"))
    hit <- !is.na(m$round.old)
    m$round[hit] <- m$round.old[hit]
    m$mechanism[hit] <- m$mechanism.old[hit]
    new <- m[, names(new), with = FALSE]
  }
  data.table::rbindlist(list(host, new), use.names = TRUE, fill = TRUE)
}

## Rotate a shattered fragment list so that the circle CLOSES on a specified
## junction: the fragment whose reference start is `entry` is placed first, the
## fragment whose reference end is `exit` is placed last, both forward, and only
## the fragments between them are reordered and reoriented. The closing
## (last -> first) junction then reproduces the preserved junction exactly.
## Returns NULL when the junction's ends are no longer fragment boundaries (it
## was cut away by an earlier lossy round), so the caller can fall back.
.sim_reclose_on <- function(frags, entry, exit, inv_p) {
  st <- GenomicRanges::start(frags); en <- GenomicRanges::end(frags)
  i_first <- which(st == entry); i_last <- which(en == exit)
  if (!length(i_first) || !length(i_last)) return(NULL)
  i_first <- i_first[1]; i_last <- i_last[1]
  if (i_first == i_last) return(NULL)          # a single fragment spans both ends

  middle <- setdiff(seq_along(frags), c(i_first, i_last))
  if (length(middle)) {
    middle <- middle[sample.int(length(middle))]
    flip <- stats::runif(length(middle)) < inv_p
    frags$inv[middle[flip]] <- 1 - frags$inv[middle[flip]]
  }
  out <- c(frags[i_first], frags[middle], frags[i_last])
  out$inv[1] <- 0                              # entry fragment read forwards
  out$inv[length(out)] <- 0                    # exit fragment read forwards
  out
}

#' Reconfigure an ecDNA while preserving its founder junction
#'
#' Shatters a circle and religates it into a NEW circle, but keeps a nominated
#' junction --- by default the oldest amplicon junction, i.e. the founder
#' circularisation duplication --- intact as the closing junction. Everything
#' between its two ends is cut, reordered and reoriented at random.
#'
#' The preservation is **imposed, not derived**. This operator exists to explore
#' the counterfactual --- what an amplicon would look like if one junction were
#' carried through repeated rebuilding --- and the resulting population shows a
#' junction hierarchy ordered by age, with the preserved junction clonal and each
#' later round a step below, on a population in which nothing is unshattered.
#'
#' @section There is no known selection for a specific junction:
#' Do not read the preservation as a mechanism. A religated circle is a circle
#' whichever junction closes it: circularity is topological, so a descendant that
#' loses the founder duplication is no less viable than one that keeps it, and
#' nothing distinguishes an arbitrary duplication between two arbitrary ends as a
#' target of selection. Under free reshuffling the founder survives a round with
#' probability about `(1 - del_p) / (2 * n_breaks)`, so its chance of remaining
#' clonal across several rounds is negligible.
#'
#' A **clonal** founder duplication in real data is therefore evidence that the
#' molecules carrying it were never extensively reshuffled --- the parsimonious
#' reading being that chromothripsis struck a minority of the population and the
#' rest was untouched. The one indirect route to genuine selection is on the
#' whole CONFIGURATION rather than the junction: retaining the founder junction
#' amounts to retaining the excised segment contiguous and in reference
#' orientation, which preserves the oncogene's regulatory neighbourhood, and
#' ecDNA is known to carry co-amplified enhancers whose spacing to the promoter
#' matters. But that argument selects for the intact arrangement as a whole, and
#' so would preserve the internal junctions too --- not the founder alone, which
#' is what this operator does. Treat `preserve` as a hypothesis to test, not as
#' biology already established.
#'
#' @param x An [ecdna_sim].
#' @param n_breaks Integer shattering breakpoints applied to the interior
#'   (default `10`).
#' @param copies Copies per cell of the reconfigured lineage (default `NULL`,
#'   meaning it takes the whole source species --- a clean sweep). Give a number
#'   to have it arise as a subclone alongside its parent.
#' @param preserve Which junction to keep: `"founder"` (the oldest amplicon
#'   junction, the default), `"none"`, or a one-row `data.frame`/`data.table`
#'   with `start1` and `start2` naming the junction's two breakends.
#' @param del_p,inv_p Fragment loss and inversion probabilities for the interior
#'   (defaults `0.05` and `0.5`). Fragments carrying a preserved oncogene are
#'   never lost.
#' @param fitness Replication rate of the new lineage relative to its parent
#'   (default `1`); see [sim_internal_deletion()].
#' @param label Name for the new species (default `"reconfigured"`).
#' @param homology Repair pathway sealing the new junctions (default `"nhej"`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim]. When the nominated junction's ends have been cut away
#'   by an earlier round it cannot be preserved; the reconfiguration still
#'   happens and the result carries the attribute `"preserved_founder" = FALSE`.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR", flank = 1.2e6), sample = "RC01", copies = 90)
#' ## three rounds of reconfiguration, each a subclone, nothing left unshattered
#' for (i in 1:3) ec <- sim_reconfigure(ec, n_breaks = 12, copies = 12, seed = i)
#' ec
#' @seealso [sim_micronucleation()], [sim_shatter_to_episomes()],
#'   [plot_sv_reconstruction()]
#' @export
sim_reconfigure <- function(x, n_breaks = 10L, copies = NULL,
                            preserve = "founder", del_p = 0.05, inv_p = 0.5,
                            fitness = 1, label = "reconfigured",
                            homology = c("nhej", "mmej", "nahr"), seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"))
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)
  rnd <- x$round + 1L

  pop <- .sim_population(x)
  src <- which.max(vapply(pop, function(v) v$copies, numeric(1)))
  parent <- pop[[src]]
  take <- if (is.null(copies)) parent$copies else copies
  if (take > parent$copies)
    stop("`copies` (", take, ") exceeds the ", parent$copies,
         " copies available in the species being reconfigured.", call. = FALSE)

  ## which junction to carry through
  keep <- if (identical(preserve, "none")) NULL
          else if (identical(preserve, "founder")) {
            amp <- x$junctions[origin == "amplicon"]
            if (nrow(amp)) amp[order(round)][1] else NULL
          } else data.table::as.data.table(preserve)[1]

  frags <- .sim_shatter(parent$fragments, n_breaks)

  ## Fragments carrying a preserved oncogene are never lost -- nor are the two
  ## fragments carrying the ends of the junction being preserved: deleting them
  ## would destroy the very junction this operator exists to carry through, and
  ## a descendant that loses its circularisation junction is not a viable circle.
  protect <- .sim_is_preserved(frags, x$preserved)
  if (!is.null(keep))
    protect <- protect |
      GenomicRanges::start(frags) == keep$start1 |
      GenomicRanges::end(frags)   == keep$start2
  lost <- !protect & stats::runif(length(frags)) < del_p
  if (all(lost)) lost[sample.int(length(lost), 1L)] <- FALSE
  frags <- frags[!lost]
  if (length(frags) < 2L) stop("too little material left to reconfigure.", call. = FALSE)

  new_frags <- if (!is.null(keep))
    .sim_reclose_on(frags, entry = keep$start1, exit = keep$start2, inv_p = inv_p)
  else NULL
  preserved_ok <- !is.null(new_frags)

  if (!preserved_ok) {
    ## fall back to a free reconfiguration (the junction's ends are gone)
    new_frags <- frags[sample.int(length(frags))]
    flip <- stats::runif(length(new_frags)) < inv_p
    new_frags$inv[flip] <- 1 - new_frags$inv[flip]
  }

  new_j <- .sim_junctions_from_fragments(new_frags, circular = TRUE,
                                         mechanism = "reconfiguration",
                                         round = rnd, homology = homology)
  x$junctions <- .sim_merge_ledger_keep(x$junctions, new_j)

  pop[[src]]$copies <- parent$copies - take
  pop <- c(pop, list(list(fragments = new_frags, copies = take,
                          fitness = fitness, label = label)))
  x$history <- c(x$history, "reconfiguration")
  x$round <- rnd
  out <- .sim_set_population(x, pop)
  attr(out, "preserved_founder") <- preserved_ok
  out
}

## Ledger merge that RETAINS every junction already recorded, rather than
## dropping those the current structure no longer traverses. Reconfiguration
## builds a population of coexisting lineages, so a junction absent from the
## newest lineage may still be carried by an older one; .sim_junction_jcn()
## works out which lineages actually carry it at emission time.
.sim_merge_ledger_keep <- function(old, new) {
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  old <- data.table::as.data.table(old); new <- data.table::as.data.table(new)
  if (!nrow(new)) return(old)
  known <- do.call(paste, c(old[, key, with = FALSE], sep = "\r"))
  fresh <- new[!do.call(paste, c(new[, key, with = FALSE], sep = "\r")) %in% known]
  data.table::rbindlist(list(old, fresh), use.names = TRUE, fill = TRUE)
}

#' Simulate shattering into several new episomes
#'
#' The ecDNA is encapsulated in a micronucleus and shattered, and the surviving
#' fragments religate into **several separate circles** rather than one. Each
#' daughter is a new episome with its own circularisation junction, and the
#' population that results is a mixture of them.
#'
#' This is the branch [sim_micronucleation()] does not cover: that operator
#' religates the fragments back into a single circle, whereas here the fragments
#' partition, so one amplicon becomes a family of smaller ones.
#'
#' Two consequences fall out that are worth simulating for:
#'
#' * **The founder boundary duplication is usually destroyed.** It survives only
#'   if some daughter happens to rejoin the two ends it joined, which is unlikely
#'   once the circle has been cut in many places. An amplicon that has been
#'   through fission therefore keeps a high, flat copy number over whatever its
#'   surviving daughters retain, but has no duplication-orientation junction
#'   spanning it --- so [call_simple_excision()] no longer recognises it, even
#'   though it is still, genuinely, a set of episomes.
#' * **Copy number becomes patchy.** Sequence that ends up in a daughter that was
#'   lost, or in one that failed selection, drops back to the host baseline, so
#'   the flat amplicon of the parent breaks into blocks at different levels.
#'
#' @param x An [ecdna_sim]. The most abundant circle species is the one captured
#'   in the micronucleus.
#' @param n_breaks Integer number of shattering breakpoints (default `20`).
#' @param n_circles Number of daughter circles the fragments religate into
#'   (default `3`).
#' @param del_p Probability a fragment is lost outright during religation
#'   (default `0.1`).
#' @param inv_p Probability a fragment is inverted (default `0.5`).
#' @param require_oncogene Logical; drop daughters that carry no preserved locus
#'   (default `TRUE`). A circle without the oncogene confers no advantage and is
#'   diluted out. Set `FALSE` to keep every daughter, modelling the moment before
#'   selection has acted.
#' @param keep_parent Logical; does the intact parent circle survive alongside
#'   the daughters (default `FALSE`)? `FALSE` models the whole ecDNA population
#'   passing through the micronucleus, so the daughters inherit and split its
#'   copies. `TRUE` models a single molecule being shattered while the rest of the
#'   population carries on, giving a parent that keeps its founder junction plus
#'   daughters arising as minority species.
#' @param copies Copies per cell of each daughter when `keep_parent = TRUE`
#'   (default `1`, a single molecule). Ignored when `keep_parent = FALSE`, where
#'   the parent's copies are split evenly across the surviving daughters.
#' @param fitness Replication rate of the daughters relative to the parent
#'   (default `1`). Daughters are shorter, so a value above `1` is the expected
#'   direction; see [sim_internal_deletion()].
#' @param homology Repair pathway sealing the new junctions (default `"nhej"`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim] whose population is the surviving daughters (plus the
#'   parent when `keep_parent = TRUE`). Its `history` gains
#'   `"micronuclear_fission"`.
#' @references Shoshani, O. *et al.* Chromothripsis drives the evolution of gene
#'   amplification in cancer. *Nature* **591**, 137-141 (2021).
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "FIS01", copies = 60)
#' fis <- sim_shatter_to_episomes(ec, n_breaks = 20, n_circles = 3, seed = 1)
#' fis                                   # several daughter species
#' ## the founder boundary duplication is gone
#' fis$junctions[origin == "amplicon" & mechanism == "simple_excision"]
#' @seealso [sim_micronucleation()], [sim_episome()], [call_simple_excision()]
#' @export
sim_shatter_to_episomes <- function(x, n_breaks = 20L, n_circles = 3L,
                                    del_p = 0.1, inv_p = 0.5,
                                    require_oncogene = TRUE, keep_parent = FALSE,
                                    copies = 1, fitness = 1,
                                    homology = c("nhej", "mmej", "nahr"),
                                    seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), n_circles >= 1L)
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)
  rnd <- x$round + 1L

  pop <- .sim_population(x)
  src <- which.max(vapply(pop, function(v) v$copies, numeric(1)))
  parent <- pop[[src]]

  ## shatter, lose a share of the fragments, then deal what survives into
  ## `n_circles` piles -- each pile religates into its own circle
  frags <- .sim_shatter(parent$fragments, n_breaks)
  if (length(frags) < n_circles)
    stop("too few fragments (", length(frags), ") to form ", n_circles,
         " circles; raise `n_breaks`.", call. = FALSE)

  ## fragments carrying a preserved oncogene are never lost outright: an event
  ## that destroys the oncogene leaves nothing selectable, so the products we
  ## observe are conditioned on it surviving (as in .sim_reassemble())
  lost <- !.sim_is_preserved(frags, x$preserved) &
    stats::runif(length(frags)) < del_p
  if (all(lost)) lost[sample.int(length(lost), 1L)] <- FALSE
  frags <- frags[!lost]

  ## every daughter gets at least one fragment, the rest are dealt at random
  assign <- c(seq_len(n_circles),
              sample.int(n_circles, max(0L, length(frags) - n_circles), replace = TRUE))
  assign <- assign[sample.int(length(assign))][seq_len(length(frags))]

  daughters <- Filter(Negate(is.null), lapply(seq_len(n_circles), function(k) {
    f <- frags[assign == k]
    if (!length(f)) return(NULL)
    f <- f[sample.int(length(f))]                      # random order
    flip <- stats::runif(length(f)) < inv_p
    f$inv[flip] <- 1 - f$inv[flip]                     # random orientation
    f
  }))

  ## selection: only circles carrying a preserved oncogene persist
  if (require_oncogene && !is.null(x$preserved) && length(x$preserved)) {
    keep <- vapply(daughters, function(f)
      any(.sim_is_preserved(f, x$preserved)), logical(1))
    if (!any(keep))
      stop("no daughter circle retained a preserved locus; lower `n_circles` ",
           "or `del_p`.", call. = FALSE)
    daughters <- daughters[keep]
  }

  ## abundances
  new_pop <- if (keep_parent) {
    c(pop, lapply(seq_along(daughters), function(k)
      list(fragments = daughters[[k]], copies = copies, fitness = fitness,
           label = paste0("episome", k))))
  } else {
    each <- parent$copies / length(daughters)
    lapply(seq_along(daughters), function(k)
      list(fragments = daughters[[k]], copies = each, fitness = fitness,
           label = paste0("episome", k)))
  }

  ## junctions: the union over the surviving daughters. A junction the parent
  ## already carried keeps its original provenance, so the founder duplication is
  ## still attributed to the excision that made it -- IF any daughter rejoined
  ## those two ends. Usually none does, and it simply disappears from the ledger.
  new_j <- data.table::rbindlist(lapply(daughters, function(f)
    .sim_junctions_from_fragments(f, circular = TRUE,
                                  mechanism = "micronuclear_fission",
                                  round = rnd, homology = homology)),
    use.names = TRUE, fill = TRUE)
  if (keep_parent)
    new_j <- data.table::rbindlist(list(
      .sim_junctions_from_fragments(parent$fragments, circular = TRUE,
                                    mechanism = "micronuclear_fission",
                                    round = rnd, homology = homology),
      new_j), use.names = TRUE, fill = TRUE)
  if (nrow(new_j)) new_j <- .sim_collapse_junctions(new_j)

  x$junctions <- .sim_merge_ledger(x$junctions, new_j)
  x$history <- c(x$history, "micronuclear_fission")
  x$round <- rnd
  .sim_set_population(x, new_pop)
}

#' Simulate chimeric ecDNA formation by co-encapsulation of two episomes
#'
#' Two independently born episomes are captured in the **same** micronucleus,
#' shattered together, and religated into a single chimeric circle carrying
#' fragments of both --- the route by which an ecDNA acquires sequence from a
#' second, non-homologous locus.
#'
#' When the two episomes come from different chromosomes the chimera necessarily
#' carries **inter-chromosomal junctions with both breakends inside amplified
#' copy number**, which is precisely the signature [call_chimeric_amplicon()]
#' detects. It also puts two oncogenes on one circle, the co-amplification
#' (e.g. *EGFR* with *CDK4*, or *MYC* with a distant enhancer) that is common in
#' patient amplicons and cannot arise from a single excision.
#'
#' @section Circle populations:
#' A shattering event happens inside ONE micronucleus, to ONE molecule, and the
#' surviving product is what sweeps. So when the amplicon is a mixture of circle
#' species (see [sim_internal_deletion()]), this operator acts on the most
#' abundant species and the population collapses back to one --- the shattered
#' descendant. Minority species present before the event are not carried through.
#'
#' @param x,y Two [ecdna_sim] objects to co-encapsulate. `x` supplies the sample
#'   name, host background and copy level unless overridden.
#' @param n_breaks Integer breakpoints applied to the pooled fragment set
#'   (default `10`).
#' @param del_p,inv_p,dup_p,max_dup Fragment loss, inversion and reduplication,
#'   as in [sim_micronucleation()].
#' @param copies Copies per cell of the chimera; defaults to the mean of the two
#'   parents, the chimera replacing both.
#' @param sample Sample name for the chimera (defaults to `x`'s).
#' @param homology Repair pathway sealing the new junctions (default `"nhej"`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim] whose fragments carry an `origin` column naming the
#'   parental episome each came from, so the chimeric contribution is traceable.
#' @references Hung, K. L. *et al.* ecDNA hubs drive cooperative intermolecular
#'   oncogene expression. *Nature* **600**, 731-736 (2021).
#' @examples
#' a <- sim_episome(seed_locus("EGFR"), sample = "SIM02", copies = 40)
#' b <- sim_episome(seed_locus("CDK4"), sample = "SIM02", copies = 30)
#' ch <- sim_fuse_episomes(a, b, seed = 1)
#' summary(ch)   # n_chr == 2, n_tra > 0
#' @seealso [sim_micronucleation()], [call_chimeric_amplicon()]
#' @export
sim_fuse_episomes <- function(x, y, n_breaks = 10L, del_p = 0.15, inv_p = 0.5,
                              dup_p = 0.2, max_dup = 3L, copies = NULL,
                              sample = NULL, homology = c("nhej", "mmej", "nahr"),
                              seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), inherits(y, "ecdna_sim"))
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)
  sample <- sample %||% x$sample
  copies <- copies %||% mean(c(x$copies, y$copies))
  rnd <- max(x$round, y$round) + 1L

  ## Label each parent's fragments so the chimeric contribution stays traceable,
  ## then pool and reassemble them as one micronuclear content.
  fx <- x$fragments; fx$origin <- paste0(x$sample, ":A")
  fy <- y$fragments; fy$origin <- paste0(y$sample, ":B")
  pooled <- c(fx, fy)
  preserved <- .sim_c_gr(x$preserved, y$preserved)

  frags <- .sim_reassemble(pooled, preserved, n_breaks = n_breaks,
                           del_p = del_p, inv_p = inv_p,
                           dup_p = dup_p, max_dup = max_dup)
  frags$origin <- as.character(frags$origin)

  new_j <- .sim_junctions_from_fragments(frags, circular = TRUE,
                                         mechanism = "chimeric_fusion",
                                         round = rnd, homology = homology)
  ## Both parents' host scars persist -- each episome left one behind.
  j <- .sim_merge_ledger(rbind(x$junctions, y$junctions), new_j)

  .new_ecdna_sim(
    sample = sample, fragments = frags, circular = TRUE, copies = copies,
    junctions = j, host = .sim_c_gr(x$host, y$host), preserved = preserved,
    excised = .sim_c_gr(x$excised, y$excised),
    history = c(x$history, paste0("fuse(", paste(y$history, collapse = "+"), ")"),
                "chimeric_fusion"),
    round = rnd
  )
}

## Concatenate two GRanges that may carry different metadata columns.
.sim_c_gr <- function(a, b) {
  if (is.null(a) || !length(a)) return(b)
  if (is.null(b) || !length(b)) return(a)
  common <- intersect(names(S4Vectors::mcols(a)), names(S4Vectors::mcols(b)))
  S4Vectors::mcols(a) <- S4Vectors::mcols(a)[, common, drop = FALSE]
  S4Vectors::mcols(b) <- S4Vectors::mcols(b)[, common, drop = FALSE]
  c(a, b)
}

#' Simulate rounds of ecDNA amplification
#'
#' Amplifies a circle **without changing its structure**: the fragment list and
#' the junction ledger are untouched, and only `copies` (ecDNA molecules per
#' tumour cell) rises, by `fold` per round.
#'
#' `fold` is a **net growth factor of copy number under selection, not a cell
#' division**. Cell divisions do not raise ecDNA copy number: with no centromere
#' the copies partition at random between daughters, so replication doubles and
#' division halves, and the expected copy number per cell is unchanged. The
#' neutral value of `fold` is therefore `1`, and `fold = 2` is an aggressive
#' selective sweep. Use [sim_segregate()] to let the level follow from a number
#' of generations and a selection coefficient instead of setting it by hand.
#'
#' Because junction copy number is `multiplicity * copies`, and simulated read
#' support is proportional to junction copy number, the effect on the emitted
#' data is that the **boundary duplication gains read support round after round**
#' while no new junction appears. That is the signature of a clean episome
#' amplifying: one junction, rising `VF`, `PURPLE_JCN` and `PURPLE_CN`, flanks
#' unchanged at the host baseline.
#'
#' Contrast [sim_micronucleation()], which changes the structure (new junctions,
#' oscillating copy number) rather than only the level.
#'
#' @param x An [ecdna_sim], typically straight from [sim_episome()].
#' @param rounds Integer number of replication rounds (default `1`).
#' @param fold Multiplicative gain in copies per round (default `2`, one
#'   doubling). Non-integer values model partial selective sweeps. When the
#'   amplicon is a mixture of circle species each grows by `fold * fitness`, so a
#'   species with a replication advantage expands its share of the population.
#' @param jitter Log-normal dispersion applied to each round's fold change
#'   (default `0`, i.e. deterministic). Set e.g. `0.15` to model uneven
#'   segregation between rounds.
#' @param max_cn Ceiling on peak copy number (default `Inf`; see
#'   [sim_micronucleation()]).
#' @param keep_all Logical; return a list of [ecdna_sim] objects, one per round
#'   including the input as element 1 (default `FALSE`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim], or a list of them when `keep_all = TRUE`.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "REP01", copies = 4)
#' traj <- sim_replicate(ec, rounds = 5, fold = 2, keep_all = TRUE)
#'
#' ## structure is fixed; only the level rises
#' do.call(rbind, lapply(traj, summary))[, c("copies", "max_cn", "n_junctions")]
#' @seealso [sim_episome()], [sim_micronucleation()], [sim_to_epitracer()]
#' @export
sim_replicate <- function(x, rounds = 1L, fold = 2, jitter = 0,
                          max_cn = Inf, keep_all = FALSE, seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), rounds >= 0L, fold > 0)
  if (!is.null(seed)) set.seed(seed)

  out <- vector("list", rounds + 1L)
  out[[1]] <- x
  for (i in seq_len(rounds)) {
    prev <- out[[i]]
    f <- if (jitter > 0) fold * exp(stats::rnorm(1L, 0, jitter)) else fold
    ## every circle species replicates, each at its own rate: a species with
    ## fitness 1 doubles with `fold`, one with fitness > 1 outgrows the rest and
    ## its share of the population -- and so the read support of the junctions
    ## only it carries -- climbs round by round
    pop <- lapply(.sim_population(prev), function(v) {
      v$copies <- v$copies * f * v$fitness; v })
    nxt <- .sim_set_population(prev, pop)
    nxt$copies <- .sim_cap_copies(nxt$copies, nxt$fragments, max_cn)
    if (nxt$copies < sum(vapply(pop, function(v) v$copies, numeric(1)))) {
      ## the ceiling bit: scale every species back by the same factor
      sc <- nxt$copies / sum(vapply(pop, function(v) v$copies, numeric(1)))
      nxt <- .sim_set_population(nxt, lapply(pop, function(v) {
        v$copies <- v$copies * sc; v }))
    }
    nxt$history <- c(prev$history, "replication")
    out[[i + 1L]] <- nxt
  }
  if (keep_all) out else out[[rounds + 1L]]
}

## Remove a genomic interval from a fragment list, preserving amplicon order and
## each fragment's orientation. A forward fragment is traversed in coordinate
## order, an inverted one in reverse, so the pieces left after the cut must be
## re-ordered accordingly.
## NB: interval subtraction is done on coordinates rather than with setdiff().
## On a GRanges, `setdiff` can dispatch to S4Vectors' ELEMENT-WISE set difference
## instead of interval arithmetic (it then returns the fragment untouched, mcols
## and all), depending on which method tables are loaded -- the cause of the
## "multiple methods tables found for 'setdiff'" note this package's dependency
## stack emits. Doing the arithmetic here is unambiguous.
.sim_delete_interval <- function(frags, del) {
  dchr <- as.character(GenomeInfoDb::seqnames(del))[1]
  ds <- GenomicRanges::start(del)[1]; de <- GenomicRanges::end(del)[1]

  pieces <- lapply(seq_along(frags), function(i) {
    f <- frags[i]
    fs <- GenomicRanges::start(f)[1]; fe <- GenomicRanges::end(f)[1]
    same <- as.character(GenomeInfoDb::seqnames(f))[1] == dchr
    if (!same || de < fs || ds > fe) return(f)        # untouched by the deletion

    lo <- if (fs <= ds - 1L) c(fs, ds - 1L) else NULL # piece left of the deletion
    hi <- if (de + 1L <= fe) c(de + 1L, fe) else NULL # piece right of it
    bounds <- Filter(Negate(is.null), list(lo, hi))
    if (!length(bounds)) return(NULL)                 # fragment wholly deleted

    left <- GenomicRanges::GRanges(
      dchr, IRanges::IRanges(vapply(bounds, `[`, numeric(1), 1L),
                             vapply(bounds, `[`, numeric(1), 2L)))
    if (f$inv[1] == 1) left <- rev(left)              # amplicon order is reversed
    left$inv <- f$inv[1]
    left$origin <- f$origin[1]
    left
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(frags[0])
  do.call(c, pieces)
}

## Choose an internal deletion interval of `width` bp that lies inside the circle
## and does not touch any preserved oncogene.
.sim_pick_deletion <- function(frags, preserved, width, tries = 200L) {
  usable <- frags[GenomicRanges::width(frags) > width + 2L]
  if (!length(usable)) return(NULL)
  for (k in seq_len(tries)) {
    i <- if (length(usable) == 1L) 1L else
      sample(length(usable), 1L, prob = GenomicRanges::width(usable))
    f <- usable[i]
    off <- sample.int(GenomicRanges::width(f) - width - 1L, 1L)
    lo <- GenomicRanges::start(f) + off
    d <- GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(f)),
                                IRanges::IRanges(lo, lo + width - 1L))
    if (is.null(preserved) || !length(preserved) ||
        GenomicRanges::countOverlaps(d, preserved) == 0L) return(d)
  }
  NULL
}

#' Grow an ecDNA population by segregation and selection
#'
#' Evolves the amplicon's copy number over `generations` of cell division under
#' an explicit model of how ecDNA actually propagates: each S phase every circle
#' replicates, and at mitosis the `2k` copies partition **at random** between the
#' two daughters (`Binom(2k, 1/2)`), because ecDNA has no centromere. Cells are
#' then sampled in proportion to a fitness that rises with oncogene dosage.
#'
#' Use it in place of [sim_replicate()] when you want the copy level to *follow
#' from* how long the tumour has been evolving and how strongly the amplicon is
#' selected, rather than being set by hand.
#'
#' @section Why division rate does not raise read support:
#' Random segregation is copy-number-neutral in expectation --- replication
#' doubles, division halves, so `E[copies in daughter] = copies in parent`. In a
#' growing tumour with both daughters retained the ecDNA pool grows in step with
#' the cell count, but so does the diploid background it is measured against, and
#' a sequencing library is loaded to a fixed depth rather than a fixed number of
#' cells. Read support is therefore a **ratio** --- roughly
#' `depth * (junction copies / genome copies) * purity` --- and dividing faster
#' leaves it unchanged. With `selection = 0` this function reproduces that: copy
#' number stays flat however many generations are run, while cell-to-cell
#' variance grows and a tail of cells loses the amplicon entirely. Only
#' `selection > 0` raises the population mean, at a rate of about
#' `d log(copies)/d generation ~ selection`. Since bulk WGS sees only that mean,
#' the per-cell spread is invisible to it --- which is why `sim_replicate()`'s
#' deterministic growth is an adequate model for bulk data, and this function's
#' value is in deriving the level rather than asserting it.
#'
#' @param x An [ecdna_sim].
#' @param generations Integer number of cell divisions to run (default `20`).
#' @param selection Dosage-selection coefficient (default `0.05`). `0` is
#'   neutral; the population mean copy number then does not rise, however many
#'   generations elapse.
#' @param plateau Copy number above which extra dosage buys no further advantage
#'   (default `40`), the saturation seen in ecDNA fitness models.
#' @param n_cells Number of cells simulated (default `2000`). Larger is smoother
#'   and slower; the returned copy number is a population mean, so a few thousand
#'   cells is ample.
#' @param addicted Logical; do cells that lose the amplicon entirely die
#'   (default `TRUE`, oncogene addiction)?
#' @param max_cn Ceiling on peak copy number (default `Inf`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim] with `copies` set to the evolved population mean. The
#'   per-cell copy-number distribution is attached as the `"cells"` attribute,
#'   for single-cell or heterogeneity work; bulk rendering ignores it. Species
#'   proportions are rescaled together (co-segregation of distinct circle species
#'   within a cell is not modelled).
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SEG01", copies = 20)
#'
#' ## neutral: twenty generations of division leave the level where it started
#' summary(sim_segregate(ec, generations = 20, selection = 0, seed = 1))$copies
#'
#' ## under selection it climbs
#' summary(sim_segregate(ec, generations = 20, selection = 0.1, seed = 1))$copies
#' @seealso [sim_replicate()], [sim_episome()]
#' @export
sim_segregate <- function(x, generations = 20L, selection = 0.05, plateau = 40,
                          n_cells = 2000L, addicted = TRUE, max_cn = Inf,
                          seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), generations >= 0L, selection >= 0)
  if (!is.null(seed)) set.seed(seed)

  k <- rep(as.integer(round(x$copies)), n_cells)   # copy counts are integers
  for (g in seq_len(generations)) {
    ## S phase then mitosis: 2k copies dealt at random between two daughters
    a <- stats::rbinom(length(k), 2L * k, 0.5)
    kids <- c(a, 2L * k - a)
    if (addicted) {
      alive <- kids > 0
      if (any(alive)) kids <- kids[alive]
    }
    ## dosage-dependent fitness, saturating at `plateau`
    w <- if (selection > 0) (1 + selection)^pmin(kids, plateau) else NULL
    k <- kids[sample.int(length(kids), n_cells, replace = TRUE, prob = w)]
  }

  mean_k <- mean(k)
  if (!is.finite(mean_k) || mean_k <= 0)
    stop("the amplicon was lost from every simulated cell; raise `selection` ",
         "or lower `generations`.", call. = FALSE)

  ## rescale the population to the evolved mean, preserving species proportions
  scale <- .sim_cap_copies(mean_k, x$fragments, max_cn) / x$copies
  x <- .sim_set_population(x, lapply(.sim_population(x), function(v) {
    v$copies <- v$copies * scale; v }))
  x$history <- c(x$history, sprintf("segregate(%dg,s=%.3g)", generations, selection))
  attr(x, "cells") <- k
  x
}

#' Simulate an internal deletion arising within an ecDNA
#'
#' A single ecDNA molecule (or `copies` of them) loses an interval from inside
#' the circle and re-seals. The altered molecule is drawn **out of** the existing
#' pool, so the cell now carries a mixture: intact circles plus deleted ones,
#' each replicating from then on.
#'
#' This is what makes the *order* of events observable. The deleted circle still
#' carries the original boundary duplication, so that junction's copy number is
#' the whole population's; the new internal deletion is carried only by the
#' altered species, so its copy number --- and therefore its read support --- is
#' just that species' abundance. An internal deletion that arose late sits far
#' below the boundary junction it lies inside, which is exactly the read-support
#' stratification [plot_sv_reconstruction()] separates into waves. A deletion
#' introduced at birth and allowed to sweep would instead sit level with the
#' boundary junction, and the two events would be indistinguishable.
#'
#' The deleted interval never overlaps a preserved oncogene: a circle that loses
#' its oncogene confers no advantage and would not persist.
#'
#' @param x An [ecdna_sim].
#' @param start,end Explicit deletion boundaries. When `NULL` (the default) an
#'   interval of `width` bp is drawn at random from inside the circle, avoiding
#'   the preserved loci.
#' @param width Deletion size in bp when drawing at random (default `5e4`).
#' @param copies Number of ecDNA molecules that acquire the deletion (default
#'   `1` --- the minimal event, a single molecule). These are taken out of the
#'   intact pool, so the total copy number is unchanged at the moment it happens.
#' @param fitness Replication rate of the deleted circle relative to the intact
#'   one (default `1`, neutral). This matters more than it looks: a NEUTRAL
#'   minority species keeps a **fixed share** of the population forever, because
#'   both species double together --- one molecule in 32 stays at 3% however long
#'   it replicates, and its junction is emitted ~32-fold below the boundary. Real
#'   amplicons show internal junctions at roughly a half to a twentieth of the
#'   founder's read support, which needs the deleted circle to expand. A value
#'   above `1` does that, and is the biologically expected direction: a circle
#'   that has shed non-essential sequence but kept its oncogene is shorter,
#'   replicates faster and packages more copies per cell.
#' @param label Name for the new species in the population (default
#'   `"internal_del"`).
#' @param homology Repair pathway sealing the deletion (default `"nhej"`).
#' @param seed Optional integer passed to [set.seed()].
#' @return An [ecdna_sim] carrying two circle species. `summary()` reports the
#'   population, and [sim_replicate()] grows both.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "DEL01", copies = 2)
#' ec <- sim_replicate(ec, rounds = 4)          # 32 intact circles
#' ec <- sim_internal_deletion(ec, width = 5e4, fitness = 1.3, seed = 1)
#' ec <- sim_replicate(ec, rounds = 4)          # the shorter circle expands
#' ec                                           # two species, DEL share rising
#' @seealso [sim_replicate()], [sim_episome()], [plot_sv_reconstruction()]
#' @export
sim_internal_deletion <- function(x, start = NULL, end = NULL, width = 5e4,
                                  copies = 1, fitness = 1, label = "internal_del",
                                  homology = c("nhej", "mmej", "nahr"),
                                  seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), copies > 0)
  if (!is.null(seed)) set.seed(seed)
  homology <- match.arg(homology)

  pop <- .sim_population(x)
  ## the deletion arises in the most abundant species present
  src <- which.max(vapply(pop, function(v) v$copies, numeric(1)))
  parent <- pop[[src]]
  if (copies > parent$copies)
    stop("`copies` (", copies, ") exceeds the ", parent$copies,
         " copies of the species the deletion arises in.", call. = FALSE)

  del <- if (!is.null(start) && !is.null(end)) {
    GenomicRanges::GRanges(
      as.character(GenomeInfoDb::seqnames(parent$fragments))[1],
      IRanges::IRanges(start, end))
  } else {
    .sim_pick_deletion(parent$fragments, x$preserved, as.integer(width))
  }
  if (is.null(del))
    stop("could not place an internal deletion of width ", width,
         " bp inside the circle without hitting a preserved locus.", call. = FALSE)
  if (!is.null(x$preserved) && length(x$preserved) &&
      GenomicRanges::countOverlaps(del, x$preserved) > 0L)
    stop("the requested deletion overlaps a preserved oncogene locus.", call. = FALSE)

  frags <- .sim_delete_interval(parent$fragments, del)
  if (!length(frags)) stop("the deletion would remove the entire circle.", call. = FALSE)

  ## the altered molecules leave the parent pool and become their own species
  pop[[src]]$copies <- parent$copies - copies
  pop <- c(pop, list(list(fragments = frags, copies = copies,
                          fitness = fitness, label = label)))

  ## junctions the deleted circle now traverses, merged into the ledger; the
  ## boundary duplication is already there and is left with its original
  ## provenance, so only the new internal deletion is attributed to this step
  new_j <- .sim_junctions_from_fragments(frags, circular = TRUE,
                                         mechanism = "internal_deletion",
                                         round = x$round, homology = homology)
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  known <- do.call(paste, c(x$junctions[, key, with = FALSE], sep = "\r"))
  fresh <- new_j[!do.call(paste, c(new_j[, key, with = FALSE], sep = "\r")) %in% known]

  x$junctions <- data.table::rbindlist(list(x$junctions, fresh),
                                       use.names = TRUE, fill = TRUE)
  x$history <- c(x$history, "internal_deletion")
  .sim_set_population(x, pop)
}

#' Evolve a simulated ecDNA through several rounds
#'
#' Applies [sim_micronucleation()] `rounds` times, so one call produces the whole
#' trajectory from a clean episome to a heavily rearranged amplicon. Returns the
#' final amplicon by default, or the full trajectory (`keep_all = TRUE`) so the
#' rounds can be compared --- the natural way to ask at what point EpiTracer stops
#' calling an amplicon episomal and starts calling it chromothriptic.
#'
#' @param x An [ecdna_sim].
#' @param rounds Integer number of micronucleation + chromothripsis rounds.
#' @param keep_all Logical; return a list of [ecdna_sim] objects, one per round
#'   including the input as element 1 (default `FALSE`, returning only the final
#'   amplicon).
#' @param ... Passed to [sim_micronucleation()] (`n_breaks`, `del_p`, `inv_p`,
#'   `dup_p`, `amplify`, ...). A vector-valued `n_breaks` is recycled across
#'   rounds, so `n_breaks = c(20, 5)` models a severe first shattering followed
#'   by a milder one.
#' @param seed Optional integer passed to [set.seed()] once, before the first
#'   round.
#' @return An [ecdna_sim], or a list of them when `keep_all = TRUE`.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SIM03")
#' traj <- sim_evolve(ec, rounds = 3, n_breaks = 12, keep_all = TRUE, seed = 1)
#' do.call(rbind, lapply(traj, summary))[, c("rounds", "n_junctions", "circle_mb")]
#' @seealso [sim_micronucleation()]
#' @export
sim_evolve <- function(x, rounds = 1L, keep_all = FALSE, ..., seed = NULL) {
  stopifnot(inherits(x, "ecdna_sim"), rounds >= 0L)
  if (!is.null(seed)) set.seed(seed)
  dots <- list(...)
  out <- vector("list", rounds + 1L)
  out[[1]] <- x
  for (i in seq_len(rounds)) {
    args <- lapply(dots, function(v) if (length(v) > 1L) v[[(i - 1L) %% length(v) + 1L]] else v)
    out[[i + 1L]] <- do.call(sim_micronucleation, c(list(x = out[[i]]), args))
  }
  if (keep_all) out else out[[rounds + 1L]]
}
