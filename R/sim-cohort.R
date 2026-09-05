## ---------------------------------------------------------------------------
## Cohort layer: many simulated samples in one call, with ground-truth labels.
##
## sim_cohort() draws a cohort from a MIXTURE of trajectories -- amplicons caught
## at origin, amplicons that have been through one or several rounds of
## micronucleation and chromothripsis, and chimeras built from two co-encapsulated
## episomes -- randomising the oncogene, the copy level, the shattering intensity
## and the host ploidy across samples. The result is a single set of caller
## inputs plus a truth table, so a whole cohort can be pushed through EpiTracer in
## one call and scored (sim_benchmark()).
##
## Every sample is an independent WGS_ID, so the callers process them in isolation
## exactly as they would a real cohort.
## ---------------------------------------------------------------------------

## The subset of EpiTracer's bundled oncogene panel sim_cohort() draws seed loci
## from by default: recurrently ecDNA-amplified oncogenes spread across several
## chromosomes, so simulated chimeras join non-homologous partners.
.sim_default_genes <- c("EGFR", "CDK4", "MDM2", "MYC", "TERT", "SOX2", "PDGFRA", "MET")

#' Simulate a cohort of ecDNA amplicons with known mechanisms
#'
#' Draws `n` independent samples from a mixture of ecDNA trajectories and renders
#' them as one set of EpiTracer inputs with a ground-truth label per sample. The
#' three trajectory classes correspond to the questions the suite was built to
#' answer: what a **freshly born** episome looks like, what it looks like after
#' **one or several rounds** of micronucleation and chromothripsis, and what a
#' **chimera** of two co-encapsulated episomes looks like.
#'
#' @param n Number of samples (default `30`).
#' @param mixture Named numeric vector of class weights, normalised internally.
#'   Recognised names are `episomal` (born and unshattered), `chromothriptic`
#'   (one or more micronucleation rounds) and `chimeric` (two episomes fused).
#'   Defaults to an even split.
#' @param genes Character vector of oncogene symbols to draw seed loci from
#'   (default: a panel of recurrently ecDNA-amplified oncogenes).
#' @param rounds Integer vector of micronucleation rounds to sample from for the
#'   `chromothriptic` class (default `1:3`).
#' @param n_breaks Integer vector of chromothriptic breakpoint counts to sample
#'   from (default `8:30`).
#' @param copies Numeric length-2 range of ecDNA copies per cell at birth
#'   (default `c(15, 80)`).
#' @param host_ploidy Numeric vector of host baseline ploidies to sample from
#'   (default `c(2, 2, 2, 3)`, so a quarter of samples sit on a polysomic
#'   chromosome -- the case that defeats a ploidy-based flank test).
#' @param homology Repair pathways to sample the circularisation junction from
#'   (default `c("nhej", "nhej", "mmej")`).
#' @param noise A [sim_noise()] specification, or a numeric scale.
#' @param prefix Sample-name prefix (default `"SIM"`).
#' @param seed Optional integer passed to [set.seed()].
#' @param mc.cores Cores for [parallel::mclapply()] (default `1`).
#' @return A named list: the four caller inputs (`ecdna_gr`, `breakpoints_gr`,
#'   `cnv_gr`, `cancer_genes_gr`), a `truth` [data.table::data.table] with one
#'   row per simulated amplicon, and `sims`, the list of underlying [ecdna_sim]
#'   objects (so any sample can be re-rendered at different noise, or plotted).
#' @examples
#' co <- sim_cohort(n = 6, seed = 1)
#' co$truth[, .(WGS_ID, class, rounds, chimeric, truth_max_cn)]
#' @seealso [sim_episome()], [sim_evolve()], [sim_fuse_episomes()], [sim_benchmark()]
#' @export
sim_cohort <- function(n = 30,
                       mixture = c(episomal = 1, chromothriptic = 1, chimeric = 1),
                       genes = .sim_default_genes,
                       rounds = 1:3, n_breaks = 8:30, copies = c(15, 80),
                       host_ploidy = c(2, 2, 2, 3),
                       homology = c("nhej", "nhej", "mmej"),
                       noise = sim_noise(), prefix = "SIM",
                       seed = NULL, mc.cores = 1L) {
  if (!is.null(seed)) set.seed(seed)
  if (is.numeric(noise) && !inherits(noise, "sim_noise")) noise <- sim_noise(noise)
  mixture <- mixture[mixture > 0]
  stopifnot(length(mixture) > 0, all(names(mixture) %in%
    c("episomal", "chromothriptic", "chimeric")))

  cls <- sample(names(mixture), n, replace = TRUE, prob = mixture / sum(mixture))
  ids <- sprintf("%s%04d", prefix, seq_len(n))

  sims <- parallel::mclapply(seq_len(n), function(i) {
    .sim_one_sample(ids[i], cls[i], genes, rounds, n_breaks, copies,
                    host_ploidy, homology)
  }, mc.cores = mc.cores)

  inp <- sim_to_epitracer(sims, noise = noise)
  inp$truth[, class := cls[match(WGS_ID, ids)]]
  inp$sims <- sims
  inp
}

## One sample of the requested trajectory class.
.sim_one_sample <- function(id, class, genes, rounds, n_breaks, copies,
                            host_ploidy, homology) {
  pick_gene <- function(exclude = character()) sample(setdiff(genes, exclude), 1L)
  cp <- function() stats::runif(1L, copies[1], copies[2])
  hp <- sample(host_ploidy, 1L)
  hm <- sample(homology, 1L)

  born <- function(g, ...) sim_episome(seed_locus(g), sample = id, copies = cp(),
                                       host_ploidy = hp, homology = hm, ...)
  switch(class,
    episomal = born(pick_gene()),
    chromothriptic = sim_evolve(born(pick_gene()),
                                rounds = sample(rounds, 1L),
                                n_breaks = sample(n_breaks, 1L)),
    chimeric = {
      g1 <- pick_gene(); g2 <- pick_gene(exclude = g1)
      sim_fuse_episomes(born(g1), born(g2), n_breaks = sample(n_breaks, 1L))
    },
    stop("unknown simulation class: ", class, call. = FALSE))
}

#' Score EpiTracer's mechanism callers against simulated ground truth
#'
#' Runs the mechanism callers over a simulated cohort and joins their per-amplicon
#' verdicts to the known truth, returning both the per-sample calls and a summary
#' of how each trajectory class was classified. This is how the simulation suite
#' is meant to be consumed: build a cohort whose mechanisms are known by
#' construction, then ask what EpiTracer says about it.
#'
#' Note the classes are **not** mutually exclusive by design --- a chimeric
#' amplicon has genuinely been through micronucleation *and* chromothripsis, so a
#' `TRUE` from both callers is correct, not a false positive. The summary
#' therefore reports the rate at which each caller fires per class rather than a
#' single confusion matrix.
#'
#' @param cohort The list returned by [sim_cohort()].
#' @param callers Character vector naming which callers to run; any of
#'   `"simple_excision"`, `"chromothripsis"`, `"chimeric"`, `"brf"`,
#'   `"bfb"`. Defaults to the first three, the mechanisms this suite generates.
#' @param mc.cores Cores passed to the callers (default `1`).
#' @return A list with `calls` (one row per amplicon: truth columns plus one
#'   logical per caller) and `summary` (per-class firing rate of each caller).
#' @examples
#' co <- sim_cohort(n = 6, seed = 1)
#' bm <- sim_benchmark(co)
#' bm$summary
#' @seealso [sim_cohort()]
#' @export
sim_benchmark <- function(cohort,
                          callers = c("simple_excision", "chromothripsis",
                                      "chimeric"),
                          mc.cores = 1L) {
  stopifnot(is.list(cohort), !is.null(cohort$truth))
  args <- list(cohort$ecdna_gr, cohort$breakpoints_gr, cohort$cnv_gr,
               cohort$cancer_genes_gr)
  fns <- list(simple_excision = list(call_simple_excision, "episomal"),
              chromothripsis  = list(call_chromothripsis,  "chromothripsis"),
              chimeric        = list(call_chimeric_amplicon, "chimeric"),
              brf             = list(call_brf,             "brf"),
              bfb             = list(call_bfb,             "bfb"))
  callers <- match.arg(callers, names(fns), several.ok = TRUE)

  calls <- data.table::copy(cohort$truth)
  for (nm in callers) {
    f <- fns[[nm]][[1]]; col <- fns[[nm]][[2]]
    res <- suppressWarnings(suppressMessages(
      do.call(f, c(args, list(mc.cores = mc.cores)))))
    ## an amplicon fires if ANY of its breakpoints carries the positive flag
    per <- if (nrow(res)) res[, .(hit = any(get(col) == "TRUE")), by = ID] else
      data.table::data.table(ID = character(), hit = logical())
    calls[per, (nm) := i.hit, on = "ID"]
    calls[is.na(get(nm)), (nm) := FALSE]
  }

  summ <- calls[, c(list(n = .N), lapply(.SD, function(v) mean(as.logical(v)))),
                by = class, .SDcols = callers]
  list(calls = calls[], summary = summ[])
}
