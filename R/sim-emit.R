## ---------------------------------------------------------------------------
## Emission: render a ground-truth `ecdna_sim` as the inputs EpiTracer consumes.
##
## Both outputs come from the SAME fragment list, so copy number and structural
## variants are consistent by construction:
##
##   copy number  = host baseline - excised allele + copies * coverage(fragments)
##   breakpoints  = the junction ledger, two breakends per junction
##
## On top of that sits a technical model of what a short-read pipeline actually
## reports -- finite tumour purity and depth, copy-number segmentation error and
## resolution limits, junction dropout and false positives, and repair-pathway
## microhomology at the breakpoints. Turning that noise off (`sim_noise(0)`)
## recovers the exact ground truth, which is what the unit tests assert against.
## ---------------------------------------------------------------------------

#' Sequencing and calling noise for simulated amplicons
#'
#' Bundles the technical parameters [sim_to_epitracer()] uses to turn ground
#' truth into the imperfect observation a short-read pipeline reports. Passing a
#' scalar to `sim_noise()` scales every noise term at once, giving a one-knob
#' difficulty sweep; passing named arguments sets terms individually.
#'
#' @param scale Numeric multiplier applied to every stochastic term (default
#'   `1`). `sim_noise(0)` is a noiseless, exact rendering of the ground truth;
#'   `sim_noise(2)` roughly doubles the segmentation error, the read-support
#'   dispersion and the junction dropout and false-positive rates.
#' @param purity Tumour purity (default `0.8`). Lowers junction variant allele
#'   fraction and read support.
#' @param depth Mean sequencing depth in x (default `60`), the usual WGS target.
#'   Read support for a junction scales with `depth * purity`.
#' @param cn_sd Copy-number segmentation error in copies, absolute component
#'   (default `0.25`; PURPLE's segment-level error is typically 0.2-0.5).
#' @param cn_cv Proportional component of the segmentation error (default
#'   `0.02`), so high-copy segments carry proportionally more error.
#' @param vf_cv **Extra-Poisson** dispersion of junction read support (default
#'   `0.1`). Read support is a fragment COUNT, so its dominant noise is counting
#'   noise: support is drawn from a negative binomial whose variance is
#'   `mu + (vf_cv * mu)^2`. A weakly supported junction is therefore noisy in
#'   relative terms (at `mu = 12`, CV about 0.30, Poisson-dominated) while a
#'   strongly supported one is tight (at `mu = 3000`, CV about 0.10, set by
#'   `vf_cv`, which stands for local coverage variation -- GC, mappability).
#'   Modelling this as a constant-CV log-normal instead would over-disperse
#'   high-support junctions by an order of magnitude and hide the fact that
#'   support is proportionate to junction copy number.
#' @param dropout Probability a true junction is missed, before the additional
#'   support-dependent dropout of weakly supported junctions (default `0.03`).
#' @param min_vf Junctions whose simulated read support falls below this are
#'   dropped as unsupported (default `4`), the usual caller floor.
#' @param fp_rate Expected number of false-positive junctions added per amplicon
#'   (default `0.5`).
#' @param min_seg_width Copy-number segments narrower than this are merged into
#'   their neighbours (default `1e4`), modelling the resolution limit of
#'   read-depth segmentation -- without it a heavily shattered circle emits
#'   hundreds of sub-kb segments no real caller would resolve.
#' @param bp_jitter Positional error on reported breakend coordinates, in bp
#'   (default `0`; set e.g. `10` to model imprecise breakpoints).
#' @return A named list of noise parameters, of class `sim_noise`.
#' @examples
#' sim_noise(0)          # exact ground truth
#' sim_noise(purity = 0.4, depth = 30)   # a low-purity, shallow sample
#' @seealso [sim_to_epitracer()]
#' @export
sim_noise <- function(scale = 1, purity = 0.8, depth = 60,
                      cn_sd = 0.25, cn_cv = 0.02, vf_cv = 0.1,
                      dropout = 0.03, min_vf = 4, fp_rate = 0.5,
                      min_seg_width = 1e4, bp_jitter = 0) {
  structure(list(
    purity = purity, depth = depth,
    cn_sd = cn_sd * scale, cn_cv = cn_cv * scale, vf_cv = vf_cv * scale,
    dropout = dropout * scale, min_vf = min_vf, fp_rate = fp_rate * scale,
    min_seg_width = min_seg_width, bp_jitter = bp_jitter * scale
  ), class = "sim_noise")
}

## Microhomology length at a junction, by the repair pathway that sealed it.
## Blunt non-homologous end joining leaves 0-1 bp; polymerase-theta-mediated end
## joining (MMEJ) leaves the 2-10 bp microhomology signature; homology-driven
## repair leaves >= 14 bp. These are the thresholds call_simple_excision() uses
## (mh_min_homology = 2, hr_min_homology = 14) to report the breakpoint-homology
## class of a circularisation junction.
.sim_homlen <- function(pathway, n) {
  vapply(pathway, function(p) switch(p,
    nhej = sample(0:1, 1L, prob = c(0.8, 0.2)),
    mmej = sample(2:10, 1L),
    nahr = sample(14:60, 1L),
    0L), numeric(1), USE.NAMES = FALSE)[seq_len(n)]
}

## Copy-number profile as disjoint segments across the host windows:
## host ploidy, minus one allele where the source chromosome gave up the excised
## interval, plus the amplicon's own contribution (copies x per-circle coverage).
## Allele-specific copy number follows from the episome deriving from ONE
## haplotype: the source allele carries the ecDNA, the other stays at baseline,
## so the amplified segments come out as major = amplicon CN, minor = baseline --
## the allele pattern real ecDNA amplicons show.
.sim_cn_segments <- function(x) {
  host <- x$host
  cov <- .sim_amplicon_cn(x)      # summed over every circle species
  exc <- x$excised

  ## Amplicon or excision may sit outside the stored host window (a chimera pulls
  ## in a second chromosome); widen the background so every locus has a baseline.
  host <- .sim_widen_host(host, c(cov, exc))

  pieces <- Filter(function(g) !is.null(g) && length(g),
                   list(.sim_bare(host), .sim_bare(cov), .sim_bare(exc)))
  seg <- GenomicRanges::disjoin(do.call(c, pieces))
  seg <- IRanges::subsetByOverlaps(seg, host)

  ## Per-segment contributions, read at each segment's midpoint so a segment is
  ## wholly inside or outside each source range (guaranteed by disjoin()).
  mid <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(seg),
                                IRanges::IRanges(GenomicRanges::start(seg),
                                                 width = 1L))
  n_amp <- .sim_lookup(mid, cov, "n", 0)      # already abundance-weighted
  is_exc <- GenomicRanges::countOverlaps(mid, exc %||% GenomicRanges::GRanges()) > 0
  pl <- .sim_lookup(mid, host, "ploidy", 2)

  ## source allele: baseline, less the excised copy, plus the ecDNA
  a_src <- pmax(0, pl / 2 - as.numeric(is_exc)) + n_amp
  a_alt <- rep(pl / 2, length(seg))

  data.table::data.table(
    seqnames = as.character(GenomeInfoDb::seqnames(seg)),
    start = GenomicRanges::start(seg), end = GenomicRanges::end(seg),
    sample = x$sample,
    copyNumber = a_src + a_alt,
    ploidy = pl,
    majorAlleleCopyNumber = pmax(a_src, a_alt),
    minorAlleleCopyNumber = pmin(a_src, a_alt)
  )
}

.sim_bare <- function(g) {
  if (is.null(g) || !length(g)) return(NULL)
  S4Vectors::mcols(g) <- NULL
  g
}

## Extend the host background to cover any range not already inside it, at the
## prevailing baseline ploidy, so a chimeric amplicon's second chromosome gets a
## flank to be compared against.
## (Built by union rather than setdiff: see the note on .sim_delete_interval --
## `setdiff` on a GRanges can silently dispatch to element-wise set difference.)
.sim_widen_host <- function(host, extra) {
  extra <- Filter(function(g) !is.null(g) && length(g), list(extra))
  if (!length(extra)) return(host)
  extra <- GenomicRanges::reduce(do.call(c, lapply(extra, .sim_bare))) + 1.5e7
  GenomicRanges::start(extra) <- pmax(1, GenomicRanges::start(extra))

  merged <- GenomicRanges::reduce(c(.sim_bare(host), extra))
  ## carry each chromosome's own baseline across (a chimera can pair a diploid
  ## chromosome with a polysomic one), falling back to the modal host ploidy
  mid <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(merged),
                                IRanges::IRanges(GenomicRanges::start(merged), width = 1L))
  merged$ploidy <- .sim_lookup(mid, host, "ploidy", stats::median(host$ploidy))
  GenomicRanges::sort(merged)
}

## Value of metadata column `col` of `y` at each position of `x`, `default`
## where there is no overlap.
.sim_lookup <- function(x, y, col, default) {
  if (is.null(y) || !length(y)) return(rep(default, length(x)))
  hit <- GenomicRanges::findOverlaps(x, y, select = "first")
  v <- as.numeric(S4Vectors::mcols(y)[[col]])[hit]
  v[is.na(v)] <- default
  v
}

## Merge copy-number segments below the caller's resolution into the neighbour
## they most resemble, then fuse runs of equal copy number. Models the fact that
## read-depth segmentation cannot resolve arbitrarily narrow segments.
.sim_merge_segments <- function(cn, min_width) {
  cn <- cn[order(cn$seqnames, cn$start)]
  out <- lapply(split(cn, cn$seqnames), function(d) {
    repeat {
      w <- d$end - d$start + 1
      i <- which(w < min_width)
      if (!length(i) || nrow(d) < 2L) break
      k <- i[which.min(w[i])]
      ## absorb into whichever neighbour it is closer to in copy number
      nb <- setdiff(c(k - 1L, k + 1L), c(0L, nrow(d) + 1L))
      j <- nb[which.min(abs(d$copyNumber[nb] - d$copyNumber[k]))]
      d$start[j] <- min(d$start[j], d$start[k]); d$end[j] <- max(d$end[j], d$end[k])
      d <- d[-k]
    }
    ## fuse adjacent segments of equal (rounded) copy number
    r <- rle(round(d$copyNumber, 3))
    grp <- rep(seq_along(r$lengths), r$lengths)
    d[, .(start = min(start), end = max(end), sample = sample[1],
          copyNumber = stats::weighted.mean(copyNumber, end - start + 1),
          ploidy = ploidy[1],
          majorAlleleCopyNumber = stats::weighted.mean(majorAlleleCopyNumber, end - start + 1),
          minorAlleleCopyNumber = stats::weighted.mean(minorAlleleCopyNumber, end - start + 1)),
      by = .(seqnames, grp)][, grp := NULL]
  })
  data.table::rbindlist(out)
}

#' Render a simulated amplicon as EpiTracer caller inputs
#'
#' Turns an [ecdna_sim] into the four GRanges [call_simple_excision()] and the
#' mechanism callers take, applying a model of what a short-read pipeline would
#' actually report: finite purity and depth, copy-number segmentation error and
#' resolution, junction dropout and false positives, and repair-pathway
#' microhomology. With `noise = sim_noise(0)` the rendering is exact.
#'
#' @param x An [ecdna_sim], or a list of them (one per sample) to render as one
#'   cohort.
#' @param noise A [sim_noise()] specification, or a numeric scale passed to it.
#' @param min_cn_ratio Amplification threshold used to define the reported
#'   amplicon footprint, matching [detect_amplicon_seeds()] (default `3`).
#' @param seed Optional integer passed to [set.seed()].
#' @return A named list of four [GenomicRanges::GRanges] --- `ecdna_gr`,
#'   `breakpoints_gr`, `cnv_gr`, `cancer_genes_gr` --- plus `truth`, a
#'   [data.table::data.table] of the ground-truth mechanism labels per sample.
#'   The first four can be passed straight to any EpiTracer caller.
#' @examples
#' ec <- sim_episome(seed_locus("EGFR"), sample = "SIM01", copies = 50)
#' inp <- sim_to_epitracer(ec, seed = 1)
#' res <- call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
#'                             inp$cnv_gr, inp$cancer_genes_gr)
#' unique(res$episomal)
#' @seealso [sim_episome()], [sim_to_plot_inputs()], [sim_cohort()]
#' @export
sim_to_epitracer <- function(x, noise = sim_noise(), min_cn_ratio = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.numeric(noise)) noise <- sim_noise(noise)
  if (inherits(x, "ecdna_sim")) x <- list(x)
  stopifnot(all(vapply(x, inherits, logical(1), "ecdna_sim")))

  parts <- lapply(x, .sim_emit_one, noise = noise, min_cn_ratio = min_cn_ratio)

  list(
    ecdna_gr        = to_granges(data.table::rbindlist(lapply(parts, `[[`, "ecdna"))),
    breakpoints_gr  = to_granges(data.table::rbindlist(lapply(parts, `[[`, "bp"), fill = TRUE)),
    cnv_gr          = to_granges(data.table::rbindlist(lapply(parts, `[[`, "cnv"))),
    cancer_genes_gr = GenomicRanges::reduce(
      do.call(c, lapply(parts, `[[`, "genes")), with.revmap = FALSE,
      ignore.strand = TRUE) %>% .sim_relabel_genes(lapply(parts, `[[`, "genes")),
    truth           = data.table::rbindlist(lapply(parts, `[[`, "truth"), fill = TRUE)
  )
}

## Keep a `gene` label on the reduced oncogene set (reduce() drops metadata).
.sim_relabel_genes <- function(reduced, originals) {
  orig <- do.call(c, originals)
  hit <- GenomicRanges::findOverlaps(reduced, orig, select = "first")
  reduced$gene <- as.character(orig$gene)[hit]
  reduced$gene[is.na(reduced$gene)] <- "GENE"
  reduced
}

.sim_emit_one <- function(x, noise, min_cn_ratio) {
  ## ---- copy number -------------------------------------------------------
  cn <- .sim_cn_segments(x)
  cn <- .sim_merge_segments(cn, noise$min_seg_width)
  clean_cn <- data.table::copy(cn)                       # noiseless, for CN lookup
  if (noise$cn_sd > 0 || noise$cn_cv > 0) {
    err <- stats::rnorm(nrow(cn), 0, noise$cn_sd + noise$cn_cv * cn$copyNumber)
    cn$copyNumber <- pmax(0, cn$copyNumber + err)
    ## keep the allele split consistent with the perturbed total
    frac <- data.table::fifelse(clean_cn$copyNumber > 0,
                                clean_cn$minorAlleleCopyNumber / clean_cn$copyNumber, 0)
    cn$minorAlleleCopyNumber <- pmax(0, cn$copyNumber * frac)
    cn$majorAlleleCopyNumber <- pmax(0, cn$copyNumber - cn$minorAlleleCopyNumber)
  }

  ## ---- amplicon footprint ------------------------------------------------
  amp <- clean_cn[copyNumber > min_cn_ratio * ploidy]
  amp_gr <- if (nrow(amp))
    GenomicRanges::reduce(to_granges(amp), min.gapwidth = 2e6) else GenomicRanges::GRanges()
  ecdna <- if (length(amp_gr)) data.table::data.table(
    seqnames = as.character(GenomeInfoDb::seqnames(amp_gr)),
    start = GenomicRanges::start(amp_gr), end = GenomicRanges::end(amp_gr),
    ID = paste0(x$sample, "_amp", seq_along(amp_gr)), WGS_ID = x$sample
  ) else data.table::data.table(
    seqnames = character(), start = numeric(), end = numeric(),
    ID = character(), WGS_ID = character())

  ## ---- breakpoints -------------------------------------------------------
  bp <- .sim_emit_breakpoints(x, clean_cn, noise)

  ## ---- oncogenes and truth ----------------------------------------------
  genes <- x$preserved
  if (is.null(genes$gene)) genes$gene <- paste0("GENE", seq_along(genes))
  s <- summary(x)

  list(cnv = cn, ecdna = ecdna, bp = bp, genes = genes,
       truth = data.table::data.table(
         WGS_ID = x$sample, ID = if (nrow(ecdna)) ecdna$ID[1] else NA_character_,
         mechanism = utils::tail(x$history, 1), history = s$history,
         rounds = x$round, chimeric = s$n_chr > 1L,
         episomal = TRUE,          # every amplicon here is circular by construction
         truth_copies = x$copies, truth_max_cn = s$max_cn,
         truth_n_junctions = s$n_junctions, truth_circle_mb = s$circle_mb))
}

## Two breakend rows per junction, carrying the PURPLE-style columns the callers
## require. Junction copy number is the circle's copy number times the number of
## times the junction is traversed per circle; read support follows from that,
## the tumour purity and the depth.
.sim_emit_breakpoints <- function(x, cn, noise) {
  j <- data.table::as.data.table(x$junctions)
  if (!nrow(j)) return(.sim_empty_bp(x$sample))

  ## Junction copy number: a host excision scar sits on a single chromosomal
  ## allele; an amplicon junction is carried only by the circle species that
  ## traverse it, so its copy number is summed over those species. A junction
  ## that arose late therefore sits below one present since birth.
  key <- c("chrom1", "start1", "strand1", "chrom2", "start2", "strand2")
  pj <- .sim_junction_jcn(x)
  jcn <- if (nrow(pj)) {
    m <- merge(j[, key, with = FALSE], pj, by = key, all.x = TRUE, sort = FALSE)
    ## preserve j's row order, which merge() does not guarantee
    m <- m[match(do.call(paste, c(j[, key, with = FALSE], sep = "\r")),
                 do.call(paste, c(m[, key, with = FALSE], sep = "\r")))]
    m$jcn
  } else rep(NA_real_, nrow(j))
  jcn[is.na(jcn)] <- (j$multiplicity * x$copies)[is.na(jcn)]
  jcn <- data.table::fifelse(j$origin == "host", 1, jcn)

  ## Expected supporting fragments. A junction is spanned by read pairs in
  ## proportion to its copy number relative to the sample's average, scaled by the
  ## haploid depth and tumour purity.
  hap_depth <- noise$depth / mean(cn$ploidy)
  vf_mu <- jcn * noise$purity * hap_depth * 0.5

  ## Support is a fragment count: draw it from a negative binomial with mean
  ## vf_mu and variance vf_mu + (vf_cv * vf_mu)^2 -- Poisson counting noise plus
  ## an extra-Poisson term for local coverage variation. Expected support stays
  ## strictly proportionate to junction copy number, so a junction carried by few
  ## circles lands below one carried by many, by exactly their copy-number ratio.
  vf <- if (noise$vf_cv > 0)
    stats::rnbinom(length(vf_mu), mu = vf_mu, size = 1 / noise$vf_cv^2)
  else round(vf_mu)
  vf <- pmax(0, vf)

  ## Variant allele fraction of the junction against the local tumour+normal mix.
  cn_at <- .sim_cn_at(j$chrom1, j$start1, cn)
  af <- pmin(1, jcn * noise$purity /
               pmax(1e-6, cn_at * noise$purity + mean(cn$ploidy) * (1 - noise$purity)))

  j[, `:=`(jcn = jcn, vf = vf, af = af)]

  ## Detection: unsupported junctions and a random dropout are lost. The boundary
  ## junction of a real episome is very highly supported and effectively never
  ## dropped, which is what makes the mechanism callable.
  keep <- j$vf >= noise$min_vf & stats::runif(nrow(j)) >= noise$dropout
  j <- j[keep]
  if (!nrow(j)) return(.sim_empty_bp(x$sample))

  ## False-positive junctions: low-support artefacts scattered in the amplicon.
  n_fp <- if (noise$fp_rate > 0) stats::rpois(1L, noise$fp_rate) else 0L
  fp <- if (n_fp > 0) .sim_false_positives(n_fp, x, cn, noise) else NULL

  j <- data.table::rbindlist(list(j, fp), use.names = TRUE, fill = TRUE)
  j[, event := sprintf("%s_SV%03d", x$sample, .I)]

  ## Melt each junction into its two breakends.
  bp <- data.table::rbindlist(list(
    j[, .(seqnames = chrom1, start = start1, bp_strand = strand1,
          event, svclass, jcn, vf, af, homology)],
    j[, .(seqnames = chrom2, start = start2, bp_strand = strand2,
          event, svclass, jcn, vf, af, homology)]
  ))
  if (noise$bp_jitter > 0)
    bp[, start := pmax(1, round(start + stats::rnorm(.N, 0, noise$bp_jitter)))]

  hom <- .sim_homlen(bp$homology, nrow(bp))
  data.table::data.table(
    seqnames = bp$seqnames, start = bp$start, end = bp$start,
    WGS_ID = x$sample, event = bp$event, svclass = bp$svclass,
    bp_strand = bp$bp_strand,
    PURPLE_AF = bp$af, PURPLE_JCN = bp$jcn, VF = as.numeric(bp$vf),
    PURPLE_CN = .sim_cn_at(bp$seqnames, bp$start, cn),
    insLen = as.integer(data.table::fifelse(bp$homology == "nhej",
                                            stats::rpois(nrow(bp), 0.5), 0L)),
    HOMLEN = as.integer(hom)
  )[order(seqnames, start)]
}

.sim_empty_bp <- function(sample) {
  data.table::data.table(
    seqnames = character(), start = numeric(), end = numeric(),
    WGS_ID = character(), event = character(), svclass = character(),
    bp_strand = character(), PURPLE_AF = numeric(), PURPLE_JCN = numeric(),
    VF = numeric(), PURPLE_CN = numeric(), insLen = integer(), HOMLEN = integer())
}

## Total copy number of the segment containing each position.
.sim_cn_at <- function(chrom, pos, cn) {
  q <- GenomicRanges::GRanges(chrom, IRanges::IRanges(pmax(1, pos), width = 1L))
  hit <- GenomicRanges::findOverlaps(q, to_granges(cn), select = "first")
  v <- cn$copyNumber[hit]
  v[is.na(v)] <- mean(cn$ploidy)
  as.numeric(v)
}

## Artefactual junctions: both breakends inside the amplified footprint (where
## mapping artefacts concentrate), at low read support.
.sim_false_positives <- function(n, x, cn, noise) {
  amp <- cn[copyNumber > 3 * ploidy]
  if (!nrow(amp)) return(NULL)
  pick <- function() {
    i <- sample(nrow(amp), 1L, prob = amp$end - amp$start + 1)
    list(chr = amp$seqnames[i],
         pos = amp$start[i] + sample.int(amp$end[i] - amp$start[i] + 1L, 1L) - 1L)
  }
  rows <- lapply(seq_len(n), function(k) {
    a <- pick(); b <- pick()
    s1 <- sample(c("+", "-"), 1L); s2 <- sample(c("+", "-"), 1L)
    jj <- .sim_junction(list(chrom = a$chr, pos = a$pos, strand = s1),
                        list(chrom = b$chr, pos = b$pos, strand = s2),
                        mechanism = "artefact", round = NA_integer_,
                        homology = "nhej", origin = "artefact")
    jj[, `:=`(jcn = 1, vf = round(noise$min_vf + stats::rexp(1L, 0.2)),
              af = 0.02)]
    jj
  })
  data.table::rbindlist(rows)
}

#' Render a simulated amplicon as plotter inputs
#'
#' The same simulated amplicon in the data-frame format [plot_sv_linear()] and
#' [plot_sv_reconstruction()] expect, so a simulated trajectory can be inspected
#' with the package's own viewers --- the quickest way to see whether a parameter
#' set reproduces the look of a patient amplicon.
#'
#' @param x An [ecdna_sim], or a list of them.
#' @param noise A [sim_noise()] specification, or a numeric scale.
#' @param seed Optional integer passed to [set.seed()].
#' @return A named list of three `data.frame`s --- `cnv_data`, `sv_data`,
#'   `wgd_data` --- matching [plot_sv_linear()]'s arguments.
#' @examples
#' ec <- sim_evolve(sim_episome(seed_locus("EGFR"), sample = "SIM04"),
#'                  rounds = 2, n_breaks = 15, seed = 1)
#' pin <- sim_to_plot_inputs(ec, seed = 1)
#' p <- plot_sv_linear(sample = "SIM04", cnv_data = pin$cnv_data,
#'                     sv_data = pin$sv_data, wgd_data = pin$wgd_data,
#'                     chromosome = "chr7")
#' @seealso [sim_to_epitracer()], [plot_sv_linear()], [plot_sv_reconstruction()]
#' @export
sim_to_plot_inputs <- function(x, noise = sim_noise(), seed = NULL) {
  inp <- sim_to_epitracer(x, noise = noise, seed = seed)
  cnv <- gr2dt(inp$cnv_gr)

  ## Breakends back to paired BEDPE: the two rows sharing an `event`.
  bp <- gr2dt(inp$breakpoints_gr)[order(event, seqnames, start)]
  sv <- bp[, if (.N == 2L) .(
    chrom1 = seqnames[1], start1 = start[1], chrom2 = seqnames[2], start2 = start[2],
    strand1 = bp_strand[1], strand2 = bp_strand[2], svclass = svclass[1],
    VF = VF[1], JCN = PURPLE_JCN[1], sample = WGS_ID[1]), by = event]
  sv[, event := NULL]

  list(
    cnv_data = as.data.frame(cnv[, .(sample, seqnames, start, end, copyNumber,
                                     ploidy, majorAlleleCopyNumber,
                                     minorAlleleCopyNumber)]),
    sv_data  = as.data.frame(sv),
    wgd_data = data.frame(sample = unique(cnv$sample), Polyploidy = "No",
                          stringsAsFactors = FALSE)
  )
}
