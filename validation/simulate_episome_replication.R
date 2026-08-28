## ---------------------------------------------------------------------------
## EpiTracer validation -- the simplest ecDNA case, step by step
##
## One episomal amplicon, born by simple excision, carrying NOTHING but its
## boundary duplication (the circularisation junction) plus the excision scar it
## left behind on the chromosome. It then replicates; after four rounds a single
## molecule acquires an internal deletion, and replication continues.
##
##   round 0     birth by simple excision      -> 1 DUP + 1 scar DEL
##   rounds 1-4  replication                   -> copies double each round
##   round 4     one molecule loses an interval-> a second circle species
##   rounds 5-10 replication to a plateau       -> the shorter circle takes share
##
## What this shows:
##   * the junction set never grows by replication alone -- only the level moves
##   * boundary read support (VF) is proportional to ecDNA copies per cell
##   * the amplicon is invisible until copy number clears min_cn_ratio * ploidy
##   * a LATE internal deletion is carried by only a fraction of the circles, so
##     it is emitted below the boundary junction it lies inside. That gap is the
##     timing signal: a deletion present from birth would sit level with the
##     boundary, and the two events would be indistinguishable.
##   * how DEEP the gap is depends on whether the deleted circle expands. Under
##     neutral growth its share is frozen and the gap stays at its founding ratio;
##     given a replication advantage the deleted subclone climbs and the gap
##     narrows into the range real amplicons show.
##
## Run from the package root:  Rscript validation/simulate_episome_replication.R
## Writes validation/output/{episome_replication.csv, episome_replication.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2); library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
q <- function(expr) suppressWarnings(suppressMessages(expr))

GENE          <- "EGFR"
START_COPIES  <- 2      # a single excised circle, on both sister chromatids
ROUNDS_PRE    <- 4      # replication before the deletion
ROUNDS_POST   <- 6      # replication after it
MAX_CN        <- 250    # ecDNA copy number plateaus under selection
FOLD          <- 2      # one doubling per round
DEL_START     <- 55300000
DEL_END       <- 55350000
DEL_COPIES    <- 1      # ONE molecule acquires it -- the minimal event
## The shortened circle keeps its oncogene but has less DNA to replicate, so it
## outgrows the intact one. This matters: a NEUTRAL minority species keeps a
## fixed share of the population forever (one circle in 32 stays at 3%), which
## would put the internal deletion ~32-fold below the boundary junction. Real
## amplicons show internal junctions at a half to a twentieth of the founder's
## read support, which needs the deleted circle to expand. 1.3 lands at ~8-fold,
## matching the DO11501T1 CDK4 amplicon (founder VF 496 vs bulk ~75).
DEL_FITNESS   <- 1.3
PURITY        <- 0.8
DEPTH         <- 60
MIN_CN_RATIO  <- 3

## Clean rendering: no artefact junctions, no dropout, so the emitted breakpoints
## are exactly the circle's own. Copy-number and read-support noise stay ON, so
## the numbers still look like a pipeline's output rather than arithmetic.
CLEAN <- sim_noise(purity = PURITY, depth = DEPTH, fp_rate = 0, dropout = 0)

theme_val <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background   = element_rect(fill = "white", colour = NA),
          panel.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 11))
}

## ---------------------------------------------------------------------------
## 1. Build the trajectory
## ---------------------------------------------------------------------------
ec <- sim_episome(seed_locus(GENE), sample = "EPI01", copies = START_COPIES,
                  host_ploidy = 2, seed = 1)
message("At birth:"); print(ec)
stopifnot(nrow(ec$junctions) == 2L, sum(ec$junctions$svclass == "DUP") == 1L)

traj <- sim_replicate(ec, rounds = ROUNDS_PRE, fold = FOLD, keep_all = TRUE)

## the deletion arises in one molecule at round ROUNDS_PRE, then both grow
after_del <- sim_internal_deletion(traj[[length(traj)]],
                                   start = DEL_START, end = DEL_END,
                                   copies = DEL_COPIES, fitness = DEL_FITNESS)
message("\nAfter the internal deletion:"); print(after_del)

traj <- c(traj, sim_replicate(after_del, rounds = ROUNDS_POST, fold = FOLD,
                              max_cn = MAX_CN, keep_all = TRUE))
## traj is now: rounds 0..4 (intact), then round 4+del, then rounds 5..10

## ---------------------------------------------------------------------------
## 2. Render each step as WGS-like inputs and call it
## ---------------------------------------------------------------------------
res <- rbindlist(lapply(seq_along(traj), function(i) {
  x <- traj[[i]]
  has_del <- "internal_deletion" %in% x$history
  step <- if (i <= ROUNDS_PRE + 1L) as.character(i - 1L)
          else if (i == ROUNDS_PRE + 2L) paste0(ROUNDS_PRE, "+del")
          else as.character(i - 2L)

  inp <- sim_to_epitracer(x, noise = CLEAN, seed = 100 + i)
  bp  <- gr2dt(inp$breakpoints_gr)
  dup <- bp[svclass == "DUP"][1]
  idel <- bp[svclass == "DEL" & start > DEL_START - 1e5 & start < DEL_END + 1e5][1]

  detected <- length(inp$ecdna_gr) > 0L
  se <- if (detected) q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                                             inp$cnv_gr, inp$cancer_genes_gr)) else NULL
  data.table(
    step, order = i, has_del,
    copies       = round(x$copies, 1),
    n_species    = length(EpiTracer:::.sim_population(x)),
    del_share    = { p <- EpiTracer:::.sim_population(x)
                     if (length(p) > 1L) round(100 * p[[2]]$copies / x$copies, 1)
                     else NA_real_ },
    n_junctions  = length(unique(bp$event)),
    amplicon_cn  = round(max(inp$cnv_gr$copyNumber), 1),
    dup_jcn = dup$PURPLE_JCN, dup_vf = dup$VF,
    del_jcn = if (!is.na(idel$start)) idel$PURPLE_JCN else NA_real_,
    del_vf  = if (!is.na(idel$start)) idel$VF else NA_real_,
    detected,
    episomal = !is.null(se) && nrow(se) > 0 && any(se$episomal == "TRUE"))
}))

fwrite(res, "validation/output/episome_replication.csv")
cat("\n=== Clean episome: replication, then a late internal deletion ===\n")
print(res[, .(step, copies, n_junctions, amplicon_cn, del_share,
              dup_jcn, dup_vf, del_jcn, del_vf, episomal)])

fc <- res[episomal == TRUE][1]
cat(sprintf("\nFirst called episomal at round %s (%.0f copies, CN %.1f).\n",
            fc$step, fc$copies, fc$amplicon_cn))
last <- res[.N]
cat(sprintf("At the end the boundary DUP sits at VF %.0f and the internal DEL at VF %.0f",
            last$dup_vf, last$del_vf))
cat(sprintf(" -- a %.0f-fold gap that dates the deletion to after the amplicon was established.\n",
            last$dup_vf / last$del_vf))
cat(sprintf("The deleted circle (fitness %.2f) has reached %.1f%% of the ecDNA population.\n",
            DEL_FITNESS, last$del_share))

## ---------------------------------------------------------------------------
## 3. Figures
## ---------------------------------------------------------------------------
ACCENT <- "#0f7d78"; WARN <- "#c2571a"; DELC <- "#3f5d9e"
res[, xo := seq_len(.N)]
lab <- res$step

vf_long <- melt(res[, .(xo, `boundary DUP` = dup_vf, `internal DEL` = del_vf)],
                id.vars = "xo", variable.name = "junction", value.name = "vf",
                na.rm = TRUE)

p1 <- ggplot(vf_long, aes(xo, vf, colour = junction)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_colour_manual(values = c(`boundary DUP` = ACCENT, `internal DEL` = DELC),
                      name = NULL) +
  scale_x_continuous(breaks = res$xo, labels = lab) +
  scale_y_log10() +
  labs(title = "A  Read support dates each junction",
       subtitle = "the late deletion never catches the boundary it sits inside",
       x = "Replication round", y = "Read support (VF, log)") +
  theme_val()

p2 <- ggplot(res, aes(xo, amplicon_cn)) +
  geom_line(colour = ACCENT, linewidth = 0.9) + geom_point(size = 2.4, colour = ACCENT) +
  geom_hline(yintercept = MIN_CN_RATIO * 2, linetype = "dashed", colour = WARN) +
  annotate("text", x = 1, y = MIN_CN_RATIO * 2, vjust = -0.6, hjust = 0,
           label = "detection floor (3 x ploidy)", colour = WARN, size = 3) +
  scale_x_continuous(breaks = res$xo, labels = lab) + scale_y_log10() +
  labs(title = "B  Copy number and the detection floor",
       x = "Replication round", y = "Peak copy number (log)") +
  theme_val()

p3 <- ggplot(res, aes(xo, n_junctions)) +
  geom_col(fill = "#c9ced4", width = 0.65) +
  scale_x_continuous(breaks = res$xo, labels = lab) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(title = "C  Replication adds no junctions; the deletion adds exactly one",
       x = "Replication round", y = "Junctions called") +
  theme_val()

p4 <- ggplot(res[!is.na(del_share)], aes(xo, del_share)) +
  geom_col(fill = DELC, width = 0.6) +
  scale_x_continuous(breaks = res$xo, labels = lab) +
  labs(title = "D  The deleted subclone expands",
       subtitle = sprintf("shorter circle, fitness %.2f; a neutral one would sit flat",
                          DEL_FITNESS),
       x = "Replication round", y = "Share of the ecDNA population (%)") +
  theme_val()

## The amplicon itself, before and after the deletion
pre  <- sim_to_plot_inputs(traj[[ROUNDS_PRE + 1L]], noise = CLEAN, seed = 7)
post <- sim_to_plot_inputs(traj[[length(traj)]],    noise = CLEAN, seed = 7)
win <- matrix(c(54.2e6, 56.0e6), nrow = 1)

pdf("validation/output/episome_replication.pdf", width = 11, height = 8, bg = "white")
print((p1 | p2) / (p3 | p4) + plot_annotation(
  title = "A clean episome: replication, then a late internal deletion",
  theme = theme(plot.background = element_rect(fill = "white", colour = NA))))
print(q(plot_sv_linear(sample = "EPI01", cnv_data = pre$cnv_data,
                       sv_data = pre$sv_data, wgd_data = pre$wgd_data,
                       chromosome = "chr7", chromosome_range = win)))
print(q(plot_sv_linear(sample = "EPI01", cnv_data = post$cnv_data,
                       sv_data = post$sv_data, wgd_data = post$wgd_data,
                       chromosome = "chr7", chromosome_range = win)))
invisible(dev.off())

message("Wrote validation/output/episome_replication.pdf")
