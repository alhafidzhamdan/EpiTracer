## ---------------------------------------------------------------------------
## EpiTracer validation -- simulated ecDNA trajectories
##
## Uses the package's own simulation suite (R/sim-*.R) to ask three questions
## that a hand-written scenario table cannot answer, because copy number and
## structural variants here are derived from ONE ground-truth fragment list and
## are therefore mutually consistent:
##
##   (1) At what point along an ecDNA's life does EpiTracer stop calling it a
##       simple excision and start calling it chromothriptic? A clean episome is
##       shattered in a micronucleus round after round; each round is rendered as
##       WGS-like inputs and pushed through every caller.
##   (2) Do two episomes fused in one micronucleus produce the chimeric,
##       translocation-bearing amplicon call_micronucleation() is built for?
##   (3) How do the calls hold up as tumour purity and depth fall?
##
## Run from the package root:  Rscript validation/simulate_trajectories.R
## Writes validation/output/{trajectory_metrics.csv, chimera_metrics.csv,
## purity_sweep.csv, simulate_trajectories.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2); library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

set.seed(2024)
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
q <- function(expr) suppressWarnings(suppressMessages(expr))

## Every figure on an explicit white canvas, so the PDF does not inherit a
## transparent background when dropped into a manuscript.
theme_val <- function(...) {
  theme_minimal(base_size = 11, ...) +
    theme(plot.background  = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 11))
}

N_REPLICATES <- 20      # independent trajectories per round
MAX_ROUNDS   <- 4
GENE         <- "EGFR"

## ---------------------------------------------------------------------------
## (1) Trajectory: episome -> repeated micronucleation + chromothripsis
## ---------------------------------------------------------------------------
message("Simulating ", N_REPLICATES, " trajectories over ", MAX_ROUNDS, " rounds ...")

traj <- rbindlist(lapply(seq_len(N_REPLICATES), function(rep) {
  ec <- sim_episome(seed_locus(GENE), sample = sprintf("TRJ%02d", rep),
                    copies = runif(1, 25, 70), seed = 1000 + rep)
  steps <- sim_evolve(ec, rounds = MAX_ROUNDS, n_breaks = sample(10:25, 1),
                      keep_all = TRUE, seed = 2000 + rep)
  rbindlist(lapply(seq_along(steps), function(i) {
    x <- steps[[i]]; x$sample <- sprintf("TRJ%02d_R%d", rep, i - 1L)
    inp <- sim_to_epitracer(x, seed = 3000 + rep * 10 + i)
    if (!length(inp$ecdna_gr)) return(NULL)
    a <- list(inp$ecdna_gr, inp$breakpoints_gr, inp$cnv_gr, inp$cancer_genes_gr)
    se <- q(do.call(call_simple_excision, a))
    ct <- q(do.call(call_chromothripsis, a))
    mn <- q(do.call(call_micronucleation, a))
    br <- q(do.call(call_brf, a))
    data.table(
      replicate = rep, round = i - 1L,
      truth_junctions = nrow(x$junctions),
      called_sv = length(unique(gr2dt(inp$breakpoints_gr)$event)),
      max_cn = summary(x)$max_cn,
      episomal        = any(se$episomal == "TRUE"),
      chromothripsis  = any(ct$chromothripsis == "TRUE"),
      micronucleation = any(mn$micronucleation == "TRUE"),
      brf             = any(br$brf == "TRUE"),
      oscillations = max(ct$cn_oscillations, na.rm = TRUE),
      joins_p = suppressWarnings(max(ct$sv_type_pval, na.rm = TRUE)))
  }), fill = TRUE)
}), fill = TRUE)

traj_summary <- traj[, .(
  n = .N,
  episomal = mean(episomal), chromothripsis = mean(chromothripsis),
  micronucleation = mean(micronucleation), brf = mean(brf),
  median_junctions = median(truth_junctions),
  median_called_sv = median(called_sv),
  median_oscillations = median(oscillations)), by = round][order(round)]

fwrite(traj_summary, "validation/output/trajectory_metrics.csv")
cat("\n=== Mechanism calls by round ===\n"); print(traj_summary)

## ---------------------------------------------------------------------------
## (2) Chimeric ecDNA: two episomes co-encapsulated in one micronucleus
## ---------------------------------------------------------------------------
message("Simulating chimeric episome fusions ...")

PAIRS <- list(c("EGFR", "CDK4"), c("MYC", "TERT"), c("MDM2", "SOX2"),
              c("PDGFRA", "MET"))

chim <- rbindlist(lapply(seq_len(N_REPLICATES), function(rep) {
  pg <- PAIRS[[(rep - 1) %% length(PAIRS) + 1]]
  id <- sprintf("CHM%02d", rep)
  a <- sim_episome(seed_locus(pg[1]), sample = id, copies = runif(1, 25, 60),
                   seed = 4000 + rep)
  b <- sim_episome(seed_locus(pg[2]), sample = id, copies = runif(1, 25, 60),
                   seed = 5000 + rep)
  ch <- sim_fuse_episomes(a, b, n_breaks = sample(8:20, 1), seed = 6000 + rep)
  inp <- sim_to_epitracer(ch, seed = 7000 + rep)
  if (!length(inp$ecdna_gr)) return(NULL)
  ar <- list(inp$ecdna_gr, inp$breakpoints_gr, inp$cnv_gr, inp$cancer_genes_gr)
  mn <- q(do.call(call_micronucleation, ar))
  ct <- q(do.call(call_chromothripsis, ar))
  se <- q(do.call(call_simple_excision, ar))
  data.table(replicate = rep, genes = paste(pg, collapse = "+"),
             n_chr = summary(ch)$n_chr, n_tra = summary(ch)$n_tra,
             micronucleation = any(mn$micronucleation == "TRUE"),
             chromothripsis  = any(ct$chromothripsis == "TRUE"),
             episomal        = any(se$episomal == "TRUE"))
}), fill = TRUE)

chim_summary <- chim[, .(n = .N, micronucleation = mean(micronucleation),
                         chromothripsis = mean(chromothripsis),
                         episomal = mean(episomal),
                         median_tra = median(n_tra)), by = genes]
fwrite(chim, "validation/output/chimera_metrics.csv")
cat("\n=== Chimeric ecDNA ===\n"); print(chim_summary)

## ---------------------------------------------------------------------------
## (3) Purity / depth sweep on a freshly born episome
## ---------------------------------------------------------------------------
message("Sweeping purity and depth ...")

GRID <- CJ(purity = c(0.2, 0.3, 0.5, 0.7, 0.9), depth = c(30, 60, 90))

sweep <- rbindlist(lapply(seq_len(nrow(GRID)), function(k) {
  p <- GRID$purity[k]; d <- GRID$depth[k]
  hits <- vapply(seq_len(N_REPLICATES), function(rep) {
    ec <- sim_episome(seed_locus(GENE), sample = sprintf("PS%02d", rep),
                      copies = runif(1, 20, 60), seed = 8000 + rep)
    inp <- sim_to_epitracer(ec, noise = sim_noise(purity = p, depth = d),
                            seed = 9000 + rep * 100 + k)
    if (!length(inp$ecdna_gr)) return(FALSE)
    se <- q(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr,
                                 inp$cnv_gr, inp$cancer_genes_gr))
    nrow(se) > 0 && any(se$episomal == "TRUE")
  }, logical(1))
  data.table(purity = p, depth = d, sensitivity = mean(hits), n = N_REPLICATES)
}))
fwrite(sweep, "validation/output/purity_sweep.csv")
cat("\n=== Purity / depth sweep (episomal sensitivity) ===\n"); print(sweep)

## ---------------------------------------------------------------------------
## Figures -- every diagnostic is shown as a plot, not only tabulated
## ---------------------------------------------------------------------------
PAL <- c(episomal = "#0f7d78", chromothripsis = "#c2571a",
         micronucleation = "#3f5d9e", brf = "#8a8f98")

long <- melt(traj_summary,
             id.vars = "round",
             measure.vars = c("episomal", "chromothripsis", "micronucleation", "brf"),
             variable.name = "caller", value.name = "rate")

p1 <- ggplot(long, aes(round, rate, colour = caller)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "A  Mechanism calls along the ecDNA trajectory",
       subtitle = sprintf("%d replicate %s episomes, shattered round after round",
                          N_REPLICATES, GENE),
       x = "Rounds of micronucleation + chromothripsis",
       y = "Amplicons called") +
  theme_val()

p2 <- ggplot(traj, aes(factor(round), truth_junctions)) +
  geom_boxplot(fill = "#eef1f4", colour = "#42474d", outlier.size = 0.7) +
  labs(title = "B  Junction burden accumulates", x = "Round",
       y = "True junctions on the circle") +
  theme_val()

p3 <- ggplot(traj, aes(factor(round), oscillations)) +
  geom_boxplot(fill = "#eef1f4", colour = "#42474d", outlier.size = 0.7) +
  geom_hline(yintercept = 3, linetype = "dashed", colour = "#c2571a") +
  labs(title = "C  Copy-number oscillation",
       subtitle = "dashed line: ShatterSeek threshold used by call_chromothripsis()",
       x = "Round", y = "Copy-number turning points") +
  theme_val()

p4 <- ggplot(chim_summary, aes(genes, micronucleation)) +
  geom_col(fill = PAL[["micronucleation"]], width = 0.65) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "D  Chimeric episomes are called micronucleated",
       subtitle = "two episomes fused in one micronucleus",
       x = NULL, y = "Called micronucleated") +
  theme_val() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p5 <- ggplot(sweep, aes(purity, sensitivity, colour = factor(depth))) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("30" = "#c2571a", "60" = "#0f7d78", "90" = "#3f5d9e"),
                      name = "Depth (x)") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "E  Episomal sensitivity vs tumour purity",
       x = "Tumour purity", y = "Episomes recovered") +
  theme_val()

## The simulated amplicon itself, at birth and after three rounds -- the visual
## check that a trajectory looks like a patient amplicon.
ec0 <- sim_episome(seed_locus(GENE), sample = "SHOWCASE", copies = 45, seed = 99)
ec3 <- sim_evolve(ec0, rounds = 3, n_breaks = 18, seed = 99)
pin0 <- sim_to_plot_inputs(ec0, seed = 99)
pin3 <- sim_to_plot_inputs(ec3, seed = 99)

pdf("validation/output/simulate_trajectories.pdf", width = 11, height = 8,
    bg = "white")
print((p1 | p2) / (p3 | p4) + plot_annotation(
  title = "EpiTracer: simulated ecDNA trajectories",
  theme = theme(plot.background = element_rect(fill = "white", colour = NA))))
print(p5)
print(q(plot_sv_linear(sample = "SHOWCASE", cnv_data = pin0$cnv_data,
                       sv_data = pin0$sv_data, wgd_data = pin0$wgd_data,
                       chromosome = "chr7",
                       chromosome_range = matrix(c(54e6, 56.5e6), nrow = 1))))
print(q(plot_sv_linear(sample = "SHOWCASE", cnv_data = pin3$cnv_data,
                       sv_data = pin3$sv_data, wgd_data = pin3$wgd_data,
                       chromosome = "chr7",
                       chromosome_range = matrix(c(54e6, 56.5e6), nrow = 1))))
invisible(dev.off())

message("Wrote validation/output/simulate_trajectories.pdf")
