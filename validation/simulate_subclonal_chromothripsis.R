## ---------------------------------------------------------------------------
## EpiTracer validation -- when does the founder boundary DUP survive internal
## chromothripsis and stand high above lower-VF internal SVs?
##
## This is the DO11441T1 EGFR pattern: a boundary duplication at VF 1462 / JCN 91
## with a scatter of internal SVs at VF 49-70 / JCN 2.2-2.7. The question is what
## has to be true for an amplicon to look like that.
##
## The answer follows from one fact: every junction carried by the SAME molecule
## has the same junction copy number, and therefore the same read support. A
## boundary at 20x the read support of the internal SVs cannot come from a single
## shattered circle -- on one molecule they would all be level. It requires the
## internal junctions to be SUBCLONAL: only a fraction of the ecDNA population
## was micronucleated and shattered, while the intact founder circles carried on.
##
## That makes the ratio quantitative:
##
##     boundary VF / internal VF  ~=  intact circles / shattered circles
##
## so the observed ratio ESTIMATES the fraction of the ecDNA population that has
## been through a micronucleus. This script sweeps that fraction and the
## shattering intensity, and checks the estimate against the two real amplicons
## in example_data/.
##
## Run from the package root:
##   Rscript validation/simulate_subclonal_chromothripsis.R
## Writes validation/output/{subclonal_chromothripsis.csv,
## subclonal_chromothripsis.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2); library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
q <- function(expr) suppressWarnings(suppressMessages(expr))

TOTAL   <- 95         # intact founder circles per cell (DO11441T1-like)
FLANK   <- 1.2e6      # a 2.6 Mb circle
HOST_PL <- 3          # DO11441T1 carries a whole-chr7 gain
SHATTERED <- c(1, 2, 3, 5, 8, 12, 20, 30, 45)   # circles that went through the micronucleus
BREAKS    <- c(5L, 10L, 20L, 40L)
N_REP     <- 8
CLEAN <- sim_noise(fp_rate = 0, dropout = 0)

## Real anchors, read straight from the cohort data.
REAL <- local({
  sv <- as.data.table(readRDS("example_data/all_348_SV_bedpe.rds"))
  f <- function(s, chr, lo, hi) {
    e <- sv[sample == s & chrom1 == chr & start1 > lo & start1 < hi]
    if (!nrow(e)) return(NULL)
    e <- e[order(-VF)]
    data.table(sample = s, founder_vf = e$VF[1], founder_jcn = e$JCN[1],
               internal_vf = e$VF[2], internal_jcn = e$JCN[2],
               vf_ratio = e$VF[1] / e$VF[2], jcn_ratio = e$JCN[1] / e$JCN[2])
  }
  rbindlist(list(f("DO11441T1", "7", 52e6, 56e6),
                 f("DO11501T1", "12", 57e6, 60e6)), fill = TRUE)
})
cat("=== Real amplicons ===\n"); print(REAL)
cat(sprintf("\nImplied shattered fraction (1 / jcn_ratio): %s\n",
            paste(sprintf("%s %.1f%%", REAL$sample, 100 / REAL$jcn_ratio),
                  collapse = "; ")))

theme_val <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background   = element_rect(fill = "white", colour = NA),
          panel.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 11))
}

## ---------------------------------------------------------------------------
## Sweep: how many circles were shattered x how hard
## ---------------------------------------------------------------------------
message("Sweeping subclonal shattering ...")

grid <- CJ(shattered = SHATTERED, n_breaks = BREAKS, rep = seq_len(N_REP))

res <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  k <- grid$shattered[i]; nb <- grid$n_breaks[i]; r <- grid$rep[i]
  ec <- sim_episome(seed_locus("EGFR", flank = FLANK), sample = "SUB",
                    copies = TOTAL, host_ploidy = HOST_PL, seed = 1000 + r)
  founder <- ec$junctions[origin == "amplicon"][1]

  ## a minority of circles goes through the micronucleus and religates;
  ## the intact founder population carries on unchanged
  sub <- try(sim_shatter_to_episomes(ec, n_breaks = nb, n_circles = 1L,
                                     keep_parent = TRUE, copies = k,
                                     seed = 2000 + r * 100 + nb), silent = TRUE)
  if (inherits(sub, "try-error")) return(NULL)

  inp <- sim_to_epitracer(sub, noise = CLEAN, seed = 3000 + i)
  bp  <- gr2dt(inp$breakpoints_gr)
  amp <- bp[PURPLE_CN > 3 * HOST_PL]
  if (!nrow(amp)) return(NULL)

  ## the founder junction, identified by its exact breakend coordinates
  is_f <- amp$start %in% c(founder$start1, founder$start2) & amp$svclass == "DUP"
  if (!any(is_f)) return(NULL)
  fvf <- max(amp[is_f]$VF); fjcn <- max(amp[is_f]$PURPLE_JCN)
  int <- amp[!is_f]
  if (!nrow(int)) return(NULL)

  data.table(shattered = k, n_breaks = nb, rep = r,
             shattered_frac = k / (TOTAL + k),
             founder_vf = fvf, founder_jcn = fjcn,
             internal_vf_max = max(int$VF), internal_jcn_max = max(int$PURPLE_JCN),
             n_internal_sv = length(unique(int$event)),
             vf_ratio = fvf / max(int$VF),
             founder_is_top = fvf >= max(int$VF))
}), fill = TRUE)

fwrite(res, "validation/output/subclonal_chromothripsis.csv")

summ <- res[, .(vf_ratio = median(vf_ratio),
                internal_vf = median(internal_vf_max),
                n_internal_sv = median(n_internal_sv),
                founder_top = mean(founder_is_top)),
            by = .(shattered, shattered_frac, n_breaks)][order(shattered, n_breaks)]

cat("\n=== Founder : internal read-support ratio ===\n")
print(dcast(summ, shattered + shattered_frac ~ n_breaks, value.var = "vf_ratio"))
cat("\nfounder junction has the highest VF in",
    sprintf("%.0f%%", 100 * mean(res$founder_is_top)), "of runs\n")

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------
ACCENT <- "#0f7d78"; WARN <- "#c2571a"

p1 <- ggplot(summ, aes(100 * shattered_frac, vf_ratio, colour = factor(n_breaks))) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  geom_hline(yintercept = REAL$jcn_ratio, linetype = "dashed", colour = WARN) +
  annotate("text", x = 1.2, y = REAL$jcn_ratio, hjust = 0, vjust = -0.5,
           label = REAL$sample, colour = WARN, size = 3) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_viridis_d(name = "breakpoints", option = "D", end = 0.85) +
  labs(title = "A  The VF ratio measures the shattered fraction",
       subtitle = "ratio ~= intact circles / shattered circles, independent of how hard it shattered",
       x = "Share of ecDNA population shattered (%, log)",
       y = "Founder VF / top internal VF (log)") +
  theme_val()

p2 <- ggplot(summ, aes(n_breaks, n_internal_sv, colour = factor(shattered))) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_colour_viridis_d(name = "circles\nshattered", option = "C", end = 0.85) +
  labs(title = "B  Shattering intensity sets the junction COUNT",
       subtitle = "how many internal SVs, not how far below the founder they sit",
       x = "Shattering breakpoints", y = "Internal SVs called") +
  theme_val()

p3 <- ggplot(res, aes(factor(shattered), internal_vf_max)) +
  geom_boxplot(fill = "#eef1f4", colour = "#42474d", outlier.size = 0.7) +
  geom_hline(yintercept = median(res$founder_vf), linetype = "dashed", colour = ACCENT) +
  annotate("text", x = 0.7, y = median(res$founder_vf), hjust = 0, vjust = -0.5,
           label = "founder VF", colour = ACCENT, size = 3) +
  scale_y_log10() +
  labs(title = "C  Internal SV support rises with the shattered fraction",
       subtitle = "the founder stays put -- it rides every circle",
       x = "Circles shattered (of 95)", y = "Top internal SV VF (log)") +
  theme_val()

demo <- local({
  ec <- sim_episome(seed_locus("EGFR", flank = FLANK), sample = "SUBDEMO",
                    copies = TOTAL, host_ploidy = HOST_PL, seed = 1)
  s <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 1L,
                               keep_parent = TRUE, copies = 5, seed = 7)
  sim_to_plot_inputs(s, noise = CLEAN, seed = 7)
})

pdf("validation/output/subclonal_chromothripsis.pdf", width = 11, height = 8, bg = "white")
print((p1 | p2) / p3 + plot_annotation(
  title = "A high-VF boundary DUP above low-VF internal SVs requires subclonal shattering",
  theme = theme(plot.background = element_rect(fill = "white", colour = NA))))
print(q(plot_sv_linear(sample = "SUBDEMO", cnv_data = demo$cnv_data,
                       sv_data = demo$sv_data, wgd_data = demo$wgd_data,
                       chromosome = "chr7",
                       chromosome_range = matrix(c(53.4e6, 56.9e6), nrow = 1))))
invisible(dev.off())

message("Wrote validation/output/subclonal_chromothripsis.pdf")
