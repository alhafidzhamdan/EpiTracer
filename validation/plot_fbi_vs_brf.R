## ---------------------------------------------------------------------------
## Fold-back-inversion (FBI) burden vs breakage-fusion (BRF) burden per amplicon
## across the cohort. Both are hallmarks of intrachromosomal breakage-fusion-
## bridge-like amplification, so they should correlate positively.
##   FBI burden = n_foldbacks         (call_bfb)
##   BRF burden = n_parallel_pairs    (call_brf, adjacent parallel breakpoints)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

d <- fread("validation/output/cohort_original_caller_per_chr.tsv")
d[is.na(n_foldbacks), n_foldbacks := 0L]
d[is.na(n_parallel_pairs), n_parallel_pairs := 0L]
d[, class := fifelse(episomal == TRUE, "episomal ecDNA", "non-episomal")]

## correlations (Spearman is rank-based and robust to the heavy skew; Pearson is
## reported on log1p, i.e. the slope of the log-log trend line drawn below)
sp <- suppressWarnings(cor.test(d$n_foldbacks, d$n_parallel_pairs, method = "spearman"))
pe <- suppressWarnings(cor.test(log1p(d$n_foldbacks), log1p(d$n_parallel_pairs)))
lab <- sprintf("Spearman rho = %.2f (p = %.2g)\nPearson r (log-log) = %.2f\nn = %d amplicons",
               unname(sp$estimate), sp$p.value, unname(pe$estimate), nrow(d))
message(gsub("\n", " | ", lab))

brks <- c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)
p <- ggplot(d, aes(log1p(n_foldbacks), log1p(n_parallel_pairs))) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", fill = "grey80", linewidth = 0.6) +
  geom_jitter(aes(colour = class), width = 0.06, height = 0.06, alpha = 0.75, size = 1.9) +
  annotate("text", x = log1p(0.2), y = log1p(3300), hjust = 0, vjust = 1, label = lab, size = 3.4) +
  scale_x_continuous(breaks = log1p(brks), labels = brks) +
  scale_y_continuous(breaks = log1p(brks), labels = brks) +
  scale_colour_manual(values = c("episomal ecDNA" = "#1E7A6F", "non-episomal" = "#C0504D"), name = NULL) +
  labs(title = "Fold-back-inversion vs BRF burden per amplicon (GBM cohort, per-chr caller)",
       x = "Fold-back inversions (n_foldbacks)", y = "BRF parallel-breakpoint pairs (n_parallel_pairs)") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey70"),
        panel.grid.minor = element_blank())

outpng <- "validation/output/fbi_vs_brf_scatter.png"
ggsave(outpng, p, width = 8, height = 6.5, dpi = 130, bg = "white")
message("Wrote ", outpng)
