## ---------------------------------------------------------------------------
## EpiTracer validation -- an episome shattered into new episomes
##
## A clean episomal ecDNA (one boundary duplication, flat high copy number) is
## encapsulated in a micronucleus, shattered, and the fragments religate into
## SEVERAL new circles rather than one. Each daughter is a genuine episome with
## its own circularisation junction; the ones that keep the oncogene are selected.
##
## The question: what does EpiTracer see afterwards?
##
## The answer matters, because the founder boundary duplication only survives if
## some daughter happens to rejoin the exact two ends it joined. Once the circle
## has been cut in twenty places that essentially never happens -- so the amplicon
## is still, genuinely, a set of episomes, but the single spanning
## duplication-orientation junction that call_simple_excision() keys on is gone.
##
## Run from the package root:  Rscript validation/simulate_episome_fission.R
## Writes validation/output/{fission_metrics.csv, fission_by_nbreaks.csv,
## episome_fission.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2); library(patchwork)
})
devtools::load_all(".", quiet = TRUE)

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
q <- function(expr) suppressWarnings(suppressMessages(expr))

N_REP    <- 30
GENE     <- "EGFR"
COPIES   <- 95        # DO11441T1-like
FLANK    <- 1.2e6     # a 2.6 Mb circle, like the real amplicon
HOST_PL  <- 3         # DO11441T1 carries a whole-chr7 gain
CLEAN    <- sim_noise(fp_rate = 0, dropout = 0)

theme_val <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background   = element_rect(fill = "white", colour = NA),
          panel.background  = element_rect(fill = "white", colour = NA),
          legend.background = element_rect(fill = "white", colour = NA),
          plot.title = element_text(face = "bold", size = 11))
}

## Did any daughter rejoin the two ends the founder duplication joined?
founder_survives <- function(parent, child) {
  f <- parent$junctions[origin == "amplicon"]
  if (!nrow(f)) return(FALSE)
  nrow(merge(child$junctions, f[, .(chrom1, start1, chrom2, start2)],
             by = c("chrom1", "start1", "chrom2", "start2"))) > 0L
}

score <- function(x, seed) {
  inp <- sim_to_epitracer(x, noise = CLEAN, seed = seed)
  if (!length(inp$ecdna_gr)) return(NULL)
  a <- list(inp$ecdna_gr, inp$breakpoints_gr, inp$cnv_gr, inp$cancer_genes_gr)
  se <- q(do.call(call_simple_excision, a))
  ct <- q(do.call(call_chromothripsis, a))
  bp <- gr2dt(inp$breakpoints_gr)[!duplicated(event)]
  list(episomal   = nrow(se) > 0 && any(se$episomal == "TRUE"),
       boundary   = nrow(se) > 0 && any(se$duplication_at_boundary == "TRUE"),
       chromothr  = nrow(ct) > 0 && any(ct$chromothripsis == "TRUE"),
       n_sv       = nrow(bp),
       max_cn     = max(inp$cnv_gr$copyNumber),
       n_seeds    = length(inp$ecdna_gr))
}

## ---------------------------------------------------------------------------
## 1. Parent vs daughters, over replicates
## ---------------------------------------------------------------------------
message("Simulating ", N_REP, " fission events ...")

res <- rbindlist(lapply(seq_len(N_REP), function(i) {
  ec <- sim_episome(seed_locus(GENE, flank = FLANK), sample = sprintf("FIS%02d", i),
                    copies = COPIES, host_ploidy = HOST_PL, seed = 1000 + i)
  fis <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 3L, seed = 2000 + i)

  p <- score(ec,  3000 + i); d <- score(fis, 4000 + i)
  if (is.null(p) || is.null(d)) return(NULL)
  data.table(
    replicate = i,
    n_daughters = length(EpiTracer:::.sim_population(fis)),
    founder_kept = founder_survives(ec, fis),
    parent_episomal = p$episomal, parent_boundary = p$boundary,
    daughter_episomal = d$episomal, daughter_boundary = d$boundary,
    daughter_chromothripsis = d$chromothr,
    parent_sv = p$n_sv, daughter_sv = d$n_sv,
    parent_cn = p$max_cn, daughter_cn = d$max_cn)
}), fill = TRUE)

fwrite(res, "validation/output/fission_metrics.csv")

summ <- res[, .(n = .N,
                mean_daughters = mean(n_daughters),
                founder_kept = mean(founder_kept),
                parent_episomal = mean(parent_episomal),
                daughter_episomal = mean(daughter_episomal),
                daughter_boundary_dup = mean(daughter_boundary),
                daughter_chromothripsis = mean(daughter_chromothripsis),
                median_sv_parent = median(parent_sv),
                median_sv_daughter = median(daughter_sv))]
cat("\n=== Parent episome vs its fission products ===\n"); print(t(summ))

## ---------------------------------------------------------------------------
## 2. How much shattering does it take to lose the founder junction?
## ---------------------------------------------------------------------------
message("Sweeping shattering intensity ...")

BREAKS <- c(1L, 2L, 3L, 5L, 8L, 12L, 20L, 30L)
sweep <- rbindlist(lapply(BREAKS, function(nb) {
  out <- rbindlist(lapply(seq_len(N_REP), function(i) {
    ec <- sim_episome(seed_locus(GENE, flank = FLANK), sample = sprintf("SW%02d", i),
                      copies = COPIES, host_ploidy = HOST_PL, seed = 5000 + i)
    fis <- try(sim_shatter_to_episomes(ec, n_breaks = nb, n_circles = 2L,
                                       seed = 6000 + i * 100 + nb), silent = TRUE)
    if (inherits(fis, "try-error")) return(NULL)
    d <- score(fis, 7000 + i * 100 + nb)
    if (is.null(d)) return(NULL)
    data.table(founder_kept = founder_survives(ec, fis),
               episomal = d$episomal, boundary = d$boundary)
  }), fill = TRUE)
  if (!nrow(out)) return(NULL)
  data.table(n_breaks = nb, n = nrow(out),
             founder_kept = mean(out$founder_kept),
             still_episomal = mean(out$episomal),
             boundary_dup = mean(out$boundary))
}), fill = TRUE)

fwrite(sweep, "validation/output/fission_by_nbreaks.csv")
cat("\n=== Founder junction survival vs shattering intensity ===\n"); print(sweep)

## ---------------------------------------------------------------------------
## 3. Figures
## ---------------------------------------------------------------------------
ACCENT <- "#0f7d78"; WARN <- "#c2571a"

calls <- data.table(
  stage = factor(rep(c("parent episome", "fission products"), each = 3),
                 levels = c("parent episome", "fission products")),
  caller = factor(rep(c("episomal", "boundary DUP", "chromothripsis"), 2),
                  levels = c("episomal", "boundary DUP", "chromothripsis")),
  rate = c(summ$parent_episomal, mean(res$parent_boundary), 0,
           summ$daughter_episomal, summ$daughter_boundary_dup,
           summ$daughter_chromothripsis))

p1 <- ggplot(calls, aes(caller, rate, fill = stage)) +
  geom_col(position = position_dodge(0.7), width = 0.62) +
  scale_fill_manual(values = c(`parent episome` = ACCENT,
                               `fission products` = WARN), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "A  Fission destroys the episomal call",
       subtitle = "the products are still circles -- the spanning duplication is not",
       x = NULL, y = "Amplicons called") +
  theme_val()

p2 <- ggplot(sweep, aes(n_breaks)) +
  geom_line(aes(y = founder_kept, colour = "founder DUP retained"), linewidth = 0.9) +
  geom_point(aes(y = founder_kept, colour = "founder DUP retained"), size = 2.3) +
  geom_line(aes(y = still_episomal, colour = "still called episomal"), linewidth = 0.9) +
  geom_point(aes(y = still_episomal, colour = "still called episomal"), size = 2.3) +
  scale_colour_manual(values = c(`founder DUP retained` = ACCENT,
                                 `still called episomal` = WARN), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_log10(breaks = BREAKS) +
  labs(title = "B  How much shattering it takes",
       x = "Shattering breakpoints (log)", y = "Rate") +
  theme_val()

p3 <- ggplot(melt(res[, .(replicate, parent = parent_sv, daughters = daughter_sv)],
                  id.vars = "replicate", variable.name = "stage", value.name = "n_sv"),
             aes(stage, n_sv, fill = stage)) +
  geom_boxplot(width = 0.55, outlier.size = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c(parent = ACCENT, daughters = WARN)) +
  labs(title = "C  Junction burden", x = NULL, y = "Junctions called") +
  theme_val()

p4 <- ggplot(res, aes(factor(n_daughters))) +
  geom_bar(fill = WARN, width = 0.6) +
  labs(title = "D  Daughters surviving selection",
       subtitle = "of 3 attempted; those without the oncogene are lost",
       x = "Surviving circle species", y = "Replicates") +
  theme_val()

## One worked example, parent and products
ec  <- sim_episome(seed_locus(GENE, flank = FLANK), sample = "FISDEMO",
                   copies = COPIES, host_ploidy = HOST_PL, seed = 1)
fis <- sim_shatter_to_episomes(ec, n_breaks = 20L, n_circles = 3L, seed = 7)
pin_p <- sim_to_plot_inputs(ec,  noise = CLEAN, seed = 1)
pin_d <- sim_to_plot_inputs(fis, noise = CLEAN, seed = 7)
win <- matrix(c(53.4e6, 56.9e6), nrow = 1)

pdf("validation/output/episome_fission.pdf", width = 11, height = 8, bg = "white")
print((p1 | p2) / (p3 | p4) + plot_annotation(
  title = "An episome shattered into new episomes",
  theme = theme(plot.background = element_rect(fill = "white", colour = NA))))
print(q(plot_sv_linear(sample = "FISDEMO", cnv_data = pin_p$cnv_data,
                       sv_data = pin_p$sv_data, wgd_data = pin_p$wgd_data,
                       chromosome = "chr7", chromosome_range = win)))
print(q(plot_sv_linear(sample = "FISDEMO", cnv_data = pin_d$cnv_data,
                       sv_data = pin_d$sv_data, wgd_data = pin_d$wgd_data,
                       chromosome = "chr7", chromosome_range = win)))
invisible(dev.off())

message("Wrote validation/output/episome_fission.pdf")
