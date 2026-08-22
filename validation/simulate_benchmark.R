## ---------------------------------------------------------------------------
## EpiTracer validation — simulation benchmark (sensitivity / specificity / PPV)
##
## Generates focal amplicons with KNOWN origin and scores call_simple_excision():
##   * episomal           — boundary DUP (highest VF), diploid flanks, excision scar
##   * gained_flanks      — same boundary DUP + scar, but amplified flanks
##                          (a complex / chromothriptic context, NOT episomal)
##   * internal_higher_vf — an internal DUP outranks the boundary DUP (NOT episomal)
##   * no_boundary_dup    — no duplication at the amplicon boundary (NOT episomal)
##
## Each amplicon is an independent sample, so the caller processes them in
## isolation. Truth labels are known by construction, giving a confusion matrix.
##
## Run from the package root:  Rscript validation/simulate_benchmark.R
## Writes validation/output/{metrics.csv, confusion.csv, benchmark.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
})
devtools::load_all(".", quiet = TRUE)

set.seed(1)                 # reproducible
N_PER_CLASS <- 150          # amplicons per truth class

## Build one amplicon of a given type on chr7 around the EGFR locus. Each gets a
## unique sample id so the caller treats them independently.
build_amplicon <- function(id, type) {
  amp_cn   <- round(runif(1, 20, 80))          # amplicon copy number
  bnd_vf   <- round(runif(1, 600, 1500))       # boundary DUP read support (VF)
  flank_cn <- if (type == "gained_flanks") round(runif(1, 7, 15)) else 2L
  s <- 55000000L; e <- 55500000L               # amplicon boundaries

  ecdna <- data.table(seqnames = "chr7", start = s, end = e,
                      ID = id, WGS_ID = id)

  ## breakpoints: boundary DUP pair + a low-VF deletion pair (the excision scar)
  bp <- data.table(
    seqnames = "chr7",
    start    = c(s - 100000L, s, e, e + 100000L),
    end      = c(s - 100000L, s, e, e + 100000L),
    WGS_ID   = id,
    event    = c("DEL1", "DUP1", "DUP1", "DEL1"),
    svclass  = c("DEL",  "DUP",  "DUP",  "DEL"),
    PURPLE_AF  = c(0.4, 0.9, 0.9, 0.4),
    PURPLE_JCN = c(1, 40, 40, 1),
    VF         = c(round(bnd_vf * 0.2), bnd_vf, bnd_vf, round(bnd_vf * 0.2)),
    PURPLE_CN  = c(flank_cn, amp_cn, amp_cn, flank_cn),
    insLen = 0L, HOMLEN = 0L
  )

  if (type == "internal_higher_vf") {
    ## an internal DUP with VF above the boundary -> boundary is not the highest
    bp <- rbind(bp, data.table(
      seqnames = "chr7", start = 55250000L, end = 55250000L, WGS_ID = id,
      event = "DUP2", svclass = "DUP",
      PURPLE_AF = 0.95, PURPLE_JCN = 50, VF = bnd_vf + 400L,
      PURPLE_CN = amp_cn, insLen = 0L, HOMLEN = 0L))
  }
  if (type == "no_boundary_dup") {
    ## replace the boundary duplications with (non-DUP) inversions
    bp[event == "DUP1", `:=`(event = "INV1", svclass = "INV")]
  }

  cnv <- data.table(
    seqnames = "chr7",
    start = c(40000000L, s, e + 1L),
    end   = c(s - 1L, e, 70000000L),
    sample = id,
    copyNumber            = c(flank_cn, amp_cn, flank_cn),
    ploidy                = 2,
    majorAlleleCopyNumber = c(flank_cn - 1, amp_cn - 1, flank_cn - 1),
    minorAlleleCopyNumber = c(1, 1, 1)
  )

  truth <- type == "episomal"
  list(ecdna = ecdna, bp = bp, cnv = cnv,
       truth = data.table(ID = id, type = type, truth_episomal = truth))
}

## Assemble the cohort: N episomal + N split across the three non-episomal types
types <- c(rep("episomal", N_PER_CLASS),
           rep(c("gained_flanks", "internal_higher_vf", "no_boundary_dup"),
               each = ceiling(N_PER_CLASS / 3)))
ids   <- sprintf("SIM%04d", seq_along(types))
amps  <- Map(build_amplicon, ids, types)

ecdna_gr <- makeGRangesFromDataFrame(rbindlist(lapply(amps, `[[`, "ecdna")),
                                     keep.extra.columns = TRUE)
breakpoints_gr <- makeGRangesFromDataFrame(rbindlist(lapply(amps, `[[`, "bp")),
                                           keep.extra.columns = TRUE)
cnv_gr <- makeGRangesFromDataFrame(rbindlist(lapply(amps, `[[`, "cnv")),
                                   keep.extra.columns = TRUE)
cancer_genes_gr <- GRanges("chr7", IRanges(55019017, 55211628), gene = "EGFR")
truth <- rbindlist(lapply(amps, `[[`, "truth"))

## Run EpiTracer on the whole simulated cohort
message("Calling ", length(ecdna_gr), " simulated amplicons ...")
res <- suppressWarnings(suppressMessages(
  call_simple_excision(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
                      ext = 1e7, mc.cores = 1)))

## Per-amplicon prediction: episomal if any breakpoint of the amplicon is TRUE
pred <- res[, .(pred_episomal = any(episomal == "TRUE")), by = ID]
ev <- merge(truth, pred, by = "ID", all.x = TRUE)
ev[is.na(pred_episomal), pred_episomal := FALSE]     # undetected -> not episomal

## Confusion matrix + metrics
TP <- ev[truth_episomal & pred_episomal, .N]
FN <- ev[truth_episomal & !pred_episomal, .N]
TN <- ev[!truth_episomal & !pred_episomal, .N]
FP <- ev[!truth_episomal & pred_episomal, .N]

metrics <- data.table(
  metric = c("Sensitivity (recall)", "Specificity", "Precision (PPV)",
             "NPV", "Accuracy", "F1", "N amplicons"),
  value  = c(TP / (TP + FN), TN / (TN + FP), TP / (TP + FP),
             TN / (TN + FN), (TP + TN) / (TP + TN + FP + FN),
             2 * TP / (2 * TP + FP + FN), nrow(ev))
)
confusion <- data.table(TP = TP, FP = FP, FN = FN, TN = TN)

## Per-type recovery (how each non-episomal mechanism is handled)
by_type <- ev[, .(n = .N, called_episomal = sum(pred_episomal)), by = type]

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(metrics,   "validation/output/metrics.csv")
fwrite(confusion, "validation/output/confusion.csv")
fwrite(by_type,   "validation/output/by_type.csv")

cat("\n=== Confusion matrix ===\n"); print(confusion)
cat("\n=== Metrics ===\n"); print(metrics[, .(metric, value = round(value, 3))])
cat("\n=== Per truth class ===\n"); print(by_type)

## Figure: confusion matrix tile + metrics bar (explicit white background)
cm <- data.table(
  truth = c("episomal", "episomal", "not episomal", "not episomal"),
  pred  = c("episomal", "not episomal", "episomal", "not episomal"),
  n     = c(TP, FN, FP, TN)
)
p_cm <- ggplot(cm, aes(pred, truth, fill = n)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = n), size = 6) +
  scale_fill_gradient(low = "#eef1f4", high = "#0f7d78", guide = "none") +
  labs(title = "Confusion matrix", x = "EpiTracer call", y = "Simulated truth") +
  theme_minimal(base_size = 12)

mb <- metrics[metric != "N amplicons"]
mb[, metric := factor(metric, levels = rev(metric))]
p_mb <- ggplot(mb, aes(value, metric)) +
  geom_col(fill = "#0f7d78", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", value)), hjust = -0.15, size = 4) +
  scale_x_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
  labs(title = "Classification performance", x = NULL, y = NULL) +
  theme_minimal(base_size = 12)

fig <- patchwork::wrap_plots(p_cm, p_mb, widths = c(1, 1.2)) +
  patchwork::plot_annotation(
    title = sprintf("EpiTracer simulation benchmark (n = %d amplicons)", nrow(ev)),
    theme = ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA)))
ggsave("validation/output/benchmark.pdf", fig, width = 11, height = 4.5,
       bg = "white")
cat("\nWrote validation/output/{metrics.csv, confusion.csv, by_type.csv, benchmark.pdf}\n")
