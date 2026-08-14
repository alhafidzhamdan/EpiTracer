## ---------------------------------------------------------------------------
## EpiTracer validation — difficulty sweep (robustness to noise)
##
## The clean benchmark (simulate_benchmark.R) separates the classes perfectly.
## Real copy-number and read-support estimates are noisy, so here we add
## increasing Gaussian noise to segment copy number and log-normal noise to the
## breakpoint variant fraction, and track how sensitivity / specificity / PPV
## degrade. A tool that fails gracefully (rather than cliff-edging) is the point.
##
## Run from the package root:  Rscript validation/simulate_sweep.R
## Writes validation/output/{sweep_metrics.csv, sweep.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2)
})
devtools::load_all(".", quiet = TRUE)

set.seed(1)
N_PER_CLASS <- 60                                  # episomal amplicons per noise level
## Copy-number noise SD in copies. PURPLE segment-level CN error is typically
## ~0.2-0.5 copies, so this grid spans the realistic regime and a little beyond.
NOISE       <- c(0, 0.1, 0.25, 0.5, 0.75, 1.0)     # VF CV = 0.1 * noise

jitter_cn <- function(x, sd) pmax(0, x + rnorm(length(x), 0, sd))

build_amplicon <- function(id, type, cn_sd, vf_cv) {
  amp_cn   <- round(runif(1, 20, 80))
  bnd_vf   <- round(runif(1, 600, 1500))
  flank_cn <- if (type == "gained_flanks") round(runif(1, 7, 15)) else 2
  s <- 55000000L; e <- 55500000L

  ## noisy copy numbers
  fcn <- jitter_cn(rep(flank_cn, 2), cn_sd); acn <- jitter_cn(amp_cn, cn_sd)
  ## noisy VF (log-normal multiplicative); the two breakends of one SV share a
  ## VF, as in real data, so noise is drawn per SV rather than per breakend.
  vfn <- function(v) round(v * exp(rnorm(1L, 0, vf_cv)))
  scar_vf <- vfn(bnd_vf * 0.2); dup_vf <- vfn(bnd_vf)

  ecdna <- data.table(seqnames = "chr7", start = s, end = e, ID = id, WGS_ID = id)
  bp <- data.table(
    seqnames = "chr7", start = c(s - 1e5L, s, e, e + 1e5L), end = c(s - 1e5L, s, e, e + 1e5L),
    WGS_ID = id, event = c("DEL1","DUP1","DUP1","DEL1"), svclass = c("DEL","DUP","DUP","DEL"),
    PURPLE_AF = c(0.4,0.9,0.9,0.4), PURPLE_JCN = c(1,40,40,1),
    VF = c(scar_vf, dup_vf, dup_vf, scar_vf),
    PURPLE_CN = c(fcn[1], acn, acn, fcn[2]), insLen = 0L, HOMLEN = 0L)
  if (type == "internal_higher_vf")
    bp <- rbind(bp, data.table(seqnames="chr7", start=55250000L, end=55250000L, WGS_ID=id,
      event="DUP2", svclass="DUP", PURPLE_AF=0.95, PURPLE_JCN=50, VF=vfn(bnd_vf + 400),
      PURPLE_CN=acn, insLen=0L, HOMLEN=0L))
  if (type == "no_boundary_dup") bp[event=="DUP1", `:=`(event="INV1", svclass="INV")]

  cnv <- data.table(seqnames="chr7", start=c(40000000L, s, e+1L), end=c(s-1L, e, 70000000L),
    sample=id, copyNumber=c(fcn[1], acn, fcn[2]), ploidy=2,
    majorAlleleCopyNumber=pmax(1, c(fcn[1], acn, fcn[2]) - 1), minorAlleleCopyNumber=1)

  list(ecdna=ecdna, bp=bp, cnv=cnv,
       truth=data.table(ID=id, truth_episomal = type == "episomal"))
}

run_level <- function(noise) {
  cn_sd <- noise; vf_cv <- 0.1 * noise
  types <- c(rep("episomal", N_PER_CLASS),
             rep(c("gained_flanks","internal_higher_vf","no_boundary_dup"),
                 each = ceiling(N_PER_CLASS/3)))
  ids  <- sprintf("N%03d_%04d", round(noise * 100), seq_along(types))
  amps <- Map(function(i,t) build_amplicon(i, t, cn_sd, vf_cv), ids, types)

  ecdna_gr <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"ecdna")), keep.extra.columns=TRUE)
  bp_gr    <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"bp")),    keep.extra.columns=TRUE)
  cnv_gr   <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"cnv")),   keep.extra.columns=TRUE)
  genes    <- GRanges("chr7", IRanges(55019017, 55211628), gene="EGFR")
  truth    <- rbindlist(lapply(amps,`[[`,"truth"))

  res  <- suppressWarnings(suppressMessages(call_episomal_ecdna(ecdna_gr, bp_gr, cnv_gr, genes, ext=1e7)))
  pred <- res[, .(pred = any(episomal == "TRUE")), by = ID]
  ev   <- merge(truth, pred, by="ID", all.x=TRUE); ev[is.na(pred), pred := FALSE]

  TP<-ev[truth_episomal & pred,.N]; FN<-ev[truth_episomal & !pred,.N]
  TN<-ev[!truth_episomal & !pred,.N]; FP<-ev[!truth_episomal & pred,.N]
  data.table(noise=noise,
             Sensitivity = TP/(TP+FN), Specificity = TN/(TN+FP),
             PPV = TP/(TP+FP), F1 = 2*TP/(2*TP+FP+FN))
}

message("Sweeping noise levels: ", paste(NOISE, collapse=", "), " ...")
sweep <- rbindlist(lapply(NOISE, run_level))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(sweep, "validation/output/sweep_metrics.csv")
cat("\n=== Metrics vs noise ===\n"); print(sweep)

long <- melt(sweep, id.vars="noise", variable.name="metric", value.name="value")
p <- ggplot(long, aes(noise, value, colour=metric)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values=c(Sensitivity="#0f7d78", Specificity="#c0851e",
                               PPV="#7048a6", F1="#0d1117")) +
  scale_y_continuous(limits=c(min(0.8, min(long$value)), 1.001)) +
  labs(title="EpiTracer robustness to copy-number / VF noise",
       subtitle=sprintf("%d amplicons per class per level", N_PER_CLASS),
       x="copy-number noise (SD, copies)", y="metric", colour=NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.background = element_rect(fill="white", colour=NA))
ggsave("validation/output/sweep.pdf", p, width=8, height=4.6, bg="white")
cat("\nWrote validation/output/{sweep_metrics.csv, sweep.pdf}\n")
