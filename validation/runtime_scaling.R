## ---------------------------------------------------------------------------
## EpiTracer validation — runtime / scaling
##
## Wall-clock of call_simple_excision() as the number of amplicons grows, to show
## the tool runs reliably at cohort scale (for the manuscript's availability /
## implementation note). Uses the same synthetic amplicon builder as the
## benchmark; each amplicon is an independent sample.
##
## Run from the package root:  Rscript validation/runtime_scaling.R
## Writes validation/output/{scaling.csv, scaling.pdf}.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges); library(data.table); library(ggplot2)
})
devtools::load_all(".", quiet = TRUE)

set.seed(1)
SIZES <- c(10, 25, 50, 100, 200, 400)   # number of amplicons
CORES <- 1                               # single core (portable, worst case)

one_amplicon <- function(id) {
  amp_cn <- round(runif(1, 20, 80)); bnd_vf <- round(runif(1, 600, 1500))
  s <- 55000000L; e <- 55500000L
  list(
    ecdna = data.table(seqnames="chr7", start=s, end=e, ID=id, WGS_ID=id),
    bp = data.table(seqnames="chr7", start=c(s-1e5L,s,e,e+1e5L), end=c(s-1e5L,s,e,e+1e5L),
      WGS_ID=id, event=c("DEL1","DUP1","DUP1","DEL1"), svclass=c("DEL","DUP","DUP","DEL"),
      PURPLE_AF=c(0.4,0.9,0.9,0.4), PURPLE_JCN=c(1,40,40,1),
      VF=c(round(bnd_vf*0.2),bnd_vf,bnd_vf,round(bnd_vf*0.2)),
      PURPLE_CN=c(2,amp_cn,amp_cn,2), insLen=0L, HOMLEN=0L),
    cnv = data.table(seqnames="chr7", start=c(40000000L,s,e+1L), end=c(s-1L,e,70000000L),
      sample=id, copyNumber=c(2,amp_cn,2), ploidy=2,
      majorAlleleCopyNumber=c(1,amp_cn-1,1), minorAlleleCopyNumber=1))
}

genes <- GRanges("chr7", IRanges(55019017, 55211628), gene="EGFR")

bench_size <- function(n) {
  amps <- lapply(sprintf("A%04d", seq_len(n)), one_amplicon)
  ecdna_gr <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"ecdna")), keep.extra.columns=TRUE)
  bp_gr    <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"bp")),    keep.extra.columns=TRUE)
  cnv_gr   <- makeGRangesFromDataFrame(rbindlist(lapply(amps,`[[`,"cnv")),   keep.extra.columns=TRUE)
  t <- system.time(suppressWarnings(suppressMessages(
    call_simple_excision(ecdna_gr, bp_gr, cnv_gr, genes, ext=1e7, mc.cores=CORES))))[["elapsed"]]
  data.table(amplicons = n, seconds = round(t, 3), per_amplicon_ms = round(1000 * t / n, 1))
}

message("Timing sizes: ", paste(SIZES, collapse=", "), " amplicons (", CORES, " core) ...")
scaling <- rbindlist(lapply(SIZES, bench_size))
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(scaling, "validation/output/scaling.csv")
cat("\n=== Runtime scaling ===\n"); print(scaling)

p <- ggplot(scaling, aes(amplicons, seconds)) +
  geom_line(colour="#0f7d78", linewidth=1) + geom_point(colour="#0f7d78", size=2.4) +
  labs(title="EpiTracer runtime scales linearly in amplicons",
       subtitle=sprintf("single core; ~%.0f ms per amplicon", mean(scaling$per_amplicon_ms)),
       x="amplicons", y="wall-clock (s)") +
  theme_minimal(base_size=12) +
  theme(plot.background = element_rect(fill="white", colour=NA))
ggsave("validation/output/scaling.pdf", p, width=7, height=4.4, bg="white")
cat("\nWrote validation/output/{scaling.csv, scaling.pdf}\n")
