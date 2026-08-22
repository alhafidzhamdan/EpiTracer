## ---------------------------------------------------------------------------
## EpiTracer validation — engineered ecDNA time course (ground-truth episome)
##
## The strongest positive control: an ecDNA formed by a DESIGNED excision-
## circularisation event, sampled over time. Pradella et al. engineered a
## Cre-loxP cassette in mouse neural progenitor cells (NPCs) that excises to
## circularise the Myc/Pvt1 locus; WGS was taken at 1-5 weeks after Cre
## induction, with a p53-Cre arm as control. EpiTracer should call the locus
## episomal exactly when the circle forms.
##
## DATA (mm10): AmpliconRepository "Pradella et al. engineered murine Myc ecDNA"
##   (project 6a5ea2166714848c8db0f348). Download + extract, then point this
##   script at results/samples. Samples: loxp15_cre_{1..5}w (engineered),
##   p53_cre_{1..5}w (control).
##
## USAGE (from the package root):
##   Rscript validation/engineered_ecdna_timecourse.R /path/to/results/samples
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(GenomicRanges); library(data.table); library(ggplot2)
})
source("validation/aa_to_epitracer.R")

args <- commandArgs(trailingOnly = TRUE)
samples_root <- if (length(args)) args[1] else stop("pass the results/samples folder")

## engineered locus (mm10): Myc + Pvt1 on chr15
onc <- GRanges("chr15", IRanges(c(61985391, 62038970), c(61990374, 62197684)),
               gene = c("Myc", "Pvt1"))
LOCUS <- list(chr = "chr15", start = 60525967, end = 62480542)

run_one <- function(sid) {
  sd <- file.path(samples_root, sid)
  gp  <- list.files(sd, "_amplicon[0-9]+_graph\\.txt$", full.names = TRUE, recursive = TRUE)
  cnvb <- list.files(sd, "CNV_CALLS\\.bed$", full.names = TRUE, recursive = TRUE)
  cnvb <- if (length(cnvb)) cnvb[1] else NULL
  ## locus max CN (from genome-wide CN)
  locus_cn <- NA_real_
  if (!is.null(cnvb)) {
    bed <- fread(cnvb, header = FALSE)
    hit <- bed[V1 == LOCUS$chr & V2 < LOCUS$end & V3 > LOCUS$start]
    locus_cn <- if (nrow(hit)) max(hit[[ncol(hit)]]) else 2
  }
  ## EpiTracer episomal call
  episomal <- FALSE
  for (g in gp) {
    inp <- tryCatch(aa_to_epitracer_inputs(g, sub("_graph\\.txt$", "_cycles.txt", g), sid,
                                           cnv_bed = cnvb), error = function(e) NULL)
    if (is.null(inp) || !length(inp$cnv_gr)) next
    r <- tryCatch(call_simple_excision(inp$ecdna_gr, inp$breakpoints_gr, inp$cnv_gr, onc),
                  error = function(e) NULL)
    if (!is.null(r) && nrow(as.data.frame(r)) &&
        any(as.data.frame(r)$episomal == "TRUE")) episomal <- TRUE
  }
  data.table(sample = sid, locus_maxCN = locus_cn, episomal = episomal)
}

arms <- list(`loxP-Cre (engineered ecDNA)` = "loxp15_cre_%dw",
             `p53-Cre (control)`           = "p53_cre_%dw")
res <- rbindlist(lapply(names(arms), function(a)
  rbindlist(lapply(1:5, function(w) {
    d <- run_one(sprintf(arms[[a]], w)); d[, `:=`(arm = a, week = w)]; d
  }))))

dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "validation/output/engineered_timecourse.tsv", sep = "\t")
print(res[order(arm, week), .(arm, week, locus_maxCN = round(locus_maxCN, 1), episomal)])

res[, call := ifelse(episomal, "EpiTracer: episomal", "not episomal")]
p <- ggplot(res, aes(week, locus_maxCN, group = arm)) +
  geom_line(aes(colour = arm), linewidth = 0.8) +
  geom_point(aes(colour = arm, shape = call, fill = call), size = 4, stroke = 1.1) +
  scale_shape_manual(values = c("EpiTracer: episomal" = 21, "not episomal" = 1), name = NULL) +
  scale_fill_manual(values = c("EpiTracer: episomal" = "#0f7d78", "not episomal" = "white"), name = NULL) +
  scale_colour_manual(values = c("loxP-Cre (engineered ecDNA)" = "#0f7d78",
                                 "p53-Cre (control)" = "grey55"), name = NULL) +
  scale_x_continuous(breaks = 1:5, labels = paste0(1:5, "w")) +
  labs(title = "EpiTracer detects the engineered Myc/Pvt1 episome as it forms",
       subtitle = "Pradella et al. mouse NPC time course (mm10, public AmpliconArchitect); filled = EpiTracer episomal call",
       x = "weeks after Cre induction", y = "max copy number at Myc/Pvt1 locus (chr15)") +
  theme_bw(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "top", legend.box = "vertical",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey30"))
ggsave("validation/output/engineered_timecourse.png", p, width = 8.4, height = 5.2, dpi = 130, bg = "white")
message("Wrote validation/output/engineered_timecourse.{tsv,png}")
