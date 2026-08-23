## ---------------------------------------------------------------------------
## Whole-genome circos for every called chromoplexy cycle, with the cycle's own
## junctions highlighted (bold) and the rest of the sample's rearrangements
## greyed. One page per cycle in a single PDF. Circos is the natural view for
## chromoplexy: the closed chr_a -> chr_b -> chr_c -> chr_a ring is drawn as a
## bold loop across the inter-chromosomal link space.
##
## USAGE (from the package root):
##   Rscript validation/plot_chromoplexy_circos.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(EpiTracer); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"

sv <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))
d  <- fread("validation/output/cohort_chromoplexy_per_chr.tsv")
setorder(d, -n_junctions, -n_chromosomes, WGS_ID)
message("Rendering ", nrow(d), " chromoplexy cycles as highlighted circos")

outpdf <- "validation/output/chromoplexy_circos_montage.pdf"
grDevices::pdf(outpdf, onefile = TRUE)
for (i in seq_len(nrow(d))) {
  s  <- d$WGS_ID[i]
  ev <- strsplit(d$events[i], ",")[[1]]
  message(sprintf("  %-14s %-22s %d junctions", s, d$chromosomes[i], d$n_junctions[i]))
  tryCatch(
    plot_sv_circos(s, sv, cn, genome = "hg38",
                   highlight_events = ev, dim_unhighlighted = TRUE, outdir = NULL),
    error = function(e) message("    failed: ", conditionMessage(e)))
}
grDevices::dev.off()
message("Wrote ", outpdf)
