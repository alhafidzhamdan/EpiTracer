## ---------------------------------------------------------------------------
## EpiTracer validation — per-locus ecDNA mechanism benchmark across CCLE
##
## Runs the per-locus classifier (classify_loci.R) on every amplicon of every
## CCLE cell line and tabulates the mechanism of each major amplified locus
## (episomal / BFB / chimeric-translocation / complex). This supersedes the
## per-amplicon episomal count with a per-locus one, so a chimeric amplicon that
## fuses an episome to a BFB (e.g. NCI-H2170) contributes an episomal locus AND
## a non-episomal one, rather than a single misleading "episomal" call.
##
## USAGE (from the package root):
##   Rscript validation/cell_line_loci_benchmark.R /path/to/results/samples
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(GenomicRanges); library(data.table)
})
source("validation/aa_to_epitracer.R")
source("validation/classify_loci.R")

args <- commandArgs(trailingOnly = TRUE)
samples_root <- if (length(args)) args[1] else stop("pass the CCLE results/samples folder")

## hg38 oncogene panel, to annotate each locus
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
onc_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
label_genes <- function(chr, lo, up) {
  h <- onc[onc$chr == chr & onc$start <= up & onc$end >= lo, "gene"]
  if (length(h)) paste(h, collapse = ",") else ""
}

sample_dirs <- list.dirs(samples_root, recursive = FALSE)
message("Classifying loci across ", length(sample_dirs), " samples ...")

all_loci <- rbindlist(lapply(sample_dirs, function(sd) {
  wgs <- basename(sd)
  cnv <- list.files(sd, "CNV_CALLS\\.bed$", full.names = TRUE, recursive = TRUE)
  cnv <- if (length(cnv)) cnv[1] else NULL
  graphs <- list.files(sd, "_amplicon[0-9]+_graph\\.txt$", full.names = TRUE, recursive = TRUE)
  rbindlist(lapply(graphs, function(g) {
    amp <- sub(".*_(amplicon[0-9]+)_.*", "\\1", basename(g))
    tryCatch(classify_loci(g, wgs, amp, cnv_bed = cnv), error = function(e) NULL)
  }), fill = TRUE)
}), fill = TRUE)

if (nrow(all_loci)) all_loci[, oncogenes := mapply(label_genes, chr, start, end)]
dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
fwrite(all_loci, "validation/output/cell_line_loci.tsv", sep = "\t")

cat("\n=== major amplified loci by mechanism (across all CCLE amplicons) ===\n")
print(all_loci[, .N, by = mechanism][order(-N)])

epi <- all_loci[mechanism == "episomal"]
cat("\n=== EPISOMAL loci:", nrow(epi), "in", uniqueN(epi$sample), "cell lines",
    "(", sum(epi$fused), "fused ) ===\n")
print(epi[order(-max_cn), .(sample, chr, maxCN = round(max_cn),
      bDUP = round(boundary_dup_cn), foldbacks, fused, oncogenes)][1:min(30, nrow(epi))])
cat("\nWrote validation/output/cell_line_loci.tsv\n")
