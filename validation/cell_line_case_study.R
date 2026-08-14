## ---------------------------------------------------------------------------
## EpiTracer validation — known-ecDNA cell-line case study (TEMPLATE)
##
## Ground-truth demonstration: run EpiTracer on cell lines whose ecDNA status at
## a focal oncogene has been established orthogonally (FISH / metaphase imaging /
## AmpliconArchitect), and confirm it (a) flags the amplicon as episomal and
## (b) reconstructs the circle. Good candidates:
##   * EGFR  — e.g. GBM39, PC3 (glioblastoma / prostate ecDNA lines)
##   * MYCN  — e.g. IMR-32, Kelly, NGP (neuroblastoma ecDNA lines)
##   * MYC   — e.g. COLO320-DM (double-minute line)
##
## Fill in the file paths below (PURPLE SV/CN + AmpliconArchitect amplicons for
## each line, coerced to the columns documented in ?call_episomal_ecdna) and run.
## This script is a template — it does not run as-is without those inputs.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer)
  library(GenomicRanges)
})

## ---- 1. Inputs (EDIT THESE) ------------------------------------------------
## Each object is a GRanges with the metadata columns listed in
## ?call_episomal_ecdna. `WGS_ID` / `sample` identify the cell line.
##
## ecdna_gr        <- readRDS("data/cell_lines/ecdna_amplicons.rds")   # AmpliconArchitect
## breakpoints_gr  <- readRDS("data/cell_lines/sv_breakpoints.rds")    # PURPLE (or coerced)
## cnv_gr          <- readRDS("data/cell_lines/cn_segments.rds")       # PURPLE (or coerced)
## cancer_genes_gr <- readRDS("data/cell_lines/cancer_genes.rds")      # needs $gene
##
## truth <- data.frame(                       # orthogonal ground truth
##   WGS_ID = c("GBM39", "IMR32", "COLO320"),
##   gene   = c("EGFR",  "MYCN",  "MYC"),
##   ecdna  = TRUE                            # confirmed ecDNA by imaging/AA
## )

## ---- 2. Call episomal ecDNA ------------------------------------------------
## calls <- call_episomal_ecdna(ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
##                              ext = 1e7, mc.cores = 4)
##
## Per-amplicon summary: is each confirmed-ecDNA line called episomal?
## library(data.table)
## per_amp <- as.data.table(calls)[, .(
##   episomal      = any(episomal == "TRUE"),
##   excision_scar = any(has_excision_scar == "TRUE")
## ), by = .(WGS_ID, ID)]
## print(per_amp)

## ---- 3. Visualise each confirmed line --------------------------------------
## dir.create("validation/output/cell_lines", recursive = TRUE, showWarnings = FALSE)
## for (s in unique(truth$WGS_ID)) {
##   plot_sv_linear(
##     sample   = s,
##     cnv_data = as.data.frame(cnv_gr[cnv_gr$sample == s]),
##     sv_data  = sv_df_for(s),          # your PURPLE bedpe -> chrom1/start1/... frame
##     outdir   = "validation/output/cell_lines"
##   )
##   plot_sv_reconstruction(              # wave-by-wave amplicon reconstruction
##     sample = s, ...
##   )
## }

## ---- 4. Report -------------------------------------------------------------
## Tabulate, for the manuscript: line | oncogene | ecDNA (orthogonal) |
## EpiTracer episomal call | excision scar | figure. Every confirmed-episomal
## line should be recovered; note any that are ecDNA-but-not-episomal (e.g.
## chromothriptic hubs) — those are correctly NOT called episomal and make the
## point that EpiTracer resolves *mechanism*, not merely ecDNA vs not.

message("Template only — edit the input paths in section 1 before running.")
