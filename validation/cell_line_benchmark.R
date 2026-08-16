## ---------------------------------------------------------------------------
## EpiTracer validation — public known-ecDNA cell-line benchmark
##
## Runs EpiTracer directly on public AmpliconArchitect (AA) reconstructions and
## scores it against orthogonal ground truth. Two evidence tiers:
##
##   GOLD  — cell lines whose ecDNA status at a focal oncogene is established by
##           metaphase FISH (validation/cell_line_panel.tsv). Measures whether
##           EpiTracer detects the amplicon on ecDNA+ lines and, crucially, does
##           NOT flag episomal on ecDNA-negative (HSR / no-amplicon) lines.
##
##   SILVER — AmpliconClassifier's own ecDNA+/BFB+ call on every downloaded
##            amplicon. Measures detection concordance at scale and reports the
##            episomal fraction among AC-called ecDNA (the mechanistic layer AC
##            does not provide).
##
## DATA (see validation/README.md, "Public cell-line benchmark"):
##   AA outputs for 329 CCLE cell lines are public on AmpliconRepository
##   (https://ampliconrepository.org, CC-BY-4.0). Download per-sample archives
##   into a directory, one sub-folder per cell line, each containing the AA
##   *_graph.txt / *_cycles.txt and the AmpliconClassifier
##   *_amplicon_classification_profiles.tsv.
##
## USAGE (from the package root):
##   Rscript validation/cell_line_benchmark.R selftest          # parser + caller smoke test (no download)
##   Rscript validation/cell_line_benchmark.R /path/to/aa_data  # real benchmark
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer)
  library(GenomicRanges)
  library(data.table)
})
source("validation/aa_to_epitracer.R")

## --- minimal focal-oncogene coordinates (hg38) ------------------------------
oncogenes_gr <- GRanges(
  c("chr8", "chr7", "chr2", "chr10"),
  IRanges(c(127735434, 55019017, 15940550, 121478330),
          c(127742951, 55211628, 15947007, 121598458)),
  gene = c("MYC", "EGFR", "MYCN", "FGFR2"))

## --- run EpiTracer on one AA amplicon folder --------------------------------
call_one_sample <- function(sample_dir, sample_id, flip_dup_del = FALSE) {
  ## recursive: AmpliconRepository puts graphs in <sample>_reconstruction_results/
  gpaths <- list.files(sample_dir, pattern = "_amplicon[0-9]+_graph\\.txt$",
                       full.names = TRUE, recursive = TRUE)
  ## genome-wide CN (recommended) — AmpliconSuite CNVkit calls, if present
  cnv_bed <- list.files(sample_dir, pattern = "(CNV_CALLS|\\.cnv)\\.bed$",
                        full.names = TRUE, recursive = TRUE)
  cnv_bed <- if (length(cnv_bed)) cnv_bed[1] else NULL
  out <- lapply(gpaths, function(gp) {
    cp <- sub("_graph\\.txt$", "_cycles.txt", gp)
    if (!file.exists(cp)) return(NULL)
    inp <- aa_to_epitracer_inputs(gp, cp, sample_id, flip_dup_del = flip_dup_del,
                                  cnv_bed = cnv_bed)
    if (!length(inp$cnv_gr)) return(NULL)
    res <- tryCatch(
      call_episomal_ecdna(ecdna_gr = inp$ecdna_gr,
                          breakpoints_gr = inp$breakpoints_gr,
                          cnv_gr = inp$cnv_gr,
                          cancer_genes_gr = oncogenes_gr),
      error = function(e) { message("  ! ", sample_id, ": ", conditionMessage(e)); NULL })
    if (is.null(res) || !nrow(as.data.frame(res))) return(NULL)
    d <- as.data.table(as.data.frame(res))
    d[, amplicon := basename(gp)]
    d[]
  })
  rbindlist(Filter(Negate(is.null), out), fill = TRUE)
}

## --- per-amplicon episomal summary ------------------------------------------
summarise_calls <- function(calls) {
  if (!nrow(calls)) return(calls)
  calls[, .(
    episomal      = any(episomal == "TRUE"),
    excision_scar = any(has_excision_scar == "TRUE"),
    max_vf        = suppressWarnings(max(VF, na.rm = TRUE))
  ), by = .(WGS_ID, amplicon)]
}

## ===========================================================================
## SELF-TEST — verifies the parser + caller path on the bundled fixture
## ===========================================================================
run_selftest <- function() {
  fx <- "validation/fixtures/COLO320DM"
  stopifnot(dir.exists(fx))
  message("[selftest] parsing AA fixture ...")
  inp <- aa_to_epitracer_inputs(
    file.path(fx, "COLO320DM_amplicon1_graph.txt"),
    file.path(fx, "COLO320DM_amplicon1_cycles.txt"),
    "COLO320DM")

  stopifnot(length(inp$cnv_gr) == 3L)
  stopifnot("DUP" %in% inp$breakpoints_gr$svclass)
  ## the highest-read-support junction must be the DUP (convention sanity check)
  top <- inp$breakpoints_gr[which.max(inp$breakpoints_gr$VF)]
  message(sprintf("[selftest] top junction: svclass=%s VF=%s (expect DUP)",
                  top$svclass, top$VF))
  stopifnot(identical(as.character(top$svclass), "DUP"))

  res <- call_episomal_ecdna(ecdna_gr = inp$ecdna_gr,
                             breakpoints_gr = inp$breakpoints_gr,
                             cnv_gr = inp$cnv_gr,
                             cancer_genes_gr = oncogenes_gr)
  d <- as.data.table(as.data.frame(res))
  message("[selftest] caller returned ", nrow(d), " rows; ",
          "episomal=", paste(unique(d$episomal), collapse = ","),
          " scar=", paste(unique(d$has_excision_scar), collapse = ","))
  stopifnot(any(d$episomal == "TRUE"))
  stopifnot(any(d$has_excision_scar == "TRUE"))
  message("[selftest] PASS — parser + episome heuristic run end-to-end on AA-format input.")
  invisible(TRUE)
}

## ===========================================================================
## REAL BENCHMARK
## ===========================================================================
run_benchmark <- function(data_root, flip_dup_del = FALSE) {
  panel <- fread("validation/cell_line_panel.tsv")
  ## AmpliconRepository archives extract to results/samples/<SAMPLE>/ — use that
  ## if present, otherwise treat each immediate sub-folder as a sample.
  samples_root <- if (dir.exists(file.path(data_root, "results", "samples")))
    file.path(data_root, "results", "samples") else data_root
  sample_dirs <- list.dirs(samples_root, recursive = FALSE)
  if (!length(sample_dirs)) stop("no sample sub-folders under ", samples_root)
  message("Found ", length(sample_dirs), " sample folders.")

  ## --- run EpiTracer on every sample ---
  all_calls <- rbindlist(lapply(sample_dirs, function(sd) {
    sid <- basename(sd)
    message("EpiTracer: ", sid)
    call_one_sample(sd, sid, flip_dup_del = flip_dup_del)
  }), fill = TRUE)
  epi <- summarise_calls(all_calls)

  ## --- SILVER: AmpliconClassifier labels (consolidated file, or per-sample) ---
  ac_file <- list.files(data_root, pattern = "_amplicon_classification_profiles\\.tsv$",
                        full.names = TRUE, recursive = TRUE)
  if (length(ac_file)) {
    ac <- as.data.table(read_ac_profiles(ac_file[1]))
    if ("sample_name" %in% names(ac)) ac[, WGS_ID := sample_name]
  } else {
    ac <- data.table()
  }

  dir.create("validation/output", showWarnings = FALSE, recursive = TRUE)
  fwrite(epi, "validation/output/cell_line_episomal_calls.tsv", sep = "\t")
  if (nrow(ac)) fwrite(ac, "validation/output/cell_line_ac_labels.tsv", sep = "\t")

  ## --- SILVER concordance (if AC ecDNA flag present) ---
  ac_flag <- grep("ecDNA", names(ac), value = TRUE)
  if (nrow(ac) && length(ac_flag)) {
    ac[, ac_ecdna := grepl("Positive", get(ac_flag[1]), ignore.case = TRUE)]
    silver <- ac[, .(ac_ecdna = any(ac_ecdna)), by = WGS_ID][
      epi[, .(episomal = any(episomal)), by = WGS_ID], on = "WGS_ID"]
    message("\n--- SILVER: EpiTracer episomal vs AmpliconClassifier ecDNA+ (per sample) ---")
    print(table(AC_ecDNA = silver$ac_ecdna, EpiTracer_episomal = silver$episomal, useNA = "ifany"))
  }

  ## --- GOLD: FISH-validated panel ---
  gold <- merge(panel, epi[, .(episomal = any(episomal)), by = WGS_ID],
                by.x = "ccle_name", by.y = "WGS_ID", all.x = FALSE)
  if (nrow(gold)) {
    message("\n--- GOLD: FISH ground truth vs EpiTracer episomal ---")
    print(gold[, .(cell_line, oncogene, ecdna_status, EpiTracer_episomal = episomal)])
    pos <- gold[ecdna_status == "positive"]
    neg <- gold[ecdna_status == "negative"]
    message(sprintf("Detection on ecDNA+ lines: %d/%d called episomal",
                    sum(pos$episomal, na.rm = TRUE), nrow(pos)))
    message(sprintf("False-episomal on ecDNA- lines: %d/%d",
                    sum(neg$episomal, na.rm = TRUE), nrow(neg)))
  } else {
    message("\n(No panel cell lines matched the downloaded CCLE names — check ccle_name mapping.)")
  }

  message("\nWrote validation/output/cell_line_{episomal_calls,ac_labels}.tsv")
  invisible(list(epi = epi, ac = ac))
}

## --- entry point ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || args[1] == "selftest") {
  run_selftest()
} else {
  run_benchmark(args[1], flip_dup_del = isTRUE("flip" %in% args))
}
