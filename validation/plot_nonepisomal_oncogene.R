## ---------------------------------------------------------------------------
## Montage of every NON-EPISOMAL amplicon that carries a known oncogene -- the
## focal ecDNA candidates that FAILED the episomal (simple-excision) call.
##
## For each amplicon we draw its focal locus (the detected seed region, padded)
## on the seed chromosome, with amp + LOH events and the oncogene labelled. The
## panels are tiled into a paginated grid and written to one PDF.
##
## USAGE (from the package root):
##   Rscript validation/plot_nonepisomal_oncogene.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"

sv <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))
pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
sv[, `:=`(chrom1 = pfx(chrom1), chrom2 = pfx(chrom2))]
cn[, seqnames := pfx(as.character(seqnames))]
cl <- load_chrom_lengths("hg38")

d <- fread("validation/output/cohort_original_caller_per_chr.tsv")
d[, mech := fifelse(bfb == TRUE, "BFB", fifelse(tba == TRUE, "TBA",
        fifelse(micronucleation == TRUE, "micronucleation",
        fifelse(brf == TRUE, "BRF-only", "unclassified"))))]
ne <- d[episomal == FALSE & genes != ""]

## total SV burden per sample (one bedpe row = one SV)
sv_burden <- sv[, .(total_sv = .N), by = sample]
ne[sv_burden, total_sv := i.total_sv, on = c(WGS_ID = "sample")]
ne[is.na(total_sv), total_sv := 0L]

## Sort by fold-back-inversion (FBI) burden, then BRF burden -- most rearranged first.
setorder(ne, -n_foldbacks, -n_parallel_pairs, WGS_ID, ID)
message("Non-episomal oncogene amplicons: ", nrow(ne), " in ", uniqueN(ne$WGS_ID), " samples")

per_chr_ploidy <- function(cs) {
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  p <- base$b[match(cs$seqnames, base$seqnames)]; p[is.na(p)] <- 2; p
}

## per-sample cache of seeds + cnv_data (so the global FBI sort order is preserved)
cache <- new.env(parent = emptyenv())
sample_ctx <- function(s) {
  key <- make.names(s)
  if (!is.null(cache[[key]])) return(cache[[key]])
  cs <- cn[sample == s]
  gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s, copyNumber = cs$copyNumber,
    ploidy = per_chr_ploidy(cs), majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  sd <- tryCatch(detect_amplicon_seeds(gr), error = function(e) NULL)
  cnv_data <- as.data.frame(cs); cnv_data$ploidy <- per_chr_ploidy(cs)
  ctx <- list(cs = cs, sd = sd, cnv_data = cnv_data, vs = as.data.frame(sv[sample == s]))
  cache[[key]] <- ctx; ctx
}

plots <- list()
for (j in seq_len(nrow(ne))) {
  s <- ne$WGS_ID[j]; amp_id <- ne$ID[j]; gene <- ne$genes[j]; mech <- ne$mech[j]
  ctx <- sample_ctx(s); sd <- ctx$sd; if (is.null(sd) || !length(sd)) next
  hit <- which(sd$ID == amp_id | sd$ID == sub(paste0("^", s, "_"), "", amp_id))
  if (!length(hit)) { message("  no seed: ", amp_id); next }
  seed <- sd[hit[1]]; chr <- as.character(seqnames(seed)); w <- width(seed)
  lo <- max(1, start(seed) - max(1e6, w)); hi <- min(as.numeric(cl[[chr]]), end(seed) + max(1e6, w))
  ## peak copy number over the amplicon seed
  cnmax <- suppressWarnings(max(ctx$cs[seqnames == chr & end >= start(seed) & start <= end(seed)]$copyNumber))
  if (!is.finite(cnmax)) cnmax <- NA_real_
  gh <- strsplit(gene, ",")[[1]]
  ## mechanism tag block
  tag <- c(if (ne$tba[j])            sprintf("TBA(%d)", ne$n_boundary_tra[j]),
           if (ne$micronucleation[j]) "mnc",
           if (ne$brf[j])             "BRF",
           if (isTRUE(ne$excision_scar[j])) "scar")
  tag <- if (length(tag)) paste(tag, collapse = " ") else mech
  title <- sprintf("%s  %s  %s", s, gene, chr)
  sub   <- sprintf("FBI %d | BRF pairs %s | SVs %d | CNmax %s | %s",
                   ne$n_foldbacks[j], format(ne$n_parallel_pairs[j], big.mark = ","),
                   ne$total_sv[j], if (is.na(cnmax)) "NA" else round(cnmax), tag)
  p <- tryCatch(
    plot_sv_linear(sample = s, cnv_data = ctx$cnv_data, sv_data = ctx$vs, genome = "hg38",
                   chromosome = chr, chromosome_range = matrix(c(lo, hi), nrow = 1),
                   events = c("amp", "loh"), genes_to_highlight = gh, save = FALSE, verbose = FALSE),
    error = function(e) { message("  ", amp_id, " plot failed: ", conditionMessage(e)); NULL })
  if (!is.null(p)) plots[[length(plots) + 1]] <- p +
    ggplot2::labs(title = title, subtitle = sub) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 9, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey25"))
  message(sprintf("%-16s %-16s FBI=%d BRF=%d SV=%d", s, gene, ne$n_foldbacks[j], ne$n_parallel_pairs[j], ne$total_sv[j]))
}
message("panels: ", length(plots))

## paginated montage: 3 columns x 6 rows per page, one PDF
if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 3; nrow <- 6; per <- ncol * nrow
  outpdf <- "validation/output/nonepisomal_oncogene_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 7.5 * ncol, height = 3.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Non-episomal oncogene amplicons %d-%d of %d — sorted by fold-back-inversion burden (focal ecDNA candidates that failed the episomal call)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("patchwork not available; montage skipped")
