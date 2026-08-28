## ---------------------------------------------------------------------------
## Montage of every CHROMOTHRIPTIC episomal ecDNA in the cohort, one panel each.
##
## For each amplicon that the caller flags episomal (call_simple_excision) AND
## chromothriptic (call_chromothripsis; ShatterSeek-style hallmarks over the
## amplified footprint), draw the amplicon zoomed to its footprint, with its
## INTERNAL (shattering) structural variants highlighted bold red over the greyed
## rest of the sample's rearrangements. amp + LOH are shaded; the oncogene is
## labelled. Title carries the shattering statistics.
##
## USAGE (from the package root):
##   Rscript validation/plot_chromothripsis_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
gene_coord <- data.frame(chr = onc$chr, start = onc$start, end = onc$end,
                         strand = onc$strand, gene = onc$gene)
centromeres <- load_centromeres("hg38")
cl <- load_chrom_lengths("hg38")

## caller-side inputs (prefixed), identical to validation/call_episomal_cohort.R
build_inputs <- function(s) {
  cs <- cn[sample == s]; vs <- sv_raw[sample == s]
  if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  ploidy_col <- base$b[match(cs$seqnames, base$seqnames)]; ploidy_col[is.na(ploidy_col)] <- 2
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s,
                    copyNumber = cs$copyNumber, ploidy = ploidy_col,
                    majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
                    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- if ("name" %in% names(vs)) as.character(vs$name) else paste0("SV", seq_len(nrow(vs)))
  mk <- function(ch, pos, strand) data.frame(chr = pfx(ch), pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = strand, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN,
    VF = vs$VF, insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4))
  ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, breakpoints_gr = bgr)
}

per_chr_ploidy <- function(cs) {
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  p <- base$b[match(cs$seqnames, base$seqnames)]; p[is.na(p)] <- 2; p
}

samples <- intersect(unique(sv_raw$sample), unique(cn$sample))
if (length(args) >= 3) samples <- intersect(samples, strsplit(args[3], ",")[[1]])   # optional smoke-test subset
message("Scanning ", length(samples), " samples for chromothriptic ecDNA ...")

panels <- list()
for (s in samples) {
  inp <- build_inputs(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                                             min_width = 1e5, breakpoints = inp$breakpoints_gr),
                       error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  a <- list(ecdna_gr = ecdna_gr, breakpoints_gr = inp$breakpoints_gr,
            cnv_gr = inp$cnv_gr, cancer_genes_gr = cancer_genes_gr)
  epi <- tryCatch(do.call(call_simple_excision, c(a, list(centromeres = centromeres))), error = function(e) NULL)
  ct  <- tryCatch(do.call(call_chromothripsis, a), error = function(e) NULL)
  if (is.null(epi) || is.null(ct) || !nrow(as.data.frame(ct))) next
  de <- as.data.table(as.data.frame(epi))
  epi_ids <- de[, .(e = any(episomal == "TRUE"),
                    gene = paste(unique(na.omit(gene)), collapse = ",")), by = ID]
  cd <- as.data.table(as.data.frame(ct))
  ct_ids <- cd[, .(ct = any(chromothripsis == "TRUE"),
                   nsv = max(n_internal_sv), nosc = max(cn_oscillations),
                   pval = suppressWarnings(min(sv_type_pval, na.rm = TRUE))), by = ID]
  tgt <- merge(epi_ids[e == TRUE], ct_ids[ct == TRUE], by = "ID")
  if (!nrow(tgt)) next

  cs <- cn[sample == s]; cs$ploidy <- per_chr_ploidy(cs)
  for (i in seq_len(nrow(tgt))) {
    id <- tgt$ID[i]
    seed <- ecdna_gr[ecdna_gr$ID == id]
    chrs <- unique(as.character(GenomeInfoDb::seqnames(seed)))          # prefixed
    ## tight amplified footprint per chromosome (padded)
    cr <- list(); keep_chr <- character()
    for (ch in chrs) {
      amp <- cs[seqnames == ch & copyNumber > 3 * ploidy &
                  end >= min(start(seed)[as.character(seqnames(seed)) == ch]) &
                  start <= max(end(seed)[as.character(seqnames(seed)) == ch])]
      if (!nrow(amp)) next
      lo <- min(amp$start); hi <- max(amp$end); pad <- max(1e5, 0.15 * (hi - lo))
      cr[[ch]] <- c(max(1, lo - pad), min(as.numeric(cl[[ch]]), hi + pad)); keep_chr <- c(keep_chr, ch)
    }
    if (!length(keep_chr)) next
    crm <- do.call(rbind, cr[keep_chr])

    vs <- sv_raw[sample == s]
    cnv_data <- as.data.frame(cs[seqnames %in% keep_chr])
    gg <- if (nchar(tgt$gene[i])) tgt$gene[i] else NA
    ttl <- sprintf("%s | %s | %s | %d internal SVs, %d CN osc, joins p=%.2g",
                   s, paste(keep_chr, collapse = ","),
                   ifelse(is.na(gg) | gg == "", "(no oncogene)", gg),
                   tgt$nsv[i], tgt$nosc[i], tgt$pval[i])
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = cnv_data, sv_data = as.data.frame(vs), genome = "hg38",
                     chromosome = keep_chr, chromosome_range = crm, events = c("amp", "loh"),
                     genes_to_highlight = if (is.na(gg)) NULL else unlist(strsplit(gg, ",")),
                     gene_coord = gene_coord, save = FALSE, verbose = FALSE),
      error = function(e) { message("  panel failed ", s, "/", id, ": ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      panels[[length(panels) + 1]] <- list(
        plot = p + ggtitle(ttl) + theme(plot.title = element_text(size = 9, face = "bold")),
        gene = gg, nsv = tgt$nsv[i])
      message(ttl)
    }
  }
}
message("chromothriptic ecDNA panels: ", length(panels))

## order EGFR first, then by internal-SV burden
ord <- order(vapply(panels, function(x) !identical(x$gene, "EGFR"), logical(1)),
             -vapply(panels, function(x) x$nsv, numeric(1)))
plots <- lapply(panels[ord], `[[`, "plot")

if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 4; nrow <- 4; per <- ncol * nrow
  outpdf <- "validation/output/chromothripsis_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 8 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Chromothripsis within episomal ecDNA %d-%d of %d (ShatterSeek-style hallmarks)",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = theme(plot.title = element_text(face = "bold", size = 15)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("patchwork not available or no panels; montage skipped")
