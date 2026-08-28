## ---------------------------------------------------------------------------
## Montage of ALL episomal-catalogue amplicons, each labelled with its inferred
## FORMATION MECHANISM from the founder-boundary analysis
## (validation/nominate_founder_boundary.R):
##   simple excision (episome)   -- DUP founder (copy-gaining circularisation)
##   two-ecDNA recombination     -- DUP founder + amplified inter-locus TRA
##   inverted duplication        -- inversion founder (not simple excision)
##   unresolved                  -- no copy-gaining founder
## Amplicon zoom, low-VF passenger TRAs dropped; the panel title is coloured by
## mechanism.
##
## USAGE:
##   Rscript validation/plot_all_amplicons_montage.R \
##       example_data/all_348_SV_bedpe.rds example_data/all_353_CN_segments.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(EpiTracer); library(data.table); library(GenomicRanges); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sv_rds <- if (length(args) >= 1) args[1] else "example_data/all_348_SV_bedpe.rds"
cn_rds <- if (length(args) >= 2) args[2] else "example_data/all_353_CN_segments.rds"
sv_raw <- as.data.table(readRDS(sv_rds)); cn <- as.data.table(readRDS(cn_rds))

nom <- fread("validation/output/founder_nomination.tsv")
nom[, mech := fifelse(founder_class %in% c("h2hINV", "t2tINV"), "inverted duplication",
             fifelse(founder_class == "DUP" & interlocus_tra == TRUE & max_tra_vf >= 100, "two-ecDNA recombination",
             fifelse(founder_class == "DUP", "simple excision (episome)", "unresolved")))]
mech_col <- c("simple excision (episome)" = "#1f6fb4", "two-ecDNA recombination" = "#6a51a3",
              "inverted duplication" = "#c2181b", "unresolved" = "#f0a848")
nom[, mech_ord := match(mech, names(mech_col))]
setorder(nom, mech_ord, gene, -founder_vf)
message("amplicons to plot: ", nrow(nom))

pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
nochr <- function(x) sub("^chr", "", x)
onc <- utils::read.table(system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer"),
                         sep = "\t", col.names = c("chr", "start", "end", "strand", "gene"))
cancer_genes_gr <- GRanges(onc$chr, IRanges(onc$start, onc$end), gene = onc$gene)
gene_coord <- data.frame(chr = onc$chr, start = onc$start, end = onc$end, strand = onc$strand, gene = onc$gene)
cl <- load_chrom_lengths("hg38")

build <- function(s) {
  cs <- cn[sample == s]; vs <- sv_raw[sample == s]; if (!nrow(vs) || !nrow(cs)) return(NULL)
  base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6])))), by = seqnames]
  pcol <- base$b[match(cs$seqnames, base$seqnames)]; pcol[is.na(pcol)] <- 2; cs$ploidy <- pcol
  cnv_gr <- GRanges(cs$seqnames, IRanges(cs$start, cs$end), sample = s, copyNumber = cs$copyNumber,
                    ploidy = pcol, majorAlleleCopyNumber = cs$majorAlleleCopyNumber,
                    minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  ev <- as.character(vs$name)
  mk <- function(ch, pos, st) data.frame(chr = pfx(ch), pos = pos, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = st, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN, VF = vs$VF,
    insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen))
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GRanges(bp$chr, IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID, event = bp$event,
    svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF, PURPLE_JCN = bp$PURPLE_JCN,
    VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GRanges(bp$chr, IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4)); ov <- findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[subjectHits(ov)], queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, bgr = bgr, cs = cs, vs = vs)
}

## keep panels tagged with mechanism so we can order the montage
panels <- list()
for (s in unique(nom$WGS_ID)) {
  inp <- build(s); if (is.null(inp)) next
  ecdna_gr <- tryCatch(detect_amplicon_seeds(inp$cnv_gr, min_cn_ratio = 3, gap = 1e6,
                          min_width = 1e5, breakpoints = inp$bgr), error = function(e) NULL)
  if (is.null(ecdna_gr) || !length(ecdna_gr)) next
  cs <- inp$cs; vs <- inp$vs
  for (i in which(nom$WGS_ID == s)) {
    id <- nom$ID[i]
    seed <- ecdna_gr[ecdna_gr$ID == id]; if (!length(seed)) next
    best <- NULL
    for (ch in unique(as.character(seqnames(seed)))) {
      amp <- cs[seqnames == ch & copyNumber > 3 * ploidy &
                  end >= min(start(seed)[as.character(seqnames(seed)) == ch]) &
                  start <= max(end(seed)[as.character(seqnames(seed)) == ch])]
      if (!nrow(amp)) next
      w <- max(amp$end) - min(amp$start)
      if (is.null(best) || w > best$w) best <- list(ch = ch, lo = min(amp$start), hi = max(amp$end), w = w)
    }
    if (is.null(best)) next
    ch <- best$ch; pc <- nochr(ch); lo <- best$lo; hi <- best$hi
    touch <- vs[(chrom1 == pc & start1 >= lo & start1 <= hi) | (chrom2 == pc & start2 >= lo & start2 <= hi)]
    maxvf <- if (nrow(touch)) max(touch$VF, na.rm = TRUE) else 0; thr <- 0.1 * maxvf
    vs_keep <- vs[!(svclass == "TRA" & VF < thr)]
    kt <- vs_keep[(chrom1 == pc & start1 >= lo & start1 <= hi) | (chrom2 == pc & start2 >= lo & start2 <= hi)]
    pos <- c(lo, hi, kt$start1[kt$chrom1 == pc], kt$start2[kt$chrom2 == pc])
    wlo <- min(pos); whi <- max(pos); pad <- max(2e6, 0.12 * (whi - wlo))
    cr <- matrix(c(max(1, wlo - pad), min(as.numeric(cl[[ch]]), whi + pad)), nrow = 1)

    gg <- if (nchar(nom$gene[i])) nom$gene[i] else NA
    m <- nom$mech[i]
    ttl <- sprintf("%s | %s | %s | %s (founder %s%s)",
                   s, ch, ifelse(is.na(gg) | gg == "", "(no oncogene)", gg), m,
                   nom$founder_class[i],
                   ifelse(nom$interlocus_tra[i], sprintf(", TRA VF %.0f", nom$max_tra_vf[i]), ""))
    p <- tryCatch(
      plot_sv_linear(sample = s, cnv_data = as.data.frame(cs[seqnames == ch]),
                     sv_data = as.data.frame(vs_keep), genome = "hg38",
                     chromosome = ch, chromosome_range = cr, events = c("amp", "loh"),
                     genes_to_highlight = if (is.na(gg)) NULL else unlist(strsplit(gg, ",")),
                     gene_coord = gene_coord, save = FALSE, verbose = FALSE),
      error = function(e) { message("  panel failed ", s, "/", id, ": ", conditionMessage(e)); NULL })
    if (!is.null(p)) {
      panels[[length(panels) + 1]] <- list(
        plot = p + ggtitle(ttl) + theme(plot.title = element_text(size = 7.5, face = "bold",
                                          colour = mech_col[[m]])),
        ord = nom$mech_ord[i], gene = gg, vf = nom$founder_vf[i])
      message(ttl)
    }
  }
}
plots <- lapply(panels[order(vapply(panels, function(x) x$ord, numeric(1)))], `[[`, "plot")
message("panels: ", length(plots))

if (requireNamespace("patchwork", quietly = TRUE) && length(plots)) {
  ncol <- 4; nrow <- 4; per <- ncol * nrow
  outpdf <- "validation/output/all_amplicons_mechanism_montage.pdf"
  grDevices::cairo_pdf(outpdf, width = 8 * ncol, height = 4.2 * nrow, onefile = TRUE)
  for (k in seq(1, length(plots), by = per)) {
    grp <- plots[k:min(k + per - 1, length(plots))]
    pg <- patchwork::wrap_plots(grp, ncol = ncol, nrow = nrow) +
      patchwork::plot_annotation(
        title = sprintf("Amplicon formation mechanisms %d-%d of %d  |  blue=simple excision  purple=two-ecDNA recombination  red=inverted duplication  amber=unresolved",
                        k, min(k + per - 1, length(plots)), length(plots)),
        theme = theme(plot.title = element_text(face = "bold", size = 12)))
    print(pg)
  }
  grDevices::dev.off()
  message("Wrote ", outpdf)
} else message("no panels")
