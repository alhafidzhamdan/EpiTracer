# App-specific helpers. The input readers/reformatters now live in the EpiTracer
# package itself (read_cnv, read_sv, read_purple_sv_vcf, prepare_amplicon_inputs);
# this file keeps only the mechanism-table aggregation used by the UI.

suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(EpiTracer)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a) || (is.character(a) && !nzchar(a[1]))) b else a
.pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

onc_gr <- local({
  f <- system.file("extdata", "oncogene_coord_hg38.bed", package = "EpiTracer")
  if (nzchar(f)) {
    o <- utils::read.table(f, sep = "\t", col.names = c("chr","start","end","strand","gene"))
    GRanges(o$chr, IRanges(o$start, o$end), gene = o$gene)
  } else GRanges()
})

seeds_for_sample <- function(inp) {
  tryCatch(detect_amplicon_seeds(inp$cnv_gr, breakpoints = inp$breakpoints_gr),
           error = function(e) GRanges())
}

## run every mechanism caller on one sample's inputs and return a tidy amplicon table
call_mechanisms <- function(inp, genome = "hg38") {
  args <- list(ecdna_gr = NULL, breakpoints_gr = inp$breakpoints_gr,
               cnv_gr = inp$cnv_gr, cancer_genes_gr = onc_gr)
  cen <- tryCatch(load_centromeres(genome), error = function(e) NULL)
  cl  <- tryCatch(load_chrom_lengths(genome), error = function(e) NULL)
  epi <- tryCatch(do.call(call_simple_excision, c(args, list(centromeres = cen))), error = function(e) NULL)
  if (is.null(epi) || !nrow(as.data.frame(epi))) return(data.table())
  de <- as.data.table(as.data.frame(epi))
  agg <- de[, .(episomal = any(episomal == "TRUE"),
                excision_scar = any(has_excision_scar == "TRUE"),
                junction_homology = na.omit(junction_homology_class)[1],
                genes = paste(unique(stats::na.omit(gene)), collapse = ",")), by = .(WGS_ID, ID)]
  add <- function(tbl, f) if (!is.null(tbl) && nrow(tbl)) merge(agg, f(as.data.table(as.data.frame(tbl))),
                                                                by = c("WGS_ID","ID"), all.x = TRUE) else agg
  agg <- add(tryCatch(do.call(call_brf, args), error=function(e)NULL),
             function(d) d[, .(brf = any(brf=="TRUE"), parallel_pairs = max(n_parallel_pairs)), by=.(WGS_ID,ID)])
  agg <- add(tryCatch(do.call(call_micronucleation, args), error=function(e)NULL),
             function(d) d[, .(micronucleation = any(micronucleation=="TRUE")), by=.(WGS_ID,ID)])
  agg <- add(tryCatch(do.call(call_bfb, c(args, list(centromeres=cen, chrom_lengths=cl))), error=function(e)NULL),
             function(d) d[, .(bfb = any(bfb=="TRUE"), foldbacks = max(n_foldbacks)), by=.(WGS_ID,ID)])
  agg <- add(tryCatch(do.call(call_translocation_bridge_amp, c(args, list(centromeres=cen, chrom_lengths=cl))), error=function(e)NULL),
             function(d) d[, .(tba = any(tba=="TRUE"), tb_confident = any(tb_confident=="TRUE"),
                               tb_partner = paste(unique(tb_partner_chr[tb_partner_chr!=""]), collapse=";")), by=.(WGS_ID,ID)])
  agg[, mechanism := fifelse(episomal & micronucleation, "candidate_episomal+micronucleation",
                     fifelse(micronucleation, "micronucleation/chromothripsis",
                     fifelse(episomal, "episomal_ecDNA", "unclassified")))]
  agg[]
}
