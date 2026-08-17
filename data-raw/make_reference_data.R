## ---------------------------------------------------------------------------
## Generate the bundled karyotype + oncogene references for each genome build.
##
## Produces, in inst/extdata/:
##   chr_info_{hg38,hg19,mm10}.rds      -- cytoband ideogram (UCSC cytoBand)
##   oncogene_coord_{hg38,hg19,mm10}.bed -- the default oncogene panel per build
##
## The hg38 files ship as-is (generated previously); this script regenerates the
## hg19 and mm10 equivalents from the SAME oncogene panel, so `genome=` can be
## switched without the caller supplying custom references.
##
## Sources (all public REST, no local downloads needed):
##   * Cytobands  -- UCSC REST API   https://api.genome.ucsc.edu (cytoBand track)
##   * hg19 genes -- MyGene.info      genomic_pos_hg19 field
##   * mm10 genes -- UCSC REST API    mm10 refGene track (MyGene mouse is GRCm39)
##
## Run from the package root:  Rscript data-raw/make_reference_data.R
## ---------------------------------------------------------------------------

stopifnot(requireNamespace("jsonlite", quietly = TRUE),
          requireNamespace("curl", quietly = TRUE))
outdir <- "inst/extdata"

fetch <- function(url) {
  h <- curl::new_handle(useragent = "Mozilla/5.0")
  rawToChar(curl::curl_fetch_memory(url, handle = h)$content)
}
getj <- function(url) jsonlite::fromJSON(fetch(url), simplifyVector = FALSE)

## human oncogene panel (from the existing hg38 bed) -------------------------
panel <- utils::read.table(file.path(outdir, "oncogene_coord_hg38.bed"),
                           sep = "\t", stringsAsFactors = FALSE,
                           col.names = c("chr", "start", "end", "strand", "gene"))
human_genes <- panel$gene
strand_chr <- function(s) ifelse(s == 1 | s == "+", "+", "-")
main_chr <- function(build) if (build == "mm10") c(1:19, "X", "Y") else c(1:22, "X", "Y")

## ---- hg19 gene coordinates (MyGene genomic_pos_hg19) ----------------------
message("Fetching hg19 gene coordinates from MyGene.info ...")
hg19_ensembl37 <- function(g) {   # fallback for genes absent from MyGene hg19
  d <- tryCatch(getj(sprintf(
    "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/%s?content-type=application/json", g)),
    error = function(e) NULL)
  if (is.null(d) || is.null(d$seq_region_name)) return(NULL)
  data.frame(chr = paste0("chr", d$seq_region_name), start = d$start, end = d$end,
             strand = strand_chr(d$strand), gene = g, stringsAsFactors = FALSE)
}
hg19 <- do.call(rbind, lapply(human_genes, function(g) {
  d <- getj(sprintf("https://mygene.info/v3/query?q=symbol:%s&species=human&fields=genomic_pos_hg19", g))
  hits <- d$hits
  gp <- if (length(hits)) hits[[1]]$genomic_pos_hg19 else NULL
  if (is.null(gp)) return(hg19_ensembl37(g))
  if (!is.null(gp$chr)) gp <- list(gp)              # single record -> list
  gp <- Filter(function(e) as.character(e$chr) %in% main_chr("hg19"), gp)
  if (!length(gp)) return(hg19_ensembl37(g))
  e <- gp[[which.max(vapply(gp, function(x) x$end - x$start, numeric(1)))]]
  data.frame(chr = paste0("chr", e$chr), start = e$start, end = e$end,
             strand = strand_chr(e$strand), gene = g, stringsAsFactors = FALSE)
}))

## ---- mm10 gene coordinates (mouse orthologs; UCSC mm10 refGene) -----------
message("Resolving mouse orthologs + fetching mm10 refGene from UCSC ...")
special <- c(TP53 = "Trp53")
to_mouse <- function(g) if (g %in% names(special)) special[[g]] else
  paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g))))
mouse_sym <- vapply(human_genes, to_mouse, character(1))

## mouse chromosome per ortholog (MyGene mouse, for the chromosome only)
mouse_chr <- vapply(mouse_sym, function(m) {
  d <- getj(sprintf("https://mygene.info/v3/query?q=symbol:%s&species=mouse&fields=genomic_pos", m))
  if (!length(d$hits)) return(NA_character_)
  gp <- d$hits[[1]]$genomic_pos
  if (!is.null(gp$chr)) gp <- list(gp)
  cc <- Filter(function(e) as.character(e$chr) %in% main_chr("mm10"), gp)
  if (length(cc)) as.character(cc[[1]]$chr) else NA_character_
}, character(1))

## exact mm10 coordinates from UCSC refGene, fetched once per chromosome
mm10 <- do.call(rbind, lapply(unique(stats::na.omit(mouse_chr)), function(ch) {
  tr <- getj(sprintf("https://api.genome.ucsc.edu/getData/track?genome=mm10;track=refGene;chrom=chr%s", ch))
  items <- tr[["refGene"]]; if (is.null(items)) items <- tr[[paste0("chr", ch)]]
  want <- mouse_sym[mouse_chr == ch & !is.na(mouse_chr)]
  do.call(rbind, lapply(want, function(m) {
    hit <- Filter(function(g) identical(g$name2, m), items)
    if (!length(hit)) return(NULL)
    e <- hit[[which.max(vapply(hit, function(x) x$txEnd - x$txStart, numeric(1)))]]
    data.frame(chr = paste0("chr", ch), start = e$txStart + 1L, end = e$txEnd,
               strand = e$strand, gene = m, stringsAsFactors = FALSE)
  }))
}))

## ---- cytoband ideograms (UCSC cytoBand) -----------------------------------
gie_color <- c(gneg = "white", gpos25 = "grey75", gpos33 = "grey66",
               gpos50 = "grey50", gpos66 = "grey33", gpos75 = "grey25",
               gpos100 = "grey0", acen = "red", gvar = "white", stalk = "white")
build_chr_info <- function(build) {
  message("Fetching ", build, " cytoBand from UCSC ...")
  d <- getj(sprintf("https://api.genome.ucsc.edu/getData/track?genome=%s;track=cytoBand", build))$cytoBand
  chroms <- paste0("chr", main_chr(build))
  rows <- do.call(rbind, lapply(chroms, function(ch) {
    b <- d[[ch]]; if (is.null(b)) return(NULL)
    do.call(rbind, lapply(b, function(x) data.frame(
      seqnames = ch, start = x$chromStart, end = x$chromEnd,
      name = x$name, gieStain = x$gieStain, stringsAsFactors = FALSE)))
  }))
  rows$width  <- rows$end - rows$start + 1L
  rows$strand <- "*"
  rows$y <- 1L
  col <- gie_color[rows$gieStain]; col[is.na(col)] <- "white"
  rows$color <- unname(col)
  ## rows are already in chromosome order (built per chrom in `chroms` order);
  ## keep seqnames as character to match the hg38 reference exactly.
  rows[, c("seqnames", "start", "end", "width", "strand", "name", "gieStain", "y", "color")]
}

## ---- write ----------------------------------------------------------------
write_bed <- function(df, path) utils::write.table(
  df[order(factor(df$chr, levels = paste0("chr", c(1:22, "X", "Y"))), df$start),
     c("chr", "start", "end", "strand", "gene")],
  path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write_bed(hg19, file.path(outdir, "oncogene_coord_hg19.bed"))
write_bed(mm10, file.path(outdir, "oncogene_coord_mm10.bed"))
saveRDS(build_chr_info("hg19"), file.path(outdir, "chr_info_hg19.rds"))
saveRDS(build_chr_info("mm10"), file.path(outdir, "chr_info_mm10.rds"))

message(sprintf("Done. hg19 genes: %d/%d, mm10 genes: %d/%d.",
                nrow(hg19), length(human_genes), nrow(mm10), length(human_genes)))
