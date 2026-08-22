#' Centromere regions for a bundled genome
#'
#' Returns the centromere (cytoband `acen`) span per chromosome from the bundled
#' cytoband data. Used by [call_simple_excision()] to flag and exclude
#' "chromosomal bridge" events -- amplicons whose amplified edge terminates
#' *within* a centromere (a dicentric / bridge signature) rather than forming a
#' self-contained episome.
#'
#' @param genome One of `"hg38"`, `"hg19"`, `"mm10"`.
#' @return A [GenomicRanges::GRanges] with one range per chromosome (the
#'   centromere span). Empty if the genome's cytobands carry no `acen` bands.
#' @examples
#' load_centromeres("hg38")
#' @export
load_centromeres <- function(genome = c("hg38", "hg19", "mm10")) {
  genome <- match.arg(genome)
  f <- system.file("extdata", paste0("chr_info_", genome, ".rds"), package = "EpiTracer")
  if (!nzchar(f)) return(GenomicRanges::GRanges())
  ci <- as.data.frame(readRDS(f))
  acen <- ci[ci$gieStain == "acen", , drop = FALSE]
  if (!nrow(acen)) return(GenomicRanges::GRanges())
  chrs <- unique(as.character(acen$seqnames))
  gr <- do.call(c, lapply(chrs, function(ch) {
    b <- acen[as.character(acen$seqnames) == ch, ]
    GenomicRanges::GRanges(ch, IRanges::IRanges(min(b$start), max(b$end)))
  }))
  gr
}
