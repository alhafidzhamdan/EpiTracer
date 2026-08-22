#' Chromosome lengths for a bundled genome
#'
#' Returns the length (bp) of each chromosome from the bundled cytoband data (the
#' maximum band coordinate per chromosome). Used by [call_simple_excision()] to
#' test telomere proximity for the breakage-fusion-bridge (BFB) annotation: a
#' p-telomere sits at coordinate ~0 and a q-telomere at ~the chromosome length.
#'
#' @param genome One of `"hg38"`, `"hg19"`, `"mm10"`.
#' @return A named numeric vector of chromosome lengths (names are `seqnames`,
#'   e.g. `"chr7"`). Empty if the genome's cytoband data is unavailable.
#' @examples
#' load_chrom_lengths("hg38")["chr7"]
#' @export
load_chrom_lengths <- function(genome = c("hg38", "hg19", "mm10")) {
  genome <- match.arg(genome)
  f <- system.file("extdata", paste0("chr_info_", genome, ".rds"), package = "EpiTracer")
  if (!nzchar(f)) return(stats::setNames(numeric(0), character(0)))
  ci <- as.data.frame(readRDS(f))
  chrs <- unique(as.character(ci$seqnames))
  stats::setNames(vapply(chrs, function(ch) max(ci$end[as.character(ci$seqnames) == ch]),
                         numeric(1)), chrs)
}
