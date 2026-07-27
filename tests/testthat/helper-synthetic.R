## Synthetic-data builders shared across tests.
## `flank_cn` controls the copy number of the segments flanking the amplicon:
## low (2) => non-gained flanks (episome); high (10) => gained flanks (not an
## episome, exercises the corrected boundary-flank logic).

make_episome_inputs <- function(flank_cn = 2) {
  suppressPackageStartupMessages({
    requireNamespace("GenomicRanges")
    requireNamespace("data.table")
  })

  ecdna_gr <- GenomicRanges::GRanges(
    "chr7", IRanges::IRanges(55000000, 55500000),
    ID = "S1_amp1", WGS_ID = "S1"
  )

  bp <- data.table::data.table(
    seqnames   = "chr7",
    start      = c(54900000, 55000000, 55500000, 55600000),
    end        = c(54900000, 55000000, 55500000, 55600000),
    WGS_ID     = "S1",
    event      = c("DEL1", "DUP1", "DUP1", "DEL1"),
    svclass    = c("DEL", "DUP", "DUP", "DEL"),
    PURPLE_AF  = c(0.4, 0.9, 0.9, 0.4),
    PURPLE_JCN = c(1, 40, 40, 1),
    VF         = c(200, 1000, 1000, 200),
    PURPLE_CN  = c(2, 50, 50, 2),
    insLen     = 0L,
    HOMLEN     = 0L
  )
  breakpoints_gr <- GenomicRanges::makeGRangesFromDataFrame(bp, keep.extra.columns = TRUE)

  cnv <- data.table::data.table(
    seqnames = "chr7",
    start    = c(40000000, 55000000, 55500001),
    end      = c(54999999, 55500000, 70000000),
    sample   = "S1",
    copyNumber = c(flank_cn, 50, flank_cn),
    ploidy     = 2,
    majorAlleleCopyNumber = c(flank_cn - 1, 49, flank_cn - 1),
    minorAlleleCopyNumber = c(1, 1, 1)
  )
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

  cancer_genes_gr <- GenomicRanges::GRanges(
    "chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR"
  )

  list(ecdna_gr = ecdna_gr, breakpoints_gr = breakpoints_gr,
       cnv_gr = cnv_gr, cancer_genes_gr = cancer_genes_gr)
}

make_plot_inputs <- function() {
  karyotype <- data.frame(
    chrom    = "chr7",
    start    = c(0, 50e6, 60e6),
    end      = c(50e6, 60e6, 159e6),
    gieStain = c("gneg", "gpos50", "gneg")
  )
  gene_coord <- data.frame(
    chr = "chr7", start = 55019017, end = 55211628, strand = "+", gene = "EGFR"
  )
  wgd_data <- data.frame(sample = "S1", Polyploidy = "No")
  cnv_data <- data.frame(
    sample   = "S1", seqnames = "chr7",
    start    = c(40e6, 55e6, 55.5e6 + 1),
    end      = c(54.99e6, 55.5e6, 70e6),
    copyNumber = c(2, 8, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 7, 1),
    minorAlleleCopyNumber = c(1, 1, 1)
  )
  sv_data <- data.frame(
    chrom1 = "7", start1 = 55e6, chrom2 = "7", start2 = 55.4e6,
    strand1 = "-", strand2 = "+", svclass = "DUP",
    VF = 1000, JCN = 40, sample = "S1"
  )
  list(karyotype = karyotype, gene_coord = gene_coord, wgd_data = wgd_data,
       cnv_data = cnv_data, sv_data = sv_data)
}
