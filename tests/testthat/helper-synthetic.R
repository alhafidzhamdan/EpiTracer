## Synthetic-data builders shared across tests.
## `flank_cn` = copy number of the segments IMMEDIATELY flanking the amplicon:
## low (2) => non-gained flanks (episome); high (10) => flanks gained above the
## chromosome baseline (not an episome). `chr_baseline` = copy number of the rest
## of the chromosome (the local baseline): 2 = diploid; higher models a polysomic
## chromosome (chr7 gain in glioblastoma). flank_cn == chr_baseline models a focal
## episome on a polysomic chromosome (missed by the "ploidy" flank test, recovered
## by "chromosome").

make_episome_inputs <- function(flank_cn = 2, chr_baseline = 2) {
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

  ## chromosome background at `chr_baseline`, immediate flanks at `flank_cn`,
  ## amplicon at 50. The background establishes the local per-chromosome baseline.
  seg_cn <- c(chr_baseline, flank_cn, 50, flank_cn, chr_baseline)
  cnv <- data.table::data.table(
    seqnames = "chr7",
    start    = c(1,        40000000, 55000000, 55500001, 70000001),
    end      = c(39999999, 54999999, 55500000, 70000000, 159000000),
    sample   = "S1",
    copyNumber = seg_cn,
    ploidy     = 2,
    majorAlleleCopyNumber = pmax(1, seg_cn - 1),
    minorAlleleCopyNumber = 1
  )
  cnv_gr <- GenomicRanges::makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

  cancer_genes_gr <- GenomicRanges::GRanges(
    "chr7", IRanges::IRanges(55019017, 55211628), gene = "EGFR"
  )

  list(ecdna_gr = ecdna_gr, breakpoints_gr = breakpoints_gr,
       cnv_gr = cnv_gr, cancer_genes_gr = cancer_genes_gr)
}

## Two amplified loci on different chromosomes joined by an inter-chromosomal
## TRA, for the concatenated-axis multi-locus plot.
make_multilocus_inputs <- function() {
  karyotype <- data.frame(
    chrom    = c("chr7", "chr7", "chr12", "chr12"),
    start    = c(0, 50e6, 0, 50e6),
    end      = c(50e6, 159e6, 50e6, 133e6),
    gieStain = c("gneg", "gpos50", "gneg", "gpos50")
  )
  gene_coord <- data.frame(
    chr = c("chr7", "chr12"), start = c(55019017, 57747727),
    end = c(55211628, 57756013), strand = c("+", "+"), gene = c("EGFR", "CDK4")
  )
  wgd_data <- data.frame(sample = "S1", Polyploidy = "No")
  cnv_data <- data.frame(
    sample   = "S1",
    seqnames = c("chr7", "chr7", "chr7", "chr12", "chr12", "chr12"),
    start    = c(40e6, 55e6, 56e6 + 1, 40e6, 57e6, 58e6 + 1),
    end      = c(54.99e6, 56e6, 70e6, 56.99e6, 58e6, 70e6),
    copyNumber = c(2, 20, 2, 2, 15, 2), ploidy = 2,
    majorAlleleCopyNumber = c(1, 19, 1, 1, 14, 1),
    minorAlleleCopyNumber = c(1, 1, 1, 1, 1, 1)
  )
  sv_data <- data.frame(
    chrom1 = c("7", "12"), start1 = c(55.2e6, 57.2e6),
    chrom2 = c("12", "7"), start2 = c(57.8e6, 55.8e6),
    strand1 = c("+", "-"), strand2 = c("-", "+"),
    svclass = c("TRA", "DUP"), VF = c(100, 60), JCN = c(5, 3), sample = "S1"
  )
  list(karyotype = karyotype, gene_coord = gene_coord, wgd_data = wgd_data,
       cnv_data = cnv_data, sv_data = sv_data)
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

## Synthetic small mutations for the optional SNV panel of plot_sv_linear().
## Uses the SNV column convention (`sampleID` sample, `allelic_freq` VAF, `type`
## mutation class). Five S1 SNVs inside the chr7 window (one with an artefactual
## VAF > 1), one S1 indel in-window (must be excluded), and an off-target SNV for
## a different sample (must be filtered out). With SNV-only + first-per-chromosome
## (intermutation-distance) handling, both `snv_y` modes yield 4 plotted points.
make_snv_inputs <- function() {
  data.frame(
    sampleID     = c(rep("S1", 6), "OTHER"),
    seqnames     = "chr7",
    start        = c(54.2e6, 55.10e6, 55.20e6, 55.25e6, 55.30e6, 55.40e6, 55.15e6),
    allelic_freq = c(0.30,   0.95,    0.42,    5.0,      0.18,    0.60,    0.5),  # 5.0 = artefact
    variant_cn   = c(2,      1,       3,       1,        8,       1,       1),
    type         = c("SNV",  "SNV",   "SNV",   "SNV",    "SNV",   "Indel", "SNV"),
    ref          = "A", mut = "G",
    stringsAsFactors = FALSE
  )
}
