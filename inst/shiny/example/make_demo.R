# Generate a small SYNTHETIC demo dataset (no patient data) for the EpiTracer app:
# an EGFR episomal ecDNA on chr7 and a BFB-like terminal staircase on chr8.
suppressPackageStartupMessages(library(data.table))

# NB: copy-number segments must be NON-OVERLAPPING (each start = previous end + 1).
cn <- rbindlist(list(
  data.table(seqnames = "chr7", start = c(1, 54000000, 56000001),
             end = c(53999999, 56000000, 159000000),
             copyNumber = c(2, 42, 2), minorAlleleCopyNumber = c(1, 1, 1),
             majorAlleleCopyNumber = c(1, 41, 1), ploidy = 2),
  data.table(seqnames = "chr8",
             start = c(1, 120000001, 128000001, 135000001, 140000001),
             end   = c(120000000, 128000000, 135000000, 140000000, 145138636),
             copyNumber = c(2, 6, 10, 16, 1), minorAlleleCopyNumber = c(1, 1, 1, 1, 0),
             majorAlleleCopyNumber = c(1, 5, 9, 15, 1), ploidy = 2)
))
cn[, sample := "DEMO1"]

sv <- data.table(
  chrom1 = c("7","7","8","8","8"), start1 = c(54000000, 54200000, 128000000, 135000000, 140000000),
  chrom2 = c("7","7","8","8","8"), start2 = c(56000000, 55800000, 128050000, 135050000, 140050000),
  strand1 = c("-","+","+","+","+"), strand2 = c("+","-","+","+","+"),
  svclass = c("DUP","DEL","h2hINV","h2hINV","h2hINV"),
  name = paste0("sv", 1:5), VF = c(1800, 300, 900, 700, 500),
  JCN = c(40, 2, 20, 16, 12), homlen = c(0, 0, 0, 0, 0), sample = "DEMO1")

fwrite(cn, file.path("inst/shiny/example/DEMO1_CN_segments.tsv"), sep = "\t")
fwrite(sv, file.path("inst/shiny/example/DEMO1_SV_bedpe.tsv"), sep = "\t")
cat("wrote demo:", nrow(cn), "CN segments,", nrow(sv), "SVs\n")
