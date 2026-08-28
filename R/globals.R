## Quiet R CMD check's "no visible binding for global variable" notes for the
## non-standard-evaluation column names used inside the ggplot builders and the
## data.table expressions of the simulation suite.
utils::globalVariables(c(
  "lab", "timing",
  ## simulation: junction ledger
  "svclass", "multiplicity", "mechanism", "origin", "homology",
  "chrom1", "start1", "strand1", "chrom2", "start2", "strand2",
  "jcn", "vf", "af",
  ## simulation: copy-number and breakpoint emission
  "seqnames", "start", "end", "grp", "sample", "copyNumber", "ploidy",
  "majorAlleleCopyNumber", "minorAlleleCopyNumber",
  "event", "bp_strand", "VF", "PURPLE_JCN", "WGS_ID",
  ## simulation: cohort and benchmark
  "ID", "class", "hit", "i.hit"
))
