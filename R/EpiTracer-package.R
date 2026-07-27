#' @keywords internal
"_PACKAGE"

## Quiet R CMD check NOTEs about non-standard evaluation used by
## data.table and dplyr (bare column names referenced in the code below).
utils::globalVariables(c(
  ".", "ID", "WGS_ID", "seqnames", "start", "end", "width", "strand",
  "event", "svclass", "VF", "PURPLE_CN", "PURPLE_AF", "PURPLE_JCN",
  "insLen", "HOMLEN", "homLen", "AF", "JCN",
  "sample", "copyNumber", "ploidy",
  "majorAlleleCopyNumber", "minorAlleleCopyNumber",
  "duplication_at_boundary", "duplication_at_boundary_has_highest_VF",
  "episome_region", "deletion_flanking_boundary",
  "before_prox_boundary_not_gained", "after_dist_boundary_not_gained",
  "episomal", "has_excision_scar",
  "chr", "chrom1", "chrom2", "start1", "start2", "strand1", "strand2",
  "strands", "gene", "gene_name", "pos_label", "min_coord", "max_coord",
  "Polyploidy", "y", "label", "pos", "gieStain", "color",
  "x", "xend", "yend", "lwd", "colour", "grp"
))
