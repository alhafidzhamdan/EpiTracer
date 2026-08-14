#' Classify a single ecDNA amplicon as episomal or not
#'
#' Internal worker called once per amplicon by [call_episomal_ecdna()]. It
#' annotates the structural variant breakpoints falling within (or just
#' outside) an ecDNA amplicon, then applies the episome heuristic:
#' \enumerate{
#'   \item find duplication (DUP) breakpoints that sit at the amplicon
#'     boundaries and are themselves amplified;
#'   \item require the boundary DUP to carry the highest variant fraction (VF)
#'     of any DUP in the amplicon;
#'   \item require the chromosomal segments immediately flanking the boundaries
#'     to be non-gained (consistent with a circular episome excised from an
#'     otherwise diploid region);
#'   \item flag a shared flanking deletion as a candidate excision scar.
#' }
#'
#' @param this_amplicon_id Character scalar; a single value of `ecdna_gr$ID`.
#' @param ecdna_gr,breakpoints_gr,cnv_gr,cancer_genes_gr See
#'   [call_episomal_ecdna()].
#' @param ext Integer; bp to extend the amplicon by when searching for
#'   boundary SVs.
#' @param verbose Logical; print progress/diagnostics.
#'
#' @return A [data.table::data.table] of annotated breakpoints for the
#'   amplicon (may have zero rows), or an empty data.table if the amplicon has
#'   no breakpoints in range.
#' @keywords internal
#' @importFrom data.table data.table rbindlist :=
#' @importFrom GenomicRanges trim
#' @importFrom regioneR toGRanges
#' @importFrom dplyr filter select arrange rename mutate
classify_amplicon_episomal <- function(this_amplicon_id,
                                       ecdna_gr,
                                       breakpoints_gr,
                                       cnv_gr,
                                       cancer_genes_gr,
                                       ext = 1e7,
                                       verbose = FALSE) {

  this_amplicon_gr <- ecdna_gr[ecdna_gr$ID %in% this_amplicon_id]

  ## Drop annotation columns that are re-derived below (if present):
  this_amplicon_gr$gene  <- NULL
  this_amplicon_gr$group <- NULL
  this_amplicon_gr$arm   <- NULL

  this_sample <- unique(this_amplicon_gr$WGS_ID)
  ## Base GRanges subsetting (avoids requiring plyranges for filter.GRanges):
  this_sample_cnv_gr <- cnv_gr[cnv_gr$sample == this_sample]

  this_sample_breakpoints <- breakpoints_gr[breakpoints_gr$WGS_ID == this_sample] %>%
    gr2dt() %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::select(seqnames, start, end, event, svclass,
                  AF = PURPLE_AF, JCN = PURPLE_JCN, VF, PURPLE_CN,
                  insLen, homLen = HOMLEN) %>%
    regioneR::toGRanges(genome = "hg38")

  ## Breakpoints at (or just outside) the amplicon boundaries:
  this_sample_breakpoints_ecdna_annotated <-
    (this_sample_breakpoints %$% GenomicRanges::trim((this_amplicon_gr + ext))) %>%
    gr2dt() %>%
    dplyr::filter(seqnames %in% as.character(this_amplicon_gr@seqnames)) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(ID != "") %>%
    dplyr::select(-c(strand, width))

  ## Add oncogene info:
  this_sample_breakpoints_ecdna_annotated <-
    ((this_sample_breakpoints_ecdna_annotated %>% regioneR::toGRanges()) %$% cancer_genes_gr) %>%
    gr2dt() %>%
    dplyr::select(-c(strand, width))

  ## Annotate with minor and major allele info:
  this_sample_breakpoints_ecdna_annotated <-
    ((this_sample_breakpoints_ecdna_annotated %>% regioneR::toGRanges()) %$% this_sample_cnv_gr) %>%
    gr2dt() %>%
    dplyr::select(-c(strand, width))

  if (nrow(this_sample_breakpoints_ecdna_annotated) > 0) {

    this_sample_breakpoints_ecdna_annotated$PURPLE_CN <-
      as.numeric(this_sample_breakpoints_ecdna_annotated$PURPLE_CN)
    this_sample_ploidy <- as.numeric(this_sample_breakpoints_ecdna_annotated$ploidy[1])

    ## Prep columns:
    this_sample_breakpoints_ecdna_annotated$duplication_at_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$duplication_at_boundary_has_highest_VF <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$episome_region <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$deletion_flanking_boundary <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$before_prox_boundary_not_gained <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$after_dist_boundary_not_gained <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$episomal <- "FALSE"
    this_sample_breakpoints_ecdna_annotated$has_excision_scar <- "FALSE"

    ## The same amplicon can span >1 chromosome, so classify per chromosome:
    unique_chrs <- unique(this_sample_breakpoints_ecdna_annotated$seqnames)

    this_sample_breakpoints_ecdna_annotated <- lapply(seq_along(unique_chrs), function(i) {

      this_chr <- this_sample_breakpoints_ecdna_annotated %>%
        dplyr::filter(seqnames %in% unique_chrs[i])

      ## Min/max coords of amplified regions for this chromosome:
      has_amp_region <- nrow(this_sample_cnv_gr %>% gr2dt() %>%
                               dplyr::filter(seqnames %in% unique_chrs[i]) %>%
                               dplyr::filter(copyNumber > 3 * ploidy)) > 0

      if (has_amp_region) {

        min_amp_coord <- this_sample_cnv_gr %>% gr2dt() %>%
          dplyr::filter(seqnames %in% unique_chrs[i]) %>%
          dplyr::filter(copyNumber > 3 * ploidy) %>%
          dplyr::filter(start >= min(this_amplicon_gr[as.character(GenomicRanges::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@start)) %>%
          dplyr::filter(start == min(start)) %>% .$start

        max_amp_coord <- this_sample_cnv_gr %>% gr2dt() %>%
          dplyr::filter(seqnames %in% unique_chrs[i]) %>%
          dplyr::filter(copyNumber > 3 * ploidy) %>%
          dplyr::filter(end <= max(this_amplicon_gr[as.character(GenomicRanges::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@start +
                                     this_amplicon_gr[as.character(GenomicRanges::seqnames(this_amplicon_gr)) %in% unique_chrs[i]]@ranges@width)) %>%
          dplyr::filter(end == max(end)) %>% .$end

        ## Find DUP at boundary; ensure DUP has highest VF; check flanking CN;
        ## check for excision scar.

        ## needs to be >1 DUP: a single DUP on a separate chr is unlikely to
        ## be involved in amplicon generation (outside of ecDNA)
        has_dup <- nrow(this_chr %>%
                          dplyr::filter(svclass == "DUP") %>%
                          dplyr::filter(PURPLE_CN > 3 * this_sample_ploidy)) > 1

        if (has_dup) {

          ## Proximal border:
          has_prox_border <- nrow(this_chr %>%
                                    dplyr::filter(svclass == "DUP") %>%
                                    dplyr::filter(PURPLE_CN > 3 * this_sample_ploidy) %>%
                                    dplyr::filter(start < min_amp_coord + 10000)) > 0
          if (has_prox_border) {
            prox_border <- this_chr %>%
              dplyr::filter(svclass == "DUP") %>%
              dplyr::filter(PURPLE_CN > 3 * this_sample_ploidy) %>%
              dplyr::filter(start < min_amp_coord + 10000) %>%
              dplyr::filter(start == min(start)) %>% .$event
          }

          ## Distal border:
          has_dist_border <- nrow(this_chr %>%
                                    dplyr::filter(svclass == "DUP") %>%
                                    dplyr::filter(PURPLE_CN > 3 * this_sample_ploidy) %>%
                                    dplyr::filter(end > max_amp_coord - 10000)) > 0
          if (has_dist_border) {
            dist_border <- this_chr %>%
              dplyr::filter(svclass == "DUP") %>%
              dplyr::filter(PURPLE_CN > 3 * this_sample_ploidy) %>%
              dplyr::filter(end > max_amp_coord - 10000) %>%
              dplyr::filter(start == max(start)) %>% .$event
          }

          if (has_prox_border & has_dist_border) {
            if (length(prox_border) == 1 & length(dist_border) == 1) {
              if (prox_border == dist_border) {
                this_chr[event %in% dist_border]$duplication_at_boundary <- "TRUE"
              }
            }
          }

          ## Ensure boundary DUP has highest VF:
          boundary_index <- which(this_chr$duplication_at_boundary == "TRUE")

          if (length(boundary_index) > 0) {

            ## max VF for DUP only (allow other internal SVs with higher VF)
            max_vf <- max(this_chr %>% dplyr::filter(svclass == "DUP") %>% .$VF)

            if (unique(this_chr[boundary_index]$VF) == max_vf) {

              this_chr[boundary_index]$duplication_at_boundary_has_highest_VF <- "TRUE"

              ## Classify episome:
              prox_boundary <- this_chr[boundary_index[1]]$start
              dist_boundary <- this_chr[boundary_index[2]]$start

              ## Chromosomal segments < prox_boundary and > dist_boundary
              ## should not be gained/amplified:
              before_prox_boundary_not_gained <- nrow(this_sample_cnv_gr %>% gr2dt() %>%
                                                         dplyr::filter(seqnames %in% unique_chrs[i]) %>%
                                                         dplyr::filter(end < prox_boundary) %>%
                                                         dplyr::filter(end == max(end)) %>%
                                                         dplyr::filter(copyNumber < 1.4 * ploidy)) > 0

              after_dist_boundary_not_gained <- nrow(this_sample_cnv_gr %>% gr2dt() %>%
                                                        dplyr::filter(seqnames %in% unique_chrs[i]) %>%
                                                        dplyr::filter(start > dist_boundary) %>%
                                                        dplyr::filter(start == min(start)) %>%
                                                        dplyr::filter(copyNumber < 1.4 * ploidy)) > 0

              ## NOTE: the original script tested
              ##   `after_dist_boundary_not_gained == after_dist_boundary_not_gained`
              ## which is a tautology (always TRUE). Corrected here to require
              ## BOTH flanks to be non-gained before calling an episome region.
              if (before_prox_boundary_not_gained & after_dist_boundary_not_gained) {
                this_chr$before_prox_boundary_not_gained <- "TRUE"
                this_chr$after_dist_boundary_not_gained <- "TRUE"
                this_chr[start >= prox_boundary & start <= dist_boundary,
                         episome_region := "TRUE"]
              }

              ## Define excision scar:
              prox_boundary_sv <- this_chr[boundary_index[1] - 1]$event
              dist_boundary_sv <- this_chr[boundary_index[2] + 1]$event

              if (length(prox_boundary_sv) > 0 & length(dist_boundary_sv) > 0) {
                if (!is.na(prox_boundary_sv) & !is.na(dist_boundary_sv)) {
                  if (prox_boundary_sv == dist_boundary_sv) {
                    if (unique(this_chr[event == prox_boundary_sv]$svclass == "DEL")) {
                      this_chr[event == prox_boundary_sv]$deletion_flanking_boundary <- "TRUE"
                    }
                  }
                }
              }
            }
          }
        }
      }

      ## Classify the whole amplicon (per chromosome):
      if (nrow(this_chr %>% dplyr::filter(episome_region == "TRUE")) > 0) {
        this_chr$episomal <- "TRUE"
      }
      if (nrow(this_chr %>% dplyr::filter(deletion_flanking_boundary == "TRUE")) == 2) {
        this_chr$has_excision_scar <- "TRUE"
      }

      this_chr
    }) %>% data.table::rbindlist()
  }

  this_sample_breakpoints_ecdna_annotated
}


#' Call episomal extrachromosomal DNA from WGS structural variants
#'
#' Detects ecDNA amplicons whose structure is consistent with the *episome*
#' (breakage-independent) model of formation: a circular amplicon bounded by a
#' duplication breakpoint, arising from an otherwise non-amplified chromosomal
#' region, and often leaving a deletion "excision scar" at the origin locus.
#'
#' Each amplicon called by AmpliconArchitect (supplied in `ecdna_gr`, a
#' compulsory input) is processed independently (optionally in parallel). Within
#' each amplicon the function annotates the flanking structural variant
#' breakpoints with oncogene and allele-specific copy-number context, then
#' applies the heuristic described in [classify_amplicon_episomal()].
#'
#' @param ecdna_gr A [GenomicRanges::GRanges] of ecDNA amplicon regions. Must
#'   contain metadata columns `ID` (unique amplicon identifier) and `WGS_ID`
#'   (sample identifier). Typically the AmpliconArchitect ecDNA amplicon
#'   catalogue.
#' @param breakpoints_gr A [GenomicRanges::GRanges] of PURPLE (HMF pipeline) SV
#'   breakpoints with metadata columns `WGS_ID`, `event`, `svclass` (e.g. "DUP",
#'   "DEL"), `PURPLE_AF`, `PURPLE_JCN`, `VF`, `PURPLE_CN`, `insLen`, `HOMLEN`.
#' @param cnv_gr A [GenomicRanges::GRanges] of PURPLE (HMF pipeline)
#'   allele-specific copy-number segments with metadata columns `sample`,
#'   `copyNumber`, `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`.
#' @param cancer_genes_gr A [GenomicRanges::GRanges] of cancer gene loci with a
#'   `gene` (gene symbol) metadata column, used to annotate breakpoints with the
#'   overlapping oncogene.
#' @param ext Integer; number of base pairs to extend each amplicon by when
#'   searching for boundary breakpoints (default `1e7`).
#' @param mc.cores Integer; number of cores for [parallel::mclapply()]
#'   (default `1`). Values > 1 are ignored on Windows.
#' @param verbose Logical; print per-amplicon progress (default `FALSE`).
#'
#' @return A [data.table::data.table] combining the annotated breakpoints of all
#'   amplicons, with per-breakpoint classification columns including
#'   `duplication_at_boundary`, `duplication_at_boundary_has_highest_VF`,
#'   `episome_region`, `deletion_flanking_boundary`, `episomal`, and
#'   `has_excision_scar` (all character `"TRUE"`/`"FALSE"`).
#'
#' @examples
#' \dontrun{
#' ecdna_gr        <- readRDS("ecDNA_amplicon_regions.rds")
#' breakpoints_gr  <- readRDS("SV_catalogue.rds")
#' cnv_gr          <- readRDS("CN_segments.rds")
#' cancer_genes_gr <- readRDS("cancer_genes.rds")
#'
#' episomal <- call_episomal_ecdna(
#'   ecdna_gr, breakpoints_gr, cnv_gr, cancer_genes_gr,
#'   ext = 1e7, mc.cores = 4
#' )
#' }
#'
#' @seealso [classify_amplicon_episomal()], [plot_sv_linear()]
#' @export
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
call_episomal_ecdna <- function(ecdna_gr,
                                breakpoints_gr,
                                cnv_gr,
                                cancer_genes_gr,
                                ext = 1e7,
                                mc.cores = 1L,
                                verbose = FALSE) {

  stopifnot(
    methods::is(ecdna_gr, "GRanges"),
    methods::is(breakpoints_gr, "GRanges"),
    methods::is(cnv_gr, "GRanges"),
    methods::is(cancer_genes_gr, "GRanges"),
    !is.null(ecdna_gr$ID), !is.null(ecdna_gr$WGS_ID)
  )

  ## Caller-agnostic column checks. EpiTracer works with SV/CN calls from any
  ## source (not only PURPLE / AmpliconArchitect) once coerced to these columns;
  ## fail early with an explicit message naming what is missing.
  .need <- function(gr, cols, what) {
    miss <- setdiff(cols, names(S4Vectors::mcols(gr)))
    if (length(miss))
      stop(sprintf("%s is missing required metadata column(s): %s",
                   what, paste(miss, collapse = ", ")), call. = FALSE)
  }
  .need(breakpoints_gr,
        c("WGS_ID", "event", "svclass", "PURPLE_AF", "PURPLE_JCN",
          "VF", "PURPLE_CN", "insLen", "HOMLEN"), "breakpoints_gr")
  .need(cnv_gr,
        c("sample", "copyNumber", "ploidy",
          "majorAlleleCopyNumber", "minorAlleleCopyNumber"), "cnv_gr")
  .need(cancer_genes_gr, "gene", "cancer_genes_gr")

  amplicon_ids <- unique(ecdna_gr$ID)
  if (verbose) message("Classifying ", length(amplicon_ids), " amplicons ...")

  results <- parallel::mclapply(seq_along(amplicon_ids), function(x) {
    if (verbose) message("[", x, "/", length(amplicon_ids), "] ", amplicon_ids[x])
    classify_amplicon_episomal(
      this_amplicon_id = amplicon_ids[x],
      ecdna_gr         = ecdna_gr,
      breakpoints_gr   = breakpoints_gr,
      cnv_gr           = cnv_gr,
      cancer_genes_gr  = cancer_genes_gr,
      ext              = ext,
      verbose          = verbose
    )
  }, mc.cores = mc.cores)

  data.table::rbindlist(results, fill = TRUE)
}
