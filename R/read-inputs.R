# ============================================================================
# Flexible input readers: turn generic / PURPLE / GRIDSS copy-number and
# structural-variant files into the tables EpiTracer's plotters and callers use.
# ============================================================================

#' @importFrom data.table := .I .N as.data.table fread setnames rbindlist fifelse data.table
#' @importFrom stats median
NULL

`%||%` <- function(a, b) if (is.null(a) || !length(a) || (is.character(a) && !nzchar(a[1]))) b else a

.CN_ALIASES <- list(
  sample = c("sample","sampleid","sample_id","wgs_id","tumor","tumour","id"),
  seqnames = c("seqnames","chrom","chromosome","chr","contig"),
  start = c("start","startpos","start_position","begin"),
  end = c("end","endpos","end_position","stop"),
  copyNumber = c("copynumber","cn","totalcn","total_cn","tcn","copy_number"),
  ploidy = c("ploidy"),
  majorAlleleCopyNumber = c("majorallelecopynumber","major","majorcn","major_cn","nmajor","cnmajor"),
  minorAlleleCopyNumber = c("minorallelecopynumber","minor","minorcn","minor_cn","nminor","cnminor")
)
.SV_ALIASES <- list(
  sample = c("sample","sampleid","sample_id","wgs_id","tumor","tumour"),
  chrom1 = c("chrom1","chr1","chromosome1","chra"),
  start1 = c("start1","pos1","position1","starta"),
  chrom2 = c("chrom2","chr2","chromosome2","chrb"),
  start2 = c("start2","pos2","position2","startb"),
  strand1 = c("strand1","orientation1","stranda"),
  strand2 = c("strand2","orientation2","strandb"),
  svclass = c("svclass","svtype","type","class"),
  VF = c("vf","support","read_support","tumreads","pr"),
  JCN = c("jcn","cn_change","ploidy_jcn","junction_copy_number"),
  name = c("name","event","id","sv_id"),
  homlen = c("homlen","homology","homlength","microhomology")
)

.canon <- function(dt, aliases) {
  low <- tolower(names(dt))
  for (target in names(aliases)) {
    if (target %in% names(dt)) next
    hit <- which(low %in% aliases[[target]])
    if (length(hit)) data.table::setnames(dt, names(dt)[hit[1]], target)
  }
  dt
}
.pfx <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))

## sample id from a PURPLE/generic file name, stripping the usual suffixes
.strip_sample_name <- function(name) {
  b <- basename(as.character(name %||% ""))
  b <- sub("\\.gz$", "", b, ignore.case = TRUE)
  b <- sub("\\.(rds|tsv|csv|txt|bed|bedpe|vcf)$", "", b, ignore.case = TRUE)
  b <- sub("\\.purple\\..*$", "", b, ignore.case = TRUE)
  b <- sub("\\.(svs?|cnv|copynumber|copy_number|segments?|somatic)([._-].*)?$", "", b, ignore.case = TRUE)
  if (nzchar(b)) b else "sample1"
}

## length-weighted mean copy number over autosomes -- a PURPLE-like ploidy estimate
## (the cnv.somatic.tsv carries no ploidy; it lives in .purple.purity.tsv).
.estimate_ploidy <- function(dt) {
  auto <- dt[grepl("^chr[0-9]+$", seqnames) & is.finite(copyNumber) & copyNumber < 15]
  if (!nrow(auto)) auto <- dt[is.finite(copyNumber)]
  w <- pmax(1, auto$end - auto$start)
  round(max(1.5, sum(auto$copyNumber * w) / sum(w)), 2)
}

.norm_svclass <- function(dt) {
  sc <- toupper(as.character(dt$svclass)); inter <- as.character(dt$chrom1) != as.character(dt$chrom2)
  s1 <- as.character(dt$strand1); s2 <- as.character(dt$strand2)
  out <- sc; out[inter] <- "TRA"
  keep <- sc %in% c("DEL","DUP","INS","TRA","H2HINV","T2TINV")
  out[keep & !inter] <- sc[keep & !inter]
  out[out == "H2HINV"] <- "h2hINV"; out[out == "T2TINV"] <- "t2tINV"
  amb <- !inter & !(sc %in% c("DEL","DUP","INS"))
  out[amb & s1 == "+" & s2 == "+"] <- "h2hINV"; out[amb & s1 == "-" & s2 == "-"] <- "t2tINV"
  out[amb & s1 == "+" & s2 == "-"] <- "DEL";     out[amb & s1 == "-" & s2 == "+"] <- "DUP"
  out
}

.read_table <- function(file, name) {
  if (grepl("\\.rds$", name, ignore.case = TRUE)) {
    obj <- readRDS(file); if (methods::is(obj, "GRanges")) obj <- as.data.frame(obj)
    data.table::as.data.table(obj)
  } else data.table::fread(file)
}

.as_dt <- function(x, name) {
  if (is.character(x) && length(x) == 1 && file.exists(x)) list(dt = .read_table(x, name %||% x), name = name %||% basename(x))
  else list(dt = data.table::as.data.table(x), name = name %||% "sample")
}

#' Read and reformat copy-number segments for EpiTracer
#'
#' Turns a copy-number table from (almost) any source into the columns EpiTracer's
#' plotters and callers expect: `sample`, `seqnames`, `start`, `end`, `copyNumber`,
#' `ploidy`, `majorAlleleCopyNumber`, `minorAlleleCopyNumber`. Column names are
#' matched case-insensitively, so a raw PURPLE `*.purple.cnv.somatic.tsv` works
#' directly. PURPLE's occasional negative copy numbers are clamped to 0; the allele
#' copy numbers are derived if absent; and, because a PURPLE segment file carries no
#' ploidy (it lives in `*.purple.purity.tsv`), ploidy is estimated from the segments
#' unless supplied.
#'
#' @param x A file path (TSV/CSV, `*.purple.cnv.somatic.tsv`, or a `.rds` of a
#'   GRanges/data.frame) or an in-memory data.frame/data.table.
#' @param sample Optional sample name. When given it **overrides** any `sample`
#'   column and the file name, labelling every row with this name -- use it to give
#'   CN and SV files a common name when their file-naming differs.
#' @param ploidy Optional numeric ploidy. When given it overrides the estimate.
#' @param name Optional file name used for format detection and to derive the sample
#'   when `x` is a temporary path (e.g. a Shiny upload) whose own name is opaque.
#' @return A [data.table::data.table] in EpiTracer's copy-number format.
#' @examples
#' \dontrun{
#' cnv <- read_cnv("DO13869T.purple.cnv.somatic.tsv", sample = "DO13869T")
#' }
#' @export
read_cnv <- function(x, sample = NULL, ploidy = NULL, name = NULL) {
  r <- .as_dt(x, name); dt <- .canon(r$dt, .CN_ALIASES)
  miss <- setdiff(c("seqnames","start","end","copyNumber"), names(dt))
  if (length(miss)) stop("copy-number input is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  dt[, seqnames := .pfx(as.character(seqnames))]
  dt[, `:=`(start = as.numeric(start), end = as.numeric(end), copyNumber = pmax(0, as.numeric(copyNumber)))]
  if ("minorAlleleCopyNumber" %in% names(dt)) dt[, minorAlleleCopyNumber := pmax(0, as.numeric(minorAlleleCopyNumber))]
  if ("majorAlleleCopyNumber" %in% names(dt)) dt[, majorAlleleCopyNumber := pmax(0, as.numeric(majorAlleleCopyNumber))]
  if (!"minorAlleleCopyNumber" %in% names(dt)) dt[, minorAlleleCopyNumber := pmin(1, copyNumber/2)]
  if (!"majorAlleleCopyNumber" %in% names(dt)) dt[, majorAlleleCopyNumber := pmax(0, copyNumber - minorAlleleCopyNumber)]
  if (!is.null(sample)) { .snm <- as.character(sample)[1]; dt[, sample := .snm] }
  else if (!"sample" %in% names(dt)) { .snm <- .strip_sample_name(r$name); dt[, sample := .snm] }
  dt[, sample := as.character(sample)]
  pl <- if (!is.null(ploidy) && is.finite(ploidy) && ploidy > 0) as.numeric(ploidy)
        else if ("ploidy" %in% names(dt)) as.numeric(dt$ploidy) else .estimate_ploidy(dt)
  dt[, ploidy := pl]
  dt[]
}

#' Read and reformat structural variants for EpiTracer
#'
#' Accepts a BEDPE (with `VF`/`JCN`), a PURPLE/GRIDSS breakend VCF (parsed by
#' [read_purple_sv_vcf()]), or a `.rds`, and returns the columns EpiTracer uses:
#' `chrom1`, `start1`, `chrom2`, `start2`, `strand1`, `strand2`, `svclass`, `VF`,
#' `JCN`, `name`, `homlen`, `sample`. `svclass` is normalised to
#' `DEL`/`DUP`/`INS`/`TRA`/`h2hINV`/`t2tINV`. Keeping **VF** (supporting fragments)
#' and **JCN** (junction copy number) matters -- the mechanism callers use them.
#'
#' @param x A file path (BEDPE TSV/CSV, `*.vcf`/`*.vcf.gz`, or `.rds`) or a
#'   data.frame/data.table.
#' @param sample Optional sample name; overrides any `sample` column and the file
#'   name (see [read_cnv()]).
#' @param name Optional file name for format detection / sample derivation.
#' @return A [data.table::data.table] in EpiTracer's SV (BEDPE) format.
#' @export
read_sv <- function(x, sample = NULL, name = NULL) {
  nm <- if (is.character(x) && length(x) == 1) (name %||% basename(x)) else (name %||% "sample")
  dt <- if (is.character(x) && length(x) == 1 && grepl("\\.vcf(\\.gz)?$", nm, ignore.case = TRUE)) read_purple_sv_vcf(x)
        else .as_dt(x, nm)$dt
  dt <- .canon(dt, .SV_ALIASES)
  miss <- setdiff(c("chrom1","start1","chrom2","start2","svclass"), names(dt))
  if (length(miss)) stop("SV input is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  for (col in c("strand1","strand2")) if (!col %in% names(dt)) dt[, (col) := "+"]
  if (!"VF" %in% names(dt)) dt[, VF := 100]
  if (!"JCN" %in% names(dt)) dt[, JCN := 1]
  if (!"name" %in% names(dt)) dt[, name := paste0("SV", .I)]
  if (!"homlen" %in% names(dt)) dt[, homlen := 0]
  if (!is.null(sample)) { .snm <- as.character(sample)[1]; dt[, sample := .snm] }
  else if (!"sample" %in% names(dt)) { .snm <- .strip_sample_name(nm); dt[, sample := .snm] }
  dt[, `:=`(chrom1 = .pfx(as.character(chrom1)), chrom2 = .pfx(as.character(chrom2)),  # standardise: 'chr' prefix
            start1 = as.numeric(start1), start2 = as.numeric(start2),
            VF = as.numeric(VF), JCN = as.numeric(JCN), sample = as.character(sample))]
  dt[, svclass := .norm_svclass(dt)]
  dt[]
}

#' Convert a PURPLE / GRIDSS breakend VCF to BEDPE, keeping VF and JCN
#'
#' Parses breakend (`BND`) records from a PURPLE/GRIDSS `*.vcf`/`*.vcf.gz` into the
#' paired-breakend BEDPE EpiTracer uses, **preserving `VF`** (variant fragments,
#' from the `VF` INFO field or the tumour sample's `VF` FORMAT field) and **`JCN`**
#' (from `PURPLE_JCN`). A plain coordinate-only VCF-to-BEDPE conversion drops these,
#' which degrades the mechanism calls. Mate breakends are de-duplicated and strand /
#' orientation is read from the breakend ALT notation (`N[chr:pos[` etc.), so
#' `svclass` follows the usual convention (`+/-`=DEL, `-/+`=DUP, `+/+`=h2hINV,
#' `-/-`=t2tINV, inter-chromosomal=TRA).
#'
#' @param file Path to a `*.vcf` or `*.vcf.gz`.
#' @return A [data.table::data.table] with BEDPE columns plus `VF`, `JCN`, `homlen`.
#' @export
read_purple_sv_vcf <- function(file) {
  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file) else file(file)
  lines <- readLines(con); close(con)
  body <- lines[!startsWith(lines, "#")]
  if (!length(body)) return(data.table::data.table())
  f <- strsplit(body, "\t", fixed = TRUE)
  g <- function(i) vapply(f, function(x) if (length(x) >= i) x[i] else NA_character_, "")
  chrom <- g(1); pos <- as.numeric(g(2)); id <- g(3); alt <- g(5)
  filt <- g(7); info <- g(8); format <- g(9)
  keep <- grepl("[][]", alt) & (is.na(filt) | filt %in% c("PASS",".",""))
  if (!any(keep)) return(data.table::data.table())
  iv <- function(inf, key) { has <- grepl(paste0("(?:^|;)", key, "="), inf, perl = TRUE)
    v <- sub(paste0(".*(?:^|;)", key, "=([^;]*).*"), "\\1", inf, perl = TRUE); v[!has] <- NA; v }
  fv <- function(i, key) { if (is.na(format[i])) return(NA_character_)
    ks <- strsplit(format[i], ":", fixed = TRUE)[[1]]; k <- match(key, ks)
    if (is.na(k) || length(f[[i]]) < 10) return(NA_character_)
    vals <- strsplit(f[[i]][length(f[[i]])], ":", fixed = TRUE)[[1]]
    if (length(vals) >= k) vals[k] else NA_character_ }
  rows <- lapply(which(keep), function(i) {
    m <- regmatches(alt[i], regexec("([][])([^:]+):([0-9]+)([][])", alt[i]))[[1]]
    if (length(m) < 5) return(NULL)
    o1 <- if (grepl("^[ACGTNacgtn]", alt[i])) "+" else "-"
    o2 <- if (m[2] == "[") "-" else "+"
    c1 <- chrom[i]; p1 <- pos[i]; c2 <- m[3]; p2 <- as.numeric(m[4])
    if (c2 < c1 || (c2 == c1 && p2 < p1)) { t<-c1;c1<-c2;c2<-t; t<-p1;p1<-p2;p2<-t; t<-o1;o1<-o2;o2<-t }
    vf <- suppressWarnings(as.numeric(iv(info[i], "VF"))); if (is.na(vf)) vf <- suppressWarnings(as.numeric(fv(i, "VF")))
    jcn <- suppressWarnings(as.numeric(iv(info[i], "PURPLE_JCN"))); if (is.na(jcn)) jcn <- suppressWarnings(as.numeric(iv(info[i], "JCN")))
    hom <- suppressWarnings(as.numeric(iv(info[i], "HOMLEN")))
    data.table::data.table(chrom1 = .pfx(c1), start1 = p1, chrom2 = .pfx(c2), start2 = p2,  # 'chr' prefix
      strand1 = o1, strand2 = o2, name = id[i], svclass = NA_character_,
      homlen = ifelse(is.na(hom), 0, hom), VF = vf, JCN = jcn, jkey = paste(c1,p1,c2,p2))
  })
  dt <- data.table::rbindlist(Filter(Negate(is.null), rows))
  if (!nrow(dt)) return(dt)
  dt <- dt[!duplicated(jkey)][, jkey := NULL]
  inter <- dt$chrom1 != dt$chrom2
  dt[, svclass := data.table::fifelse(inter, "TRA",
                  data.table::fifelse(strand1 == "+" & strand2 == "-", "DEL",
                  data.table::fifelse(strand1 == "-" & strand2 == "+", "DUP",
                  data.table::fifelse(strand1 == "+" & strand2 == "+", "h2hINV", "t2tINV"))))]
  dt[]
}

#' Resolve gene symbol(s) to a plotting locus
#'
#' Looks a gene symbol up in the bundled oncogene panel for `genome` and returns a
#' window suitable for [plot_sv_linear()]'s `loci` argument, so a user can ask to
#' plot e.g. `"EGFR"` or `"CDKN2A"` by name. When `cnv` is supplied the window is
#' grown to the copy-number event (amplification or homozygous deletion) overlapping
#' the gene, so the zoom fits the event; otherwise the gene span plus `flank` is used.
#'
#' @param genes Character vector of gene symbols (case-insensitive).
#' @param genome One of `"hg38"`, `"hg19"`, `"mm10"`.
#' @param cnv Optional copy-number table (from [read_cnv()]) to size the window to
#'   the overlapping event.
#' @param sample Optional sample to subset `cnv` by.
#' @param flank Bp added on each side of the resolved region (default `5e5`).
#' @param target `"amp"` sizes to an amplified segment, `"homdel"` to a homozygous
#'   deletion, `"any"` (default) to the gene span.
#' @param min_cn_ratio,homdel_thresh Amplification (`copyNumber > min_cn_ratio *
#'   ploidy`) and homozygous-deletion (`copyNumber < homdel_thresh`) thresholds.
#' @return A data.frame with `chr`, `start`, `end`, `gene` (one row per found gene).
#' @seealso [plot_sv_linear()], [read_cnv()]
#' @export
gene_locus <- function(genes, genome = "hg38", cnv = NULL, sample = NULL, flank = 5e5,
                       target = c("any","amp","homdel"), min_cn_ratio = 3, homdel_thresh = 0.5) {
  target <- match.arg(target)
  bed <- system.file("extdata", paste0("oncogene_coord_", genome, ".bed"), package = "EpiTracer")
  if (!nzchar(bed)) stop("no bundled gene coordinates for genome '", genome, "'", call. = FALSE)
  panel <- utils::read.table(bed, sep = "\t", col.names = c("chr","start","end","strand","gene"))
  want <- toupper(trimws(genes))
  hit <- panel[toupper(panel$gene) %in% want, , drop = FALSE]
  miss <- setdiff(want, toupper(hit$gene))
  if (length(miss)) warning("gene(s) not in the panel: ", paste(miss, collapse = ", "), call. = FALSE)
  if (!nrow(hit)) stop("none of the requested genes are in the bundled panel", call. = FALSE)
  cdt <- NULL
  if (!is.null(cnv)) { cdt <- data.table::as.data.table(cnv)
    if (!is.null(sample) && "sample" %in% names(cdt)) { sid <- sample; cdt <- cdt[sample == sid] } }
  do.call(rbind, lapply(seq_len(nrow(hit)), function(i) {
    g <- hit[i, ]; chr <- .pfx(as.character(g$chr)); s <- g$start; e <- g$end
    if (!is.null(cdt) && nrow(cdt)) {
      seg <- cdt[.pfx(as.character(seqnames)) == chr & end >= s & start <= e]
      if (target == "amp")    seg <- seg[copyNumber > min_cn_ratio * ploidy]
      if (target == "homdel") seg <- seg[copyNumber < homdel_thresh]
      if (nrow(seg)) { s <- min(s, min(seg$start)); e <- max(e, max(seg$end)) }
    }
    data.frame(chr = chr, start = max(1, s - flank), end = e + flank, gene = g$gene, stringsAsFactors = FALSE)
  }))
}

#' Build EpiTracer caller inputs from copy-number and SV tables
#'
#' Assembles the `cnv_gr` and `breakpoints_gr` [GenomicRanges::GRanges] that the
#' amplicon-mechanism callers ([call_simple_excision()], [call_brf()], etc.) take,
#' from the tables returned by [read_cnv()] / [read_sv()]. If the CN and SV tables
#' use different sample labels (common with PURPLE file naming), pass `sample` to
#' force a single name; the SVs are then taken from the whole SV table.
#'
#' @param cnv,sv Copy-number and SV tables (from [read_cnv()] / [read_sv()]).
#' @param sample Sample to build for. Defaults to the first CN sample.
#' @param ploidy_mode `"per_chr"` (default) sets each breakend/segment ploidy to the
#'   local per-chromosome baseline (recovers amplicons on polysomic chromosomes);
#'   `"global"` uses the CN table's own ploidy.
#' @return A list with `cnv_gr` and `breakpoints_gr`.
#' @seealso [read_cnv()], [read_sv()], [call_simple_excision()]
#' @export
prepare_amplicon_inputs <- function(cnv, sv, sample = NULL, ploidy_mode = c("per_chr","global")) {
  ploidy_mode <- match.arg(ploidy_mode)
  cnv <- data.table::as.data.table(cnv); sv <- data.table::as.data.table(sv)
  s <- sample %||% unique(cnv$sample)[1]
  cs <- cnv[sample == s]; if (!nrow(cs)) cs <- cnv
  vs <- sv[sample == s]; if (!nrow(vs)) vs <- sv           # single-sample fallback
  ploidy_col <- cs$ploidy
  if (ploidy_mode == "per_chr") {
    base <- cs[, .(b = max(2, round(stats::median(copyNumber[copyNumber < 6], na.rm = TRUE)))), by = seqnames]
    ploidy_col <- base$b[match(cs$seqnames, base$seqnames)]; ploidy_col[is.na(ploidy_col)] <- 2
  }
  cnv_gr <- GenomicRanges::GRanges(cs$seqnames, IRanges::IRanges(cs$start, cs$end), sample = s,
    copyNumber = cs$copyNumber, ploidy = ploidy_col,
    majorAlleleCopyNumber = cs$majorAlleleCopyNumber, minorAlleleCopyNumber = cs$minorAlleleCopyNumber)
  if (!nrow(vs)) {
    bgr <- GenomicRanges::GRanges(); S4Vectors::mcols(bgr) <- S4Vectors::DataFrame(
      WGS_ID = character(), event = character(), svclass = character(), bp_strand = character(),
      PURPLE_AF = numeric(), PURPLE_JCN = numeric(), VF = numeric(), insLen = numeric(),
      HOMLEN = numeric(), PURPLE_CN = numeric())
    return(list(cnv_gr = cnv_gr, breakpoints_gr = bgr))
  }
  ev <- as.character(vs$name)
  mk <- function(ch, p, st) data.frame(chr = .pfx(ch), pos = p, WGS_ID = s, event = ev,
    svclass = vs$svclass, bp_strand = st, PURPLE_AF = NA_real_, PURPLE_JCN = vs$JCN,
    VF = vs$VF, insLen = 0, HOMLEN = ifelse(is.na(vs$homlen), 0, vs$homlen), stringsAsFactors = FALSE)
  bp <- rbind(mk(vs$chrom1, vs$start1, vs$strand1), mk(vs$chrom2, vs$start2, vs$strand2))
  bgr <- GenomicRanges::GRanges(bp$chr, IRanges::IRanges(bp$pos, width = 1), WGS_ID = bp$WGS_ID,
    event = bp$event, svclass = bp$svclass, bp_strand = bp$bp_strand, PURPLE_AF = bp$PURPLE_AF,
    PURPLE_JCN = bp$PURPLE_JCN, VF = bp$VF, insLen = bp$insLen, HOMLEN = bp$HOMLEN)
  win <- GenomicRanges::GRanges(bp$chr, IRanges::IRanges(pmax(1, bp$pos - 1e4), bp$pos + 1e4))
  ov <- GenomicRanges::findOverlaps(win, cnv_gr)
  cnmax <- tapply(cnv_gr$copyNumber[S4Vectors::subjectHits(ov)], S4Vectors::queryHits(ov), max)
  pc <- rep(2, length(bgr)); pc[as.integer(names(cnmax))] <- cnmax; bgr$PURPLE_CN <- pc
  list(cnv_gr = cnv_gr, breakpoints_gr = bgr)
}
