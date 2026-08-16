## ---------------------------------------------------------------------------
## AmpliconArchitect -> EpiTracer input adapter
##
## Converts a public AmpliconArchitect (AA) reconstruction into the GRanges the
## EpiTracer caller expects, so the cell-line benchmark can run directly on
## AA outputs downloaded from AmpliconRepository (https://ampliconrepository.org,
## 329 CCLE cell lines, CC-BY-4.0) without re-running any WGS pipeline.
##
## AA "_graph.txt" gives us everything EpiTracer needs from one file:
##   * sequence edges  -> copy-number segments (cnv_gr)
##   * discordant edges -> SV breakpoints with read support (breakpoints_gr)
## and "_cycles.txt" intervals give the amplicon footprint (ecdna_gr).
##
## AA graph format (confirmed from the AmpliconSuite docs):
##   sequence   {chr:pos-}  {chr:pos+}  CN  avg_depth  size_bp  total_reads
##   discordant {chr:pos(s)}->{chr:pos(s)}  edge_CN  n_discordant_read_pairs  homlen  homseq
## Vertices are {CHROM}:{POS}{+|-}, POS 0-based.
##
## This file defines functions only; it is sourced by cell_line_benchmark.R.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
})

## --- vertex strand pair -> SV class ----------------------------------------
## The single convention that matters for the episome heuristic is which strand
## pair is the head-to-tail (tandem-duplication / circularisation) junction,
## because EpiTracer keys on a boundary `svclass == "DUP"` carrying the highest
## read support. AA's vertex orientations map to SV type as below.
##
## IMPORTANT — confirm this on a positive control before trusting the numbers.
## Run the benchmark with `selftest`, then on a known ecDNA line (e.g. COLO320DM,
## MYC): the highest-read-support junction spanning the amplicon boundary must
## come out as "DUP". If it comes out "DEL", your AA build uses the opposite
## convention — set `flip_dup_del = TRUE` in aa_to_epitracer_inputs() (a one-line
## switch) and re-run. This is the only convention-dependent step.
aa_svclass <- function(chr1, s1, chr2, s2, flip_dup_del = FALSE) {
  if (!identical(chr1, chr2)) return("TRA")
  key <- paste0(s1, s2)
  cls <- switch(key,
    "-+" = "DUP",     # head-to-tail: tandem dup / circularisation junction
    "+-" = "DEL",     # head-to-head-of-next: deletion / excision scar
    "++" = "h2hINV",
    "--" = "t2tINV",
    "BND")
  if (flip_dup_del) cls <- c(DUP = "DEL", DEL = "DUP")[cls] |> (\(x) ifelse(is.na(x), cls, x))()
  cls
}

## --- parse a vertex "chr:pos[+/-]" -----------------------------------------
.parse_vertex <- function(v) {
  m <- regmatches(v, regexec("^([^:]+):(\\d+)([+-])$", v))[[1]]
  if (length(m) != 4L) stop("unparseable AA vertex: ", v)
  list(chr = m[2], pos = as.numeric(m[3]), strand = m[4])
}

## --- read an AA _graph.txt --------------------------------------------------
read_aa_graph <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  f  <- strsplit(ln, "\t")

  seg <- lapply(f[vapply(f, function(x) x[1] == "sequence", logical(1))], function(x) {
    a <- .parse_vertex(x[2]); b <- .parse_vertex(x[3])
    data.frame(chr = a$chr,
               start = min(a$pos, b$pos), end = max(a$pos, b$pos),
               cn = as.numeric(x[4]),
               total_reads = suppressWarnings(as.numeric(x[7])),
               stringsAsFactors = FALSE)
  })
  seg <- if (length(seg)) do.call(rbind, seg) else
    data.frame(chr = character(), start = numeric(), end = numeric(),
               cn = numeric(), total_reads = numeric())

  is_disc <- vapply(f, function(x) x[1] %in% c("discordant", "concordant"), logical(1))
  brk <- lapply(f[is_disc], function(x) {
    vv <- strsplit(x[2], "->", fixed = TRUE)[[1]]
    a <- .parse_vertex(vv[1]); b <- .parse_vertex(vv[2])
    data.frame(type = x[1],
               chr1 = a$chr, pos1 = a$pos, strand1 = a$strand,
               chr2 = b$chr, pos2 = b$pos, strand2 = b$strand,
               edge_cn = as.numeric(x[3]),
               n_reads = as.numeric(x[4]),
               homlen = suppressWarnings(as.numeric(x[5])),
               stringsAsFactors = FALSE)
  })
  brk <- if (length(brk)) do.call(rbind, brk) else
    data.frame(type = character(), chr1 = character(), pos1 = numeric(),
               strand1 = character(), chr2 = character(), pos2 = numeric(),
               strand2 = character(), edge_cn = numeric(), n_reads = numeric(),
               homlen = numeric())

  list(segments = seg, breaks = brk)
}

## --- read an AA _cycles.txt (Interval lines only) ---------------------------
read_aa_cycles_intervals <- function(path) {
  ln <- readLines(path, warn = FALSE)
  iv <- ln[grepl("^Interval\\b", ln)]
  if (!length(iv)) return(data.frame(chr = character(), start = numeric(), end = numeric()))
  m <- do.call(rbind, lapply(strsplit(iv, "\t"), function(x) {
    # Interval\t{id}\t{chr}\t{start}\t{end}
    data.frame(chr = x[3], start = as.numeric(x[4]), end = as.numeric(x[5]),
               stringsAsFactors = FALSE)
  }))
  m
}

## --- assemble EpiTracer inputs from one AA amplicon --------------------------
## Returns list(ecdna_gr, breakpoints_gr, cnv_gr). Discordant junctions are
## emitted as TWO breakend rows sharing `event` (matching EpiTracer's schema,
## where a DUP appears as its two boundaries). AA gives total CN only, so
## copy number is treated as major-allele (amplicons are near-single-allele);
## minor allele is set to 0.
## `cnv_bed`: optional path to a genome-wide 4-column copy-number BED
##   (chrom, start, end, copyNumber), e.g. the AmpliconSuite CNVkit calls
##   converted to 4 columns. STRONGLY recommended for real data: the episome
##   heuristic requires the segments flanking the amplicon to be diploid, and an
##   AA amplicon graph often contains only the amplified core, so the flank test
##   is undefined without genome-wide CN. When NULL, the AA sequence edges are
##   used (fine for the self-test / amplicons whose graph includes low-CN edges).
aa_to_epitracer_inputs <- function(graph_path, cycles_path, sample_id,
                                    min_break_reads = 5, flip_dup_del = FALSE,
                                    cnv_bed = NULL) {
  g  <- read_aa_graph(graph_path)
  iv <- read_aa_cycles_intervals(cycles_path)

  ## copy-number segments — prefer genome-wide CN when supplied
  if (!is.null(cnv_bed)) {
    bed <- utils::read.delim(cnv_bed, header = FALSE, stringsAsFactors = FALSE)
    ## AmpliconSuite CNVkit calls are "chr start end CNVkit <copyNumber>"; a plain
    ## 4-column bed is "chr start end <copyNumber>". Copy number is the last column.
    cn <- as.numeric(bed[[ncol(bed)]])
    cnv_gr <- GRanges(bed[[1]], IRanges(pmax(1, as.numeric(bed[[2]])), as.numeric(bed[[3]])),
                      sample = sample_id,
                      copyNumber = cn, ploidy = 2,
                      majorAlleleCopyNumber = cn, minorAlleleCopyNumber = 0)
  } else {
    seg <- g$segments
    cnv_gr <- GRanges(seg$chr, IRanges(pmax(1, seg$start), seg$end),
                      sample = sample_id,
                      copyNumber = seg$cn, ploidy = 2,
                      majorAlleleCopyNumber = seg$cn, minorAlleleCopyNumber = 0)
  }

  ## SV breakpoints (intra-chromosomal, above a read-support floor)
  b <- g$breaks
  b <- b[b$type == "discordant" & b$n_reads >= min_break_reads, , drop = FALSE]
  if (nrow(b)) {
    b$svclass <- vapply(seq_len(nrow(b)), function(i)
      aa_svclass(b$chr1[i], b$strand1[i], b$chr2[i], b$strand2[i], flip_dup_del),
      character(1))
    b$event <- paste0(b$svclass, seq_len(nrow(b)))
    mk <- function(chr, pos, row) data.frame(
      chr = chr, pos = pos, WGS_ID = sample_id, event = row$event, svclass = row$svclass,
      PURPLE_AF = NA_real_, PURPLE_JCN = row$edge_cn, VF = row$n_reads,
      PURPLE_CN = row$edge_cn, insLen = 0,
      HOMLEN = ifelse(is.na(row$homlen), 0, abs(row$homlen)), stringsAsFactors = FALSE)
    rows <- do.call(rbind, lapply(seq_len(nrow(b)), function(i)
      rbind(mk(b$chr1[i], b$pos1[i], b[i, ]), mk(b$chr2[i], b$pos2[i], b[i, ]))))
    breakpoints_gr <- GRanges(rows$chr, IRanges(pmax(1, rows$pos), width = 1),
                              WGS_ID = rows$WGS_ID, event = rows$event,
                              svclass = rows$svclass, PURPLE_AF = rows$PURPLE_AF,
                              PURPLE_JCN = rows$PURPLE_JCN, VF = rows$VF,
                              PURPLE_CN = rows$PURPLE_CN, insLen = rows$insLen,
                              HOMLEN = rows$HOMLEN)
  } else {
    ## typed-empty: carry the required mcols so the caller's column check passes
    ## (an amplicon with no qualifying SV simply cannot be episomal).
    breakpoints_gr <- GRanges(seqnames = character(0), ranges = IRanges(integer(0), integer(0)))
    S4Vectors::mcols(breakpoints_gr) <- S4Vectors::DataFrame(
      WGS_ID = character(0), event = character(0), svclass = character(0),
      PURPLE_AF = numeric(0), PURPLE_JCN = numeric(0), VF = numeric(0),
      PURPLE_CN = numeric(0), insLen = numeric(0), HOMLEN = numeric(0))
  }

  ## amplicon footprint (falls back to the CN segment extent if no intervals)
  if (nrow(iv)) {
    ecdna_gr <- GRanges(iv$chr, IRanges(pmax(1, iv$start), iv$end))
  } else if (length(cnv_gr)) {
    ecdna_gr <- range(cnv_gr)
  } else {
    ecdna_gr <- GRanges()
  }
  if (length(ecdna_gr)) {
    ecdna_gr$ID <- paste0(sample_id, "_amp", seq_along(ecdna_gr))
    ecdna_gr$WGS_ID <- sample_id
  }

  list(ecdna_gr = ecdna_gr, breakpoints_gr = breakpoints_gr, cnv_gr = cnv_gr)
}

## --- read AmpliconClassifier classification profiles (silver-standard label) -
## _amplicon_classification_profiles.tsv: sample_name, amplicon_number,
## amplicon_decomposition_class, ecDNA+, BFB+, ...
read_ac_profiles <- function(path) {
  df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(df) <- gsub("[^A-Za-z0-9]+", "_", names(df))
  df
}
