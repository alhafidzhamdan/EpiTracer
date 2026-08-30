# Calling amplicon-formation mechanisms

[`call_simple_excision()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
asks *“is this focal amplification a simple-excision circle (an
episome)?”*. The **mechanism callers** ask the complementary question —
*by which route did the amplification form?* Each one takes the **same
inputs** as the episomal caller (`breakpoints_gr`, `cnv_gr`,
`cancer_genes_gr`, and an optional `ecdna_gr`) and returns the input
breakpoints annotated with one mechanism’s flag. They are
**independent**: an amplicon can carry more than one flag (e.g. an
episome that has since shattered is both `episomal` *and*
`chromothripsis`). Run several and join their outputs by `WGS_ID` + `ID`
to assemble a combined mechanism table.

Every example below builds a **small synthetic amplicon** that exhibits
one mechanism’s structural signature — the same minimal inputs the
package’s tests use — then runs the caller. Real inputs come from PURPLE
/ GRIDSS-style calls; see
[`vignette("epitracer")`](https://alhafidzhamdan.github.io/EpiTracer/articles/epitracer.md)
for the input format. A shared cancer-gene locus (EGFR on chr7)
annotates the calls:

``` r

cancer_genes_gr <- GRanges("chr7", IRanges(55019017, 55211628), gene = "EGFR")
```

## Chromothripsis within an ecDNA

Chromothripsis shatters a region and rejoins the fragments at random.
[`call_chromothripsis()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.md)
scores ShatterSeek-style hallmarks over the amplified footprint: enough
internal breakpoints, the four junction orientations occurring about
equally (random joins), and an **oscillating** copy-number profile. It
is meant to be run on amplicons already called episomal, to flag the
subset that have subsequently been through a micronucleus.

Here a boundary circularisation DUP encloses a footprint scattered with
balanced-orientation internal SVs over an oscillating CN profile:

``` r

lo <- 55000000L; hi <- 55500000L
ecdna_gr <- GRanges("chr7", IRanges(lo, hi), ID = "S1_amp1", WGS_ID = "S1")

## boundary circularisation DUP + 12 internal SVs (3 each of DEL/DUP/h2h/t2t INV)
bp <- data.table(seqnames = "chr7", start = c(lo, hi), end = c(lo, hi),
  WGS_ID = "S1", event = "bDUP", svclass = "DUP",
  PURPLE_AF = 0.9, PURPLE_JCN = 40, VF = 1000, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L)
classes <- c("DEL", "DUP", "h2hINV", "t2tINV"); k <- 0L
for (cls in classes) for (i in 1:3) {
  k <- k + 1L; a <- lo + as.integer((hi - lo) * k / 13); b <- a + 30000L
  bp <- rbind(bp, data.table(seqnames = "chr7", start = c(a, b), end = c(a, b),
    WGS_ID = "S1", event = paste0("INT", k), svclass = cls,
    PURPLE_AF = 0.9, PURPLE_JCN = 20, VF = 500, PURPLE_CN = 40, insLen = 0L, HOMLEN = 0L))
}
## oscillating copy number across the footprint (turning points), diploid flanks
osc <- c(lo, lo + 1e5L, lo + 2e5L, lo + 3e5L, lo + 4e5L)
cnv <- data.table(seqnames = "chr7",
  start = c(1L, 40000000L, osc, hi + 1L),
  end   = c(39999999L, lo - 1L, c(osc[-1] - 1L, hi), 159000000L),
  sample = "S1", copyNumber = c(2, 2, 40, 18, 42, 16, 40, 2), ploidy = 2,
  majorAlleleCopyNumber = 1, minorAlleleCopyNumber = 1)

res <- call_chromothripsis(ecdna_gr, gr(bp), gr(cnv), cancer_genes_gr, mc.cores = 1)
unique(res[, .(chromothripsis, chromothripsis_conf, n_internal_sv, cn_oscillations)])
#>    chromothripsis chromothripsis_conf n_internal_sv cn_oscillations
#>            <char>              <char>         <int>           <int>
#> 1:           TRUE                high            13               3
```

A clean simple episome (only a boundary junction, no internal
shattering) returns `chromothripsis = "FALSE"`, and a
fold-back–dominated footprint fails the random-joins test — the
discrimination the caller is built for.

## Micronucleation — two ecDNAs fused in a micronucleus

[`call_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_micronucleation.md)
flags the **two-ecDNA** signature: an amplicon joined to *another
amplified locus on a non-homologous chromosome* by a high-VF
interchromosomal translocation (both breakends amplified). Fusing
fragments of two different chromosomes into one amplicon is only
possible if both were present together — the signature of two episomal
ecDNAs co-encapsulated in a micronucleus and recombined after
shattering. (Micronucleation of a single or homologous chromosome
instead presents as *intrachromosomal* shattering — use
[`call_chromothripsis()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromothripsis.md)
for that.)

``` r

ecdna_gr <- GRanges("chr7", IRanges(55000000, 55500000), ID = "S1_amp1", WGS_ID = "S1")

## a high-VF chr7<->chr12 TRA, with the chr12 partner also amplified
bp <- data.table(
  seqnames = c("chr7", "chr12", "chr7", "chr7"),
  start = c(55200000, 57200000, 55100000, 55400000),
  end   = c(55200000, 57200000, 55100000, 55400000),
  WGS_ID = "S1", event = c("TRA1", "TRA1", "D1", "D1"),
  svclass = c("TRA", "TRA", "DEL", "DEL"),
  PURPLE_AF = 0.9, PURPLE_JCN = c(50, 50, 10, 10), VF = c(3000, 3000, 100, 100),
  PURPLE_CN = c(50, 50, 50, 50), insLen = 0L, HOMLEN = 0L)      # chr12 partner amplified
cnv <- data.table(seqnames = "chr7", start = c(1, 55000000, 55500001),
  end = c(54999999, 55500000, 70000000), sample = "S1", copyNumber = c(2, 50, 2),
  ploidy = 2, majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)

res <- call_micronucleation(ecdna_gr, gr(bp), gr(cnv), cancer_genes_gr, ext = 1e7, mc.cores = 1)
unique(res$micronucleation)
#> [1] "TRUE"
```

A translocation to a **non-amplified** partner (a passenger, not a
second circle) returns `"FALSE"`.

## Breakage–fusion–bridge (BFB)

Classical BFB is **telomeric**: a chromosome loses its telomere, sister
chromatids fuse, and the dicentric bridge breaks unevenly each cell
cycle, building a fold-back amplicon that is *terminal on one arm*, has
a *distal deletion to the telomere*, and shows a copy-number
**staircase**.
[`call_bfb()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_bfb.md)
gates on all three, which separates iterative BFB from a single
fold-back or a focal ecDNA. It needs centromere and chromosome-length
references to locate the telomere:

``` r

centromeres   <- load_centromeres("hg38")
chrom_lengths <- load_chrom_lengths("hg38")

lo <- 10e6; up <- 20e6; w <- up - lo
ecdna_gr <- GRanges("chr7", IRanges(lo, up), ID = "S1_amp1", WGS_ID = "S1")

## two fold-back inversions spread across the amplicon (staircase spread)
bp <- data.table(seqnames = "chr7", start = c(lo + 0.2 * w, lo + 0.8 * w),
  end = c(lo + 0.2 * w, lo + 0.8 * w), WGS_ID = "S1", event = c("FB1", "FB2"),
  svclass = "h2hINV", PURPLE_AF = 0.9, PURPLE_JCN = 30, VF = 500, PURPLE_CN = 40,
  insLen = 0L, HOMLEN = 0L)
## terminal deletion [0, lo) + a stepped CN staircase across the amplicon
cnv <- data.table(seqnames = "chr7",
  start = c(1, lo, lo + w/3 + 1, lo + 2*w/3 + 1, up + 1),
  end   = c(lo - 1, lo + w/3, lo + 2*w/3, up, 159000000),
  sample = "S1", copyNumber = c(1, 40, 25, 12, 2), ploidy = 2,      # distal_cn = 1 (deleted)
  majorAlleleCopyNumber = c(0, 39, 24, 11, 1), minorAlleleCopyNumber = c(0, 1, 1, 1, 1))

res <- call_bfb(ecdna_gr, gr(bp), gr(cnv), cancer_genes_gr,
                centromeres = centromeres, chrom_lengths = chrom_lengths, mc.cores = 1)
unique(res[, .(bfb, bfb_anchor, n_foldbacks)])
#>       bfb bfb_anchor n_foldbacks
#>    <char>     <char>       <int>
#> 1:   TRUE p-telomere           2
```

Drop any one hallmark — a diploid terminus (no distal deletion), a flat
plateau (no staircase), or clustered fold-backs — and the call reverts
to `"FALSE"`.

## Breakage–replication/fusion (BRF)

BRF (Mendez-Dorantes / Zhang / Pellman, *Nat Genet* 2026) leaves
**adjacent parallel breakpoints**: two same-orientation breakends, close
together, from *distinct* junctions.
[`call_brf()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_brf.md)
finds them via
[`find_parallel_breakpoints()`](https://alhafidzhamdan.github.io/EpiTracer/reference/find_parallel_breakpoints.md).

``` r

ecdna_gr <- GRanges("chr7", IRanges(55000000, 55500000), ID = "S1_amp1", WGS_ID = "S1")

## J1 @ 55.100 Mb and J2 @ 55.110 Mb: two + breakends 10 kb apart, distinct events
bp <- data.table(seqnames = "chr7",
  start = c(55100000, 55110000, 55400000, 55600000),
  end   = c(55100000, 55110000, 55400000, 55600000),
  WGS_ID = "S1", event = c("J1", "J2", "J1", "J2"), svclass = "h2hINV",
  bp_strand = "+", PURPLE_AF = 0.9, PURPLE_JCN = 30, VF = 500, PURPLE_CN = 50,
  insLen = 0L, HOMLEN = 0L)
cnv <- data.table(seqnames = "chr7", start = c(1, 55000000, 55500001),
  end = c(54999999, 55500000, 70000000), sample = "S1", copyNumber = c(2, 50, 2),
  ploidy = 2, majorAlleleCopyNumber = c(1, 49, 1), minorAlleleCopyNumber = 1)

res <- call_brf(ecdna_gr, gr(bp), gr(cnv), cancer_genes_gr, ext = 1e7, mc.cores = 1)
unique(res[, .(brf, n_parallel_pairs)])
#>       brf n_parallel_pairs
#>    <char>            <int>
#> 1:   TRUE                1
```

## Chromoplexy

Chromoplexy (Baca *et al.* 2013) is a **closed, balanced** chain of
rearrangements that cycles through several chromosomes and returns to
its start, at *low* copy number (it is balanced, not amplified). It is
genome-wide rather than amplicon-anchored, so
[`call_chromoplexy()`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_chromoplexy.md)
takes just `breakpoints_gr` and `cnv_gr` and returns **one row per
chain** it finds:

``` r

## a closed cycle chr1 -> chr5 -> chr12 -> chr1, each hop a deletion bridge, at CN 2
jun <- data.frame(
  chr1 = c("chr1", "chr5", "chr12"),  pos1 = c(10000000, 20005000, 30005000),
  chr2 = c("chr5", "chr12", "chr1"),  pos2 = c(20000000, 30000000, 10005000))
ev <- paste0("J", 1:3)
bp <- data.table(seqnames = c(jun$chr1, jun$chr2),
  start = as.integer(c(jun$pos1, jun$pos2)), end = as.integer(c(jun$pos1, jun$pos2)),
  WGS_ID = "S1", event = rep(ev, 2), svclass = "TRA", PURPLE_CN = 2)
cnv <- data.table(seqnames = c("chr1", "chr5", "chr12"), start = 1L, end = 2e8L,
  sample = "S1", copyNumber = 2, ploidy = 2,
  majorAlleleCopyNumber = 1, minorAlleleCopyNumber = 1)

res <- call_chromoplexy(gr(bp), gr(cnv))
res[, .(chromoplexy_id, topology, n_junctions, n_chromosomes, chromosomes, frac_cp)]
#>    chromoplexy_id topology n_junctions n_chromosomes     chromosomes frac_cp
#>             <int>   <char>       <int>         <int>          <char>   <num>
#> 1:              1    cycle           3             3 chr1,chr12,chr5       1
```

An **open** chain (ends not rejoined) or an **amplified** cycle (high
copy number — an amplicon, not a balanced event) is not called.

## Combining the callers

Because each caller reports one signature and leaves the others
untouched, the recommended workflow is to run the relevant callers on
the same inputs and join their per-breakpoint outputs by `WGS_ID` +
`ID`. An amplicon then carries an explicit flag for every mechanism
tested — `episomal`, `chromothripsis`, `micronucleation`, `bfb`, `brf` —
rather than being forced into a single class. See
[`?call_simple_excision`](https://alhafidzhamdan.github.io/EpiTracer/reference/call_simple_excision.md)
for the shared input contract, and
[`vignette("reconstruction")`](https://alhafidzhamdan.github.io/EpiTracer/articles/reconstruction.md)
for reading the founder junction off a VF-stratified reconstruction.
