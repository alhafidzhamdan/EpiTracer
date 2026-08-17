# EpiTracer validation

Evidence that the episome heuristic in `call_episomal_ecdna()` is **correct**,
not merely that it runs — the scientific core of the Bioinformatics Application
Note. Three complementary lines; the manuscript figure draws on all three.

This directory is analysis code, not part of the R package (it is
`.Rbuildignore`d). Run each script from the package root.

## 1. Simulation benchmark — `simulate_benchmark.R`

Quantitative sensitivity / specificity / PPV against a **known ground truth**.
Generates focal amplicons of four kinds and scores the caller:

| Simulated class | Truth | Why |
| :--- | :--- | :--- |
| `episomal` | episomal | boundary DUP (highest VF) + diploid flanks + excision scar |
| `gained_flanks` | not episomal | amplified flanks — a complex / chromothriptic context |
| `internal_higher_vf` | not episomal | an internal DUP outranks the boundary DUP |
| `no_boundary_dup` | not episomal | no duplication at the amplicon boundary |

```sh
Rscript validation/simulate_benchmark.R
```

Writes `output/{metrics.csv, confusion.csv, by_type.csv, benchmark.pdf}`.

On the clean regime (n = 300, seed 1) the caller separates the classes perfectly
(sensitivity = specificity = PPV = 1.0).

### Difficulty sweep — `simulate_sweep.R`

Real copy-number and read-support estimates are noisy, so this script adds
increasing Gaussian noise to segment copy number and log-normal noise to the
breakpoint variant fraction, and tracks how the metrics degrade — a tool that
fails gracefully rather than cliff-edging is the point.

```sh
Rscript validation/simulate_sweep.R      # -> output/{sweep_metrics.csv, sweep.pdf}
```

Result (n = 60 per class per level, seed 1): sensitivity stays **1.0 up to
0.25 copies SD** (the realistic PURPLE segment-error regime), degrading
gracefully to 0.87 at 0.5 and 0.68 at 1.0 copies SD — while **specificity and
PPV never drop below 1.0**. EpiTracer is high-precision and conservative under
noise: when it calls an amplicon episomal it is right, and noise costs recall,
not precision.

> Building this sweep **surfaced and fixed a real caller bug**: under noisy /
> tied boundary-DUP VFs the highest-VF test evaluated a length>1 condition and
> errored. Fixed to compare the maximum boundary VF (regression-tested), which
> is exactly the kind of robustness a reviewer probes.

### Runtime / scaling — `runtime_scaling.R`

Wall-clock of `call_episomal_ecdna()` as the number of amplicons grows, for the
manuscript's availability / implementation note.

```sh
Rscript validation/runtime_scaling.R     # -> output/{scaling.csv, scaling.pdf}
```

## 2. Public known-ecDNA cell-line benchmark — `cell_line_benchmark.R`

Runs EpiTracer directly on **public** AmpliconArchitect (AA) reconstructions and
scores it against orthogonal ground truth — no in-house data, fully
reproducible by a reviewer. This is the head-to-head-with-an-existing-tool
evidence the journal expects.

**Files**

| File | Role |
| :--- | :--- |
| `aa_to_epitracer.R` | Adapter: AA `_graph.txt` (sequence edges → CN; discordant edges → SV breakpoints with read support) + `_cycles.txt` intervals → EpiTracer's GRanges inputs. Also reads AmpliconClassifier profiles. |
| `cell_line_panel.tsv` | FISH-validated ground-truth panel (cell line → ecDNA status, oncogene, assay, citation), incl. the matched **COLO320DM (ecDNA⁺) / COLO320HSR (ecDNA⁻)** pair. |
| `cell_line_benchmark.R` | Driver: builds inputs, runs `call_episomal_ecdna()`, scores gold + silver. |
| `fixtures/COLO320DM/` | A synthetic AA-format amplicon used by the self-test. |

**Two evidence tiers**

- **GOLD** — the FISH panel. Does EpiTracer detect the amplicon on ecDNA⁺ lines,
  and does it correctly **not** flag episomal on ecDNA⁻ (HSR / no-amplicon) lines?
- **SILVER** — AmpliconClassifier's own `ecDNA+` / `BFB+` call on every downloaded
  amplicon (large, objective, no manual labels). Detection concordance at scale,
  plus the **episomal fraction among AC-called ecDNA** — the mechanistic layer AC
  does not provide.

**Getting the public data.** AA + AmpliconClassifier results for **329 CCLE cell
lines** are hosted on [AmpliconRepository](https://ampliconrepository.org)
(CC-BY-4.0), downloadable via its command-line API. Place one sub-folder per cell
line under a data root, each containing the AA `*_graph.txt` / `*_cycles.txt` and
the `*_amplicon_classification_profiles.tsv`. Optionally include the AmpliconSuite
CNVkit `*_CNV_CALLS.bed` (see caveat below).

```sh
Rscript validation/cell_line_benchmark.R selftest          # parser + caller smoke test (no download; must print PASS)
Rscript validation/cell_line_benchmark.R /path/to/aa_data  # real benchmark -> output/cell_line_*.tsv
```

**Two subtleties the adapter handles (both found while running real data):**

1. **SV convention must be read position-sorted.** EpiTracer keys on a boundary
   `DUP`. AA lists a discordant edge's two vertices in *either* order, so the SV
   type must be taken from the strands ordered by coordinate, not the raw listing
   order — `(lower −, higher +)` = DUP (circularisation), `(lower +, higher −)` =
   DEL. Verified on GBM39, whose amplicon-spanning circularisation junction is
   written in opposite orders in the hg19 vs hg38 builds yet resolves to DUP in
   both. `flip` on the command line is an escape hatch if a future AA build differs.
2. **Flank test uses the local chromosomal baseline, not global ploidy.** The
   episome heuristic needs the flanks to be non-focally-gained. Cancer
   chromosomes are often polysomic (chr7 in GBM is the textbook case), so a
   genome-wide ploidy of 2 wrongly flags a trisomic flank as "gained". With the
   genome-wide CNVkit `*_CNV_CALLS.bed` (auto-detected per sample), the adapter
   sets each segment's baseline to the per-chromosome median CN, so the flank test
   is calibrated to the local baseline. Without a CNV bed, the AA sequence edges
   are used as a fallback (flanks then unreliable).

**Positive control — GBM39 (canonical EGFR ecDNA).** Run on the hg38
AmpliconSuite output (which ships a CNVkit bed), GBM39's EGFR amplicon is called
**episomal**: an amplicon-spanning circularisation DUP carrying the highest read
support, with flanks at the chr7 (polysomic) baseline. It is *not* callable from
the AA graph alone — chr7's trisomy defeats a global-diploid flank test — which
is exactly why the two fixes above matter.

**Result (329 CCLE lines, AmpliconRepository, run 2026-08-16).** EpiTracer
processed 799 amplicons across 238 lines and called **11 episomal**. Every one
falls in an AmpliconClassifier **Cyclic** (ecDNA) amplicon — **0 of 516**
non-cyclic amplicons (Linear / Complex-non-cyclic / No-FSCNA) were called
episomal (**100% specificity** against the AC class). Within the 283 Cyclic
amplicons, EpiTracer isolates 11 (~3.9%) as *simple episomes* — the mechanism
sub-classification AmpliconClassifier does not make — at canonical ecDNA
oncogenes (MYC in NCI-H2170 [VF 1532] / NCI-H1792 / COR-L279, FGFR2 in KATO III,
MYCN in NCI-H526). See `output/cell_line_benchmark.png`. (`excision_scar` is
`FALSE` throughout: the scar deletion sits in the flank, which the AA amplicon
graph does not include — recoverable only with genome-wide SV calls, e.g.
PURPLE.) ~24 of ~1,140 amplicons (2%) are skipped on degenerate single-region
graphs.

`plot_ccle_episomal.R` draws a montage of **every** amplicon the benchmark
called episomal (one panel each, copy number + SV arcs at the amplicon locus,
oncogene-centred where the panel gene applies). The simple episomes show a single
circularisation arc over diploid flanks (e.g. MYC in MSTO-211H / NCI-H1792 /
NCI-H2170); busier panels carry more internal structure.

```sh
Rscript validation/plot_ccle_episomal.R /path/to/results/samples   # -> output/ccle_episomal_montage.png
```

### Per-locus mechanism classification (`classify_loci.R` + `cell_line_loci_benchmark.R`)

The per-amplicon episomal call is too coarse: an AmpliconArchitect amplicon can
fuse a genuine episome to a BFB or a translocation-driven amplicon (e.g.
NCI-H2170 co-amplifies a chr17 episome and a chr8 MYC amplicon). `classify_loci.R`
classifies **each major amplified locus** by mechanism from the AA graph:

- **episomal** — a boundary DUP whose **both ends reach the edges of the
  amplified region** (a true self-ligated circle, not a DUP that stops short of
  or overshoots the boundary) and is the locus's highest-VF DUP; flanks are
  diploid (vs the per-chromosome baseline); the locus is **not**
  pericentromeric; it is **not** BFB (< 3 near-self fold-back inversions); the
  junctions are **not** a thicket (≤ 3 competing high-copy junctions); and no
  inter-locus TRA exceeds the boundary DUP. (Size is not a criterion — episomes
  range from ~0.2 Mb to > 20 Mb.)
- **BFB** — fold-back inversions (breakage-fusion-bridge), no amplified translocation.
- **translocation-bridge** — an amplified inter-chromosomal translocation forms a
  dicentric bridge that amplifies, with few fold-backs (Lee et al., *Nature* 2023;
  the ERBB2/CCND1 breast-cancer mechanism — e.g. AU565 chr17/ERBB2, BT20 chr7/EGFR).
- **LTA** (loss-translocation-amplification) — an amplified translocation **plus**
  fold-back inversions (BFB cycles); the `arm_loss` column flags the sub-baseline
  loss on the amplicon chromosome (Espejo Valle-Inclán et al., *Cell* 2024).
- **complex** — none of the above (no clean boundary DUP, or too many junctions).

Inter-locus TRAs *below* the per-locus boundary DUPs are read as episome fusion
(the episomal loci are flagged `fused`, e.g. 5637's chr3+chr6 circles).

Calibrated against expert ground truth (5637 chr3+chr6 episomal; KATO III chr3
episomal / chr10 BFB; NCI-H2170 chr17 episomal / chr8 chimeric; NCI-H526 none).
Across the 329 CCLE lines: of **633 major amplified loci**, only **17 (2.7%) in
16 lines are simple episomes** (7 fused) — MYC, CCND2, PDGFRA among them — vs 132
translocation-bridge, 34 LTA, 97 BFB and 353 complex. See
`output/{cell_line_loci.tsv, ccle_locus_mechanisms.png, ccle_episomal_loci_montage.png}`.

```sh
Rscript validation/cell_line_loci_benchmark.R /path/to/results/samples   # -> output/cell_line_loci.tsv
```

A **template** for arbitrary processed inputs (PURPLE + AA, any pipeline) remains
in `cell_line_case_study.R`.

### Positive control — GBM39 EGFR ecDNA (`plot_gbm39_egfr.R`)

The canonical simple EGFR ecDNA line, drawn straight from its public AA
reconstruction (AmpliconRepository *AmpliconSuite benchmarking set hg38*,
project `6a5a4a970664f9111b586742`, sample GBM39). The episome signature is
visible at a glance — a single **circularisation junction (DUP) spanning the
amplified EGFR locus**, diploid flanks either side. EpiTracer calls it episomal
once the per-chromosome baseline is used (chr7 is polysomic in GBM; a
global-diploid flank test misses it — see §2 subtleties).

```sh
Rscript validation/plot_gbm39_egfr.R /path/to/results/samples/GBM39   # -> output/gbm39_egfr.png
```

### Engineered ground-truth episome — time course (`engineered_ecdna_timecourse.R`)

The strongest positive control: an ecDNA formed by a **designed** excision-
circularisation event, sampled over time. Pradella et al. engineered a Cre-loxP
cassette in mouse NPCs that circularises the **Myc/Pvt1** locus; WGS was taken
at 1-5 weeks, with a p53-Cre arm as control (AmpliconRepository *Pradella et al.
engineered murine Myc ecDNA*, project `6a5ea2166714848c8db0f348`, mm10).

EpiTracer calls the locus **episomal exactly at 4w and 5w** — when the circle
forms — and never in the control arm or the earlier timepoints, tracking the
locus copy number as it rises (0.5 → 3.6 → 38 → 111 → 126 over the five weeks).
At 3w the locus is already amplified (CN 38, AmpliconClassifier "Linear") but is
**not** called episomal: EpiTracer waits for the circle, calling an *episome*
rather than mere amplification. The 4w and 5w calls share the identical
circularisation junction (chr15 ~60.65 Mb ↔ ~62.38 Mb).

```sh
Rscript validation/engineered_ecdna_timecourse.R /path/to/results/samples  # -> output/engineered_timecourse.{tsv,png}
```

## 3. In-house WGS cohort (script to add)

Across the cohort, show EpiTracer's `episomal` calls are a well-defined **subset**
of AmpliconClassifier "ecDNA", enriched for a single high-VF boundary DUP,
non-gained flanks and excision-scar deletions versus complex amplicons —
demonstrating EpiTracer resolves a *mechanism* the standard tool does not.

## Mapping to the manuscript

- **Figure panel B** — simulation performance (`benchmark.pdf`).
- **Figure panel C** — a real episomal ecDNA from a confirmed cell line
  (`plot_sv_linear` / `plot_sv_reconstruction`).
- **Results text** — cohort concordance + enrichment (line 3).
