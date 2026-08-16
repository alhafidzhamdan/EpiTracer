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

**Two things to verify on first real run (documented in `aa_to_epitracer.R`):**

1. **Strand convention.** EpiTracer keys on a boundary `DUP`. The AA vertex
   strand → SV-type map is set for the common convention; confirm on a positive
   control (COLO320DM MYC boundary must come out `DUP`). If it comes out `DEL`,
   append `flip` to the command — a one-line switch.
2. **Genome-wide copy number.** The episome heuristic requires the segments
   *flanking* the amplicon to be diploid, but an AA amplicon graph often contains
   only the amplified core. Supply genome-wide CN (the CNVkit `*_CNV_CALLS.bed`,
   auto-detected in each sample folder) so the flank test is well-defined;
   otherwise the AA sequence edges are used as a fallback.

**Result (329 CCLE lines, AmpliconRepository, run 2026-08-16).** EpiTracer
processed 798 amplicons across 237 lines and called **13 episomal**. Every one
falls in an AmpliconClassifier **Cyclic** (ecDNA) amplicon — **0 of 516**
non-cyclic amplicons (Linear / Complex-non-cyclic / No-FSCNA) were called
episomal (100% specificity against the AC class). Within the 282 Cyclic
amplicons, EpiTracer isolates 13 (~4.6%) as *simple episomes* — the
mechanism sub-classification AmpliconClassifier does not make — at canonical
ecDNA oncogenes (MYC in NCI-H2170/NCI-H1792/COR-L279, FGFR2 in KATO III, MYCN
in NCI-H526). See `output/cell_line_benchmark.png`. (`excision_scar` is `FALSE`
throughout: the scar deletion sits in diploid flanks that the AA amplicon graph
does not include — recoverable only with genome-wide SV calls, e.g. PURPLE.)
~25 of ~1,116 amplicons (2%) are skipped on degenerate single-region graphs.

A **template** for arbitrary processed inputs (PURPLE + AA, any pipeline) remains
in `cell_line_case_study.R`.

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
