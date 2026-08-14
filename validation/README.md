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

## 2. Known-ecDNA cell lines — `cell_line_case_study.R`

Ground-truth **case study** on cell lines whose ecDNA status at a focal oncogene
is established orthogonally (FISH / metaphase imaging / AmpliconArchitect):
EGFR (e.g. GBM39), MYCN (e.g. IMR-32, Kelly), MYC (e.g. COLO320-DM). Confirms
EpiTracer flags the episomal lines and reconstructs the circle, and correctly
does **not** call chromothriptic-hub ecDNA episomal. A **template** — edit the
input paths to your processed cell-line PURPLE + AmpliconArchitect output.

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
