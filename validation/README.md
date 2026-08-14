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

On the clean regime shipped here (n = 300, seed 1) the caller separates the
classes perfectly (sensitivity = specificity = PPV = 1.0). **Next step for a
publication-grade ROC:** dial in difficulty — jitter the flank copy number
toward the amplification threshold, add VF noise and CN-segmentation error, and
mislabel a fraction — to show performance degrades gracefully rather than
cliff-edging, and to pick robust defaults for `ext` and the VF rule.

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
