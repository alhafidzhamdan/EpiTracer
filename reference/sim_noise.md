# Sequencing and calling noise for simulated amplicons

Bundles the technical parameters
[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md)
uses to turn ground truth into the imperfect observation a short-read
pipeline reports. Passing a scalar to `sim_noise()` scales every noise
term at once, giving a one-knob difficulty sweep; passing named
arguments sets terms individually.

## Usage

``` r
sim_noise(
  scale = 1,
  purity = 0.8,
  depth = 60,
  cn_sd = 0.25,
  cn_cv = 0.02,
  vf_cv = 0.1,
  dropout = 0.03,
  min_vf = 4,
  fp_rate = 0.5,
  min_seg_width = 10000,
  bp_jitter = 0
)
```

## Arguments

- scale:

  Numeric multiplier applied to every stochastic term (default `1`).
  `sim_noise(0)` is a noiseless, exact rendering of the ground truth;
  `sim_noise(2)` roughly doubles the segmentation error, the
  read-support dispersion and the junction dropout and false-positive
  rates.

- purity:

  Tumour purity (default `0.8`). Lowers junction variant allele fraction
  and read support.

- depth:

  Mean sequencing depth in x (default `60`), the usual WGS target. Read
  support for a junction scales with `depth * purity`.

- cn_sd:

  Copy-number segmentation error in copies, absolute component (default
  `0.25`; PURPLE's segment-level error is typically 0.2-0.5).

- cn_cv:

  Proportional component of the segmentation error (default `0.02`), so
  high-copy segments carry proportionally more error.

- vf_cv:

  **Extra-Poisson** dispersion of junction read support (default `0.1`).
  Read support is a fragment COUNT, so its dominant noise is counting
  noise: support is drawn from a negative binomial whose variance is
  `mu + (vf_cv * mu)^2`. A weakly supported junction is therefore noisy
  in relative terms (at `mu = 12`, CV about 0.30, Poisson-dominated)
  while a strongly supported one is tight (at `mu = 3000`, CV about
  0.10, set by `vf_cv`, which stands for local coverage variation – GC,
  mappability). Modelling this as a constant-CV log-normal instead would
  over-disperse high-support junctions by an order of magnitude and hide
  the fact that support is proportionate to junction copy number.

- dropout:

  Probability a true junction is missed, before the additional
  support-dependent dropout of weakly supported junctions (default
  `0.03`).

- min_vf:

  Junctions whose simulated read support falls below this are dropped as
  unsupported (default `4`), the usual caller floor.

- fp_rate:

  Expected number of false-positive junctions added per amplicon
  (default `0.5`).

- min_seg_width:

  Copy-number segments narrower than this are merged into their
  neighbours (default `1e4`), modelling the resolution limit of
  read-depth segmentation – without it a heavily shattered circle emits
  hundreds of sub-kb segments no real caller would resolve.

- bp_jitter:

  Positional error on reported breakend coordinates, in bp (default `0`;
  set e.g. `10` to model imprecise breakpoints).

## Value

A named list of noise parameters, of class `sim_noise`.

## See also

[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md)

## Examples

``` r
sim_noise(0)          # exact ground truth
#> $purity
#> [1] 0.8
#> 
#> $depth
#> [1] 60
#> 
#> $cn_sd
#> [1] 0
#> 
#> $cn_cv
#> [1] 0
#> 
#> $vf_cv
#> [1] 0
#> 
#> $dropout
#> [1] 0
#> 
#> $min_vf
#> [1] 4
#> 
#> $fp_rate
#> [1] 0
#> 
#> $min_seg_width
#> [1] 10000
#> 
#> $bp_jitter
#> [1] 0
#> 
#> attr(,"class")
#> [1] "sim_noise"
sim_noise(purity = 0.4, depth = 30)   # a low-purity, shallow sample
#> $purity
#> [1] 0.4
#> 
#> $depth
#> [1] 30
#> 
#> $cn_sd
#> [1] 0.25
#> 
#> $cn_cv
#> [1] 0.02
#> 
#> $vf_cv
#> [1] 0.1
#> 
#> $dropout
#> [1] 0.03
#> 
#> $min_vf
#> [1] 4
#> 
#> $fp_rate
#> [1] 0.5
#> 
#> $min_seg_width
#> [1] 10000
#> 
#> $bp_jitter
#> [1] 0
#> 
#> attr(,"class")
#> [1] "sim_noise"
```
