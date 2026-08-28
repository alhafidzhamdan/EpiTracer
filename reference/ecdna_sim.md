# Simulated ecDNA amplicon

The object returned by
[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md)
and transformed by the mechanism operators
([`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_fuse_episomes()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_fuse_episomes.md),
[`sim_evolve()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_evolve.md)).
It records the *ground truth* structure of one simulated amplicon, from
which
[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md)
renders caller inputs and
[`sim_to_plot_inputs()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_plot_inputs.md)
renders plotter inputs.

## Components

- `sample`:

  Character sample identifier.

- `fragments`:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  in **amplicon order** – the sequence of intervals a polymerase would
  traverse walking the circle once. Metadata `inv` (`0` forward, `1`
  inverted) and `origin` (the label of the parental episome a fragment
  came from, which distinguishes the two halves of a chimera).

- `circular`:

  `TRUE` for an episome; the last fragment joins the first.

- `copies`:

  Numeric ecDNA copies per tumour cell. Copy number is
  `copies * coverage(fragments)` on top of the host baseline, so this
  sets the amplicon's copy-number level while `fragments` sets its
  shape.

- `junctions`:

  A
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  ledger of every junction with `chrom1`, `start1`, `strand1`, `chrom2`,
  `start2`, `strand2`, `svclass`, `multiplicity` (traversals per
  circle), `round`, `mechanism`, `origin` (`"amplicon"` or `"host"`) and
  `homology` (the repair pathway).

- `host`:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of host-chromosome background with `ploidy`, giving the local
  copy-number baseline the flanks sit at.

- `variants`:

  `NULL` for a homogeneous amplicon. Once a single molecule is altered
  (see
  [`sim_internal_deletion()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_internal_deletion.md))
  the cell carries a MIXTURE, and this holds it as a list of
  `list(fragments, copies, label)`. `fragments` and `copies` above then
  report the most abundant species and the total. The mixture is what
  makes the ORDER of events observable: a junction that arose late is
  carried by fewer circles, so it is emitted at lower copy number and
  lower read support. Each species also carries a `fitness`, its
  replication rate relative to the founder, which is what lets a
  minority species expand rather than sit at a frozen share.

- `preserved`:

  Oncogene loci that mechanisms must not delete.

- `history`:

  Character log of the mechanism steps applied, in order.

- `round`:

  Integer count of evolutionary rounds applied so far.

## See also

[`sim_episome()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_episome.md),
[`sim_micronucleation()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_micronucleation.md),
[`sim_to_epitracer()`](https://alhafidzhamdan.github.io/EpiTracer/reference/sim_to_epitracer.md)
