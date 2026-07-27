# EpiTracer 0.0.0.9000

* Initial development version.
* `call_episomal_ecdna()`: episomal (breakage-independent) ecDNA caller,
  refactored from the original cohort analysis script into a function taking
  ecDNA / SV / CNV / cancer-gene `GRanges` as arguments.
    * **Bug fix:** the original script's boundary-flank check compared
      `after_dist_boundary_not_gained` to itself (always `TRUE`). It now
      correctly requires *both* flanks to be non-gained before calling an
      episome region.
* `plot_sv_linear()`: linear allele-specific copy-number / SV "recon" plotter.
  Reference annotation (karyotype ideogram, gene coordinates, CDS models) is now
  passed via the `karyotype`, `gene_coord`, and `cds_gr` arguments instead of
  being read from hard-coded paths. `useDingbats = FALSE` for modern R.
