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
    * **Bug fix:** the `chromosome_range` (zoom) code path was broken because
      `chr_selection`'s chromosome column was mis-named `chrs` instead of `chr`;
      zoomed views now correctly restrict to the requested window.
    * Now **returns the `ggplot` object** (with any written PDF path attached as
      attribute `"path"`); writing a PDF is optional (`outdir`/`save`).
    * **Robust axes:** x-axis tick spacing auto-adapts to the plotted window
      (`scale_ticks = NULL`); the copy-number axis auto-scales its breaks to the
      data (fixing label collisions at very high copy number); the read-support
      secondary axis is placed at real VF values and omitted when SV arcs are
      negligible relative to a focal amplicon.
    * **Flexible inputs:** `wgd_data` may be keyed by `sample` or `WGS_ID` (or a
      custom `wgd_sample_col`); the output filename encodes the plotted region so
      zoomed plots no longer overwrite full-chromosome ones.
    * **Polish:** gene labels are de-collided with \pkg{ggrepel} (optional,
      `repel_labels`); line geoms use `linewidth` (no more ggplot2 deprecation
      warnings); debug output is gated behind `verbose`.
