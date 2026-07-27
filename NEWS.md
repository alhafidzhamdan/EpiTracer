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
    * **General CN/SV viewer, not amplification-only.** New `events` argument
      (`"amp"`, `"gain"`, `"hetdel"`, `"homdel"`) chooses which copy-number events
      auto-detection targets; explicit `chromosome`/`loci` plot any region
      regardless of event type. Auto-detected events are **clustered** into
      separate loci (`cluster_gap`) so scattered focal events become their own
      panels. Thresholds exposed via `min_cn_ratio`/`gain_ratio`/`hetdel_ratio`/
      `homdel_thresh`.
    * `displayExon` (with `cds_gr`) re-introduced -- draws exon models for
      in-window genes.
    * Larger axis fonts and longer tick marks; ideogram 20% slimmer; a leading
      `0.0` Mb tick at a locus start is no longer drawn.
    * **Unified onto a single concatenated x-axis** (faceting retired): draws one
      locus or many side-by-side, and structural variants interconnecting separate
      amplicons are drawn as arcs spanning the loci (previously impossible under
      faceting -- see the interim `plot_amplicon_recon()`, now merged in and
      removed). Loci come from `chromosome`+`chromosome_range`, explicit `loci`, or
      are **auto-detected** from amplified segments (`copyNumber > min_cn_ratio x
      ploidy`, padded by `margin`) when neither is given.
    * Removed faceting-only arguments (`scale_ticks`, `yend_left`/`yend_right`,
      `interchrom_arcs`, `displayExon`/`cds_gr`); axis tops/ticks are now automatic.
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
    * The **read-support (SV) secondary axis is always shown** alongside the
      allele-specific copy-number axis, **sharing the same three tick positions**
      (0, half, full). Axis tops round to two significant figures (e.g. 132 -> 140,
      1462 -> 1500); both axis spines are black and stop at the baseline so the
      ideogram gap stays clear.
    * New `interchrom_arcs` argument. Inter-chromosomal connecting arcs reach
      across panels via a large offset that stretches the source panel under
      free-x faceting; set `interchrom_arcs = FALSE` for multi-chromosome plots
      to keep each panel bounded to its locus (breakpoints still marked).
    * **Performance:** SV arcs and breakpoint lines are now drawn as a handful
      of batched layers (grouped by curvature/linewidth via `scale_colour_identity`)
      instead of one ggplot layer per SV. A highly-rearranged amplicon with ~300
      SVs now renders in ~2 s instead of ~2 min, with identical output.
    * Writes a **PDF by default** (`format = "pdf"`); pass `format = c("pdf", "png")` or `"png"` to also/only write a high-resolution PNG (`dpi`).
    * **LOH / homozygous-deletion bars** are placed in the empty band between
      the ideogram and the copy-number baseline (`loh_position_ratio` now a
      fraction of that gap), so they no longer overlap the ideogram.
    * Gene labels **angle automatically (45 deg) when crowded** within a panel
      and stay horizontal otherwise; overridable via `gene_label_angle`. When no
      genes are specified, all default oncogenes falling in the plotted window
      are shown.
