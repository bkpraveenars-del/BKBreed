# BKBreed 0.1.0

First release.

* Experimental designs: `bk_rbd()`, `bk_frbd()`, `bk_augmented()`
  (Augmented Alpha-Lattice and Augmented RCBD).
* Biometrical genetics: `bk_variability()`, `bk_correlation()`, `bk_path()`,
  `bk_diversity()` (Mahalanobis D² with Tocher and Ward clustering),
  `bk_stability()` (Eberhart-Russell + AMMI).
* A single `bk_plot()` verb returns a publication-ready ggplot2 figure for
  every result, using the bespoke `theme_bk()` and `bk_palette()` colour system.
* Four bundled example datasets via `bk_data()`.
* Full demonstration script at `inst/examples/run_all.R`.
