# BKBreed 0.3.0

* New `bk_diallel()` — Griffing (1956) diallel analysis, Method 2 (parents +
  one set of F1s, no reciprocals), Model I (fixed): combining-ability ANOVA
  (GCA, SCA), parental GCA effects, SCA matrix, Baker's predictability ratio and
  variance components. Standard errors are computed exactly by propagating the
  per-entry error variance through the linear estimators (no tabulated
  constants).
* `bk_plot()` gains a `bk_diallel` method: GCA effect bars (default) and an SCA
  heatmap (`type = "sca"`).
* New bundled dataset `bk_data("diallel")` — 6-parent half diallel, 3 reps.

# BKBreed 0.2.0

* New `bk_lxt()` — Line x Tester combining-ability analysis: combining-ability
  ANOVA (crosses partitioned into lines, testers and line x tester), GCA effects
  for lines and testers, SCA effects per cross, proportional contributions, and
  GCA/SCA variance components with additive/dominance variances, average degree
  of dominance and predictability ratio.
* `bk_plot()` gains a `bk_lxt` method: GCA effect bars (default) and an SCA
  heatmap (`type = "sca"`).
* New bundled dataset `bk_data("lxt")` — 5 lines x 3 testers, 3 replications.

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
