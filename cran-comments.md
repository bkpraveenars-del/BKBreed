## Submission summary

This is a new submission of BKBreed, a small, dependency-light toolkit for
biometrical analysis of plant-breeding and genetics field experiments
(experimental designs, genetic variability, correlation and path analysis,
Line x Tester and Griffing diallel combining ability, D-square divergence, and
Eberhart-Russell / AMMI stability), with publication-ready ggplot2 figures.

## Test environments

* Local: Windows 11 x64, R 4.6.0 -- R CMD check --no-manual: OK
  (0 errors | 0 warnings | 0 notes)
* (before submission) win-builder devel and release via
  devtools::check_win_devel()
* (before submission) R-hub multi-platform check

## R CMD check results

0 errors | 0 warnings | 0 notes on the local run.

On first CRAN submission a single NOTE is expected:
  "New submission" (this package is not yet on CRAN).

## Notes for CRAN

* The package uses only base packages plus ggplot2 (Imports); ggrepel and
  patchwork are Suggests and are guarded with requireNamespace().
* No compiled code; examples run quickly and require no internet access.
* Bundled example datasets are simulated and used only for documentation and
  demonstration.
