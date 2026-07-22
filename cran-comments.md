## Resubmission (0.3.2)

Addresses the CRAN reviewer's (Konstanze Lauseker) requests:

* All acronyms are now explained in the Description text (e.g. GCV/PCV spelled
  out as genotypic/phenotypic coefficients of variation; GCA/SCA as general/
  specific combining ability; AMMI defined in full).
* Added method references with DOIs in the Description: Griffing (1956)
  <doi:10.1071/BI9560463> and Eberhart and Russell (1966)
  <doi:10.2135/cropsci1966.0011183X000600010011x>.
* inst/examples/run_all.R now writes figures to tempdir() rather than the
  working directory, per CRAN policy.

## Resubmission (0.3.1)

This resubmission fixes the one actionable NOTE from the CRAN incoming
pre-test: an over-long `\usage` line in `bk_diversity.Rd` has been wrapped so
that all Rd usage lines are <= 90 characters.

The remaining NOTE ("Possibly misspelled words in DESCRIPTION": AMMI,
Biometrical, Eberhart, GCA, Griffing, PCV, RBD, SCA, Tocher, biometrical,
diallel, heritability, intra) lists standard plant-breeding and biometry
terms and author surnames, all spelled correctly. They are documented in
inst/WORDLIST.

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

Local (Windows 11, R 4.6.0): 0 errors | 0 warnings | 0 notes.
win-builder R-devel: 0 errors | 0 warnings | 2 notes.

The two NOTEs on win-builder are the expected pair for a new package:

* "New submission" - BKBreed is not yet on CRAN.
* "Possibly misspelled words in DESCRIPTION" - the flagged words
  (biometrical, genotypic, phenotypic, heritability, GCV, PCV, Eberhart,
  Mahalanobis, Tocher, AMMI, GCA, SCA) are standard plant-breeding and
  biometry terms and are spelled correctly.

## Notes for CRAN

* The package uses only base packages plus ggplot2 (Imports); ggrepel and
  patchwork are Suggests and are guarded with requireNamespace().
* No compiled code; examples run quickly and require no internet access.
* Bundled example datasets are simulated and used only for documentation and
  demonstration.
