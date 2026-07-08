# BKBreed

**Colourful Biometrical Analysis for Plant Breeding and Genetics**

BKBreed is a small, accurate R package that runs the everyday analyses of a
plant-breeding programme in one call each — and returns a **publication-ready,
stunningly colourful `ggplot2` figure** for every result through a single verb:
`bk_plot()`.

It was built to be *different* from the existing agricultural-statistics
packages (e.g. `agricolae`, `AgroR`, `metan`, `aridagri`): those are broad and
console-first; BKBreed is **genetics-focused, figure-first, and covers the
augmented (unreplicated) designs that balanced-only packages cannot.**

---

## What it does

| Domain | Function | Signature figure |
|---|---|---|
| Randomised Block Design | `bk_rbd()` | ranked genotype means + CD letters |
| Factorial RBD | `bk_frbd()` | interaction line plot |
| **Augmented Alpha-Lattice / RCBD** | `bk_augmented()` | adjusted-means lollipop vs check line |
| Genetic variability (GCV, PCV, h², GA) | `bk_variability()` | GCV vs PCV bars + heritability |
| Genotypic & phenotypic correlation | `bk_correlation()` | diverging correlation heatmap |
| Path-coefficient analysis | `bk_path()` | direct/indirect effect matrix |
| Line × Tester combining ability | `bk_lxt()` | GCA effect bars + SCA heatmap |
| Griffing diallel (Method 2, Model I) | `bk_diallel()` | GCA effect bars + SCA heatmap |
| Mahalanobis D² divergence | `bk_diversity()` | D² principal-coordinate clusters |
| G×E stability (Eberhart-Russell + AMMI) | `bk_stability()` | AMMI-2 biplot & E-R plot |

Every function also returns a tidy result object with a formatted `print()`
method (ANOVA tables with SS/MS/F/p, SE, CD, CV%, letter groupings).

---

## Install

From the package folder (source):

```r
# option A — install the built tarball
install.packages("BKBreed_0.1.0.tar.gz", repos = NULL, type = "source")

# option B — install directly from GitHub
# remotes::install_github("bkpraveenars-del/BKBreed")
```

Dependencies: `ggplot2` (required); `ggrepel`, `patchwork` (optional, for nicer
labels). R ≥ 4.0.

---

## 60-second tour

```r
library(BKBreed)

## 1. RBD — one call gives ANOVA + CD + letters + a colour figure
rbd <- bk_rbd(bk_data("rbd"), trait = "grain_yield",
              gen = "genotype", rep = "rep")
rbd            # formatted ANOVA + ranked means
bk_plot(rbd)   # ranked colour bars with SE and letter groups

## 2. Genetic variability across several traits
traits <- c("grain_yield","plant_height","tillers","panicle_len","test_weight")
bk_plot(bk_variability(bk_data("rbd"), traits, "genotype", "rep"))

## 3. Correlation heatmap (genotypic)
bk_plot(bk_correlation(bk_data("rbd"), traits, "genotype", "rep"))

## 4. Path analysis on yield
bk_plot(bk_path(bk_data("rbd"),
        c("plant_height","tillers","panicle_len","test_weight"),
        dependent = "grain_yield", gen = "genotype", rep = "rep"))

## 5. D² divergence clusters
bk_plot(bk_diversity(bk_data("rbd"), traits, "genotype", "rep"))

## 6. Augmented alpha-lattice (unreplicated test entries + checks)
bk_plot(bk_augmented(bk_data("augmented"), "grain_yield", "genotype",
                     block = "block", rep = "rep",
                     checks = paste0("CHK-", 1:4)))

## 7. Multi-location stability — AMMI-2 biplot
bk_plot(bk_stability(bk_data("mlt"), "grain_yield",
        gen = "genotype", env = "environment", rep = "rep"))
```

Run everything and save all figures as PNG:

```r
source(system.file("examples", "run_all.R", package = "BKBreed"))
```

---

## Bundled example data

`bk_data("rbd")` pearl-millet multi-trait RBD · `bk_data("frbd")` N × variety
factorial · `bk_data("augmented")` augmented alpha-lattice (4 checks + 20 test
entries) · `bk_data("mlt")` 10 genotypes × 5 environments.

## Methods & conventions

- Variability uses the Singh & Chaudhary formulation
  (σ²g = (Mg − Me)/r, σ²p = σ²g + Me).
- Correlations use analysis-of-covariance mean cross-products.
- Augmented analysis estimates block effects by **intra-block least squares from
  the checks** and reports the four classical Federer standard errors. For exact
  BLUPs on strongly unbalanced data, pair with a mixed-model package.
- Stability follows Eberhart & Russell (1966) and Gauch's AMMI.

## License

GPL-3. Author: Dr. Praveen Kumar B. K., Agriculture University, Jodhpur —
Genetics & Plant Breeding.
