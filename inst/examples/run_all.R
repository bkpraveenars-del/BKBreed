# BKBreed — full demonstration. Run after installing the package.
# It exercises every analysis and saves each signature figure as PNG.
install.packages("BKBreed")
library(BKBreed)
library(ggplot2)

out <- file.path(tempdir(), "BKBreed_figures")
dir.create(out, showWarnings = FALSE)
sav <- function(p, f, w = 8, h = 5.2)
  ggsave(file.path(out, f), p, width = w, height = h, dpi = 300)

## 1. Randomised Block Design ------------------------------------------------
rbd <- bk_rbd(bk_data("rbd"), trait = "grain_yield",
              gen = "genotype", rep = "rep")
print(rbd); sav(bk_plot(rbd), "01_rbd_means.png")

## 2. Factorial RBD ----------------------------------------------------------
fr <- bk_frbd(bk_data("frbd"), trait = "grain_yield",
              factorA = "nitrogen", factorB = "variety", rep = "rep")
print(fr); sav(bk_plot(fr), "02_frbd_interaction.png")

## 3. Augmented Alpha-Lattice ------------------------------------------------
aug <- bk_augmented(bk_data("augmented"), trait = "grain_yield",
                    gen = "genotype", block = "block", rep = "rep",
                    checks = c("CHK-1","CHK-2","CHK-3","CHK-4"))
print(aug); sav(bk_plot(aug), "03_augmented_means.png", h = 7)

## 4. Genetic variability ----------------------------------------------------
traits <- c("grain_yield","plant_height","tillers","panicle_len","test_weight")
vb <- bk_variability(bk_data("rbd"), traits = traits,
                     gen = "genotype", rep = "rep")
print(vb); sav(bk_plot(vb), "04_variability.png")

## 5. Correlation ------------------------------------------------------------
cr <- bk_correlation(bk_data("rbd"), traits = traits,
                     gen = "genotype", rep = "rep")
print(cr); sav(bk_plot(cr, type = "genotypic"), "05_correlation.png", w = 6.5, h = 6)

## 6. Path analysis ----------------------------------------------------------
pa <- bk_path(bk_data("rbd"),
              traits = c("plant_height","tillers","panicle_len","test_weight"),
              dependent = "grain_yield", gen = "genotype", rep = "rep")
print(pa); sav(bk_plot(pa), "06_path.png", w = 6.5, h = 6)

## 7. D-square divergence ----------------------------------------------------
dv <- bk_diversity(bk_data("rbd"), traits = traits,
                   gen = "genotype", rep = "rep")
print(dv); sav(bk_plot(dv), "07_diversity.png")

## 8. GxE stability (MLT) ----------------------------------------------------
st <- bk_stability(bk_data("mlt"), trait = "grain_yield",
                   gen = "genotype", env = "environment", rep = "rep")
print(st)
sav(bk_plot(st, type = "ammi2"), "08_ammi2_biplot.png", w = 7.5, h = 6)
sav(bk_plot(st, type = "er"),    "09_eberhart_russell.png")

cat("\nAll analyses complete. Figures written to '", out, "/'\n", sep = "")
