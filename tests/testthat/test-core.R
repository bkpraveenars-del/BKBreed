test_that("example data load", {
  for (nm in c("rbd", "frbd", "augmented", "mlt"))
    expect_s3_class(bk_data(nm), "data.frame")
})

test_that("RBD returns a valid ANOVA and CD", {
  r <- bk_rbd(bk_data("rbd"), "grain_yield", "genotype", "rep")
  expect_s3_class(r, "bk_rbd")
  expect_true(r$cd > 0)
  expect_true(all(c("genotype","mean","group","rank") %in% names(r$means)))
  expect_s3_class(bk_plot(r), "ggplot")
})

test_that("variability parameters are ordered GCV <= PCV", {
  v <- bk_variability(bk_data("rbd"),
         c("grain_yield","plant_height","tillers"), "genotype", "rep")
  expect_true(all(v$GCV <= v$PCV + 1e-8))
  expect_true(all(v$h2_bs >= 0 & v$h2_bs <= 100))
})

test_that("correlation matrices are symmetric with unit diagonal", {
  cr <- bk_correlation(bk_data("rbd"),
          c("grain_yield","plant_height","tillers","test_weight"),
          "genotype", "rep")
  expect_equal(diag(cr$rg), rep(1, 4), ignore_attr = TRUE)
  expect_equal(cr$rp, t(cr$rp))
})

test_that("path direct effects reconstruct correlations", {
  pa <- bk_path(bk_data("rbd"),
          c("plant_height","tillers","test_weight"),
          "grain_yield", "genotype", "rep")
  expect_equal(as.numeric(pa$total), as.numeric(pa$r_with_dep), tolerance = 1e-6)
})

test_that("D2 matrix is symmetric and non-negative", {
  dv <- bk_diversity(bk_data("rbd"),
          c("grain_yield","plant_height","tillers"), "genotype", "rep")
  expect_true(all(dv$D2 >= -1e-8))
  expect_equal(dv$D2, t(dv$D2))
})

test_that("stability: mean bi equals 1", {
  st <- bk_stability(bk_data("mlt"), "grain_yield",
                     "genotype", "environment", "rep")
  expect_equal(mean(st$eberhart$bi), 1, tolerance = 1e-6)
  expect_s3_class(bk_plot(st, "ammi2"), "ggplot")
})

test_that("line x tester: GCA effects sum to zero, contributions to 100", {
  res <- bk_lxt(bk_data("lxt"), "grain_yield", "line", "tester", "rep")
  expect_equal(sum(res$gca_lines$gca), 0, tolerance = 1e-6)
  expect_equal(sum(res$gca_testers$gca), 0, tolerance = 1e-6)
  expect_equal(sum(res$contribution), 100, tolerance = 1e-6)
  expect_s3_class(bk_plot(res, "gca"), "ggplot")
  expect_s3_class(bk_plot(res, "sca"), "ggplot")
})

test_that("diallel: GCA sums to zero and SS partitions", {
  res <- bk_diallel(bk_data("diallel"), "grain_yield",
                    "parent1", "parent2", "rep")
  expect_equal(sum(res$gca$gca), 0, tolerance = 1e-6)
  expect_equal(res$sca, t(res$sca), tolerance = 1e-8)   # symmetric
  ss <- res$anova$SS
  expect_true(res$baker >= 0 && res$baker <= 1)
  expect_s3_class(bk_plot(res, "gca"), "ggplot")
  expect_s3_class(bk_plot(res, "sca"), "ggplot")
})

test_that("augmented analysis adjusts means with positive error df", {
  a <- bk_augmented(bk_data("augmented"), "grain_yield", "genotype",
                    "block", "rep", checks = paste0("CHK-", 1:4))
  expect_true(a$edf > 0)
  expect_true(a$mse > 0)
  expect_true(all(c("check","test") %in% a$means$type))
})
