#' BKBreed: Colourful Biometrical Analysis for Plant Breeding and Genetics
#'
#' A compact, colour-first toolkit for analysing plant-breeding and genetics
#' field experiments. Each function returns a tidy result object and a
#' publication-ready \code{ggplot2} figure drawn with a bespoke colour system.
#'
#' @section Experimental designs:
#' \itemize{
#'   \item \code{\link{bk_rbd}} - Randomised Block Design
#'   \item \code{\link{bk_frbd}} - two-factor Factorial RBD
#'   \item \code{\link{bk_augmented}} - Augmented Alpha-Lattice / Augmented RCBD
#' }
#'
#' @section Biometrical genetics:
#' \itemize{
#'   \item \code{\link{bk_variability}} - GCV, PCV, heritability, genetic advance
#'   \item \code{\link{bk_correlation}} - genotypic & phenotypic correlation
#'   \item \code{\link{bk_path}} - path-coefficient analysis
#'   \item \code{\link{bk_lxt}} - Line x Tester combining ability (GCA/SCA)
#'   \item \code{\link{bk_diallel}} - Griffing diallel (Method 2, Model I)
#'   \item \code{\link{bk_diversity}} - Mahalanobis D-square + Tocher / Ward
#'   \item \code{\link{bk_stability}} - Eberhart-Russell & AMMI (MLT)
#' }
#'
#' @section Figures:
#' Call \code{\link{bk_plot}} on any result object to obtain its signature
#' figure. Theme and palettes: \code{\link{theme_bk}}, \code{\link{bk_palette}}.
#'
#' @section Data:
#' Four example datasets ship with the package; load them with
#' \code{\link{bk_data}}.
#'
#' @keywords internal
"_PACKAGE"
