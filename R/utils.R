# ---- internal utilities (not exported) --------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# significance stars from a p-value
.stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.10,  ".", "ns")))))
}

# tidy one-way/two-way aov into a clean ANOVA data.frame
.tidy_aov <- function(fit) {
  a <- stats::anova(fit)
  data.frame(
    Source = trimws(rownames(a)),
    Df     = a[["Df"]],
    SS     = a[["Sum Sq"]],
    MS     = a[["Mean Sq"]],
    F      = a[["F value"]],
    p      = a[["Pr(>F)"]],
    Sig    = .stars(a[["Pr(>F)"]]),
    row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE
  )
}

# compact letter display for a set of means compared with a single LSD/CD.
# Valid for balanced designs where every pairwise comparison shares one CD.
.cld <- function(means, cd) {
  nm  <- names(means)
  ord <- order(means, decreasing = TRUE)
  m   <- means[ord]
  n   <- length(m)
  if (n == 1) { out <- setNames("a", nm); return(out) }
  # rightmost index still within cd of each start (means sorted descending)
  right <- integer(n)
  for (i in seq_len(n)) {
    j <- i
    while (j < n && (m[i] - m[j + 1]) <= cd) j <- j + 1
    right[i] <- j
  }
  keep <- rep(TRUE, n)
  for (i in 2:n) if (right[i] <= right[i - 1]) keep[i] <- FALSE
  starts <- which(keep)
  labs   <- character(n)
  for (c in seq_along(starts)) {
    s <- starts[c]; e <- right[s]
    lett <- letters[((c - 1) %% 26) + 1]
    labs[s:e] <- paste0(labs[s:e], lett)
  }
  out <- character(n)
  out[ord] <- labs
  names(out) <- nm
  out
}

# genotype/error mean squares + cross products for a two-way (rep + gen) model
# returns a list used by variability, correlation and path routines
.gen_err_stats <- function(df, trait, gen, rep) {
  df[[gen]] <- factor(df[[gen]])
  df[[rep]] <- factor(df[[rep]])
  r  <- nlevels(df[[rep]])
  f  <- stats::as.formula(paste0("`", trait, "` ~ ", rep, " + ", gen))
  a  <- stats::anova(stats::aov(f, data = df))
  gi <- match(gen, trimws(rownames(a)))
  ei <- nrow(a)
  list(
    Mg = a[["Mean Sq"]][gi], gdf = a[["Df"]][gi],
    Me = a[["Mean Sq"]][ei], edf = a[["Df"]][ei],
    r  = r, mean = mean(df[[trait]], na.rm = TRUE),
    p  = a[["Pr(>F)"]][gi]
  )
}

# mean cross product for a source, via the SS(x+y) identity
.mcp <- function(df, tx, ty, gen, rep, source = c("gen", "err")) {
  source <- match.arg(source)
  df$.z <- df[[tx]] + df[[ty]]
  gx <- .gen_err_stats(df, tx, gen, rep)
  gy <- .gen_err_stats(df, ty, gen, rep)
  gz <- .gen_err_stats(df, ".z", gen, rep)
  if (source == "gen") {
    ssx <- gx$Mg * gx$gdf; ssy <- gy$Mg * gy$gdf; ssz <- gz$Mg * gz$gdf
    (0.5 * (ssz - ssx - ssy)) / gx$gdf
  } else {
    ssx <- gx$Me * gx$edf; ssy <- gy$Me * gy$edf; ssz <- gz$Me * gz$edf
    (0.5 * (ssz - ssx - ssy)) / gx$edf
  }
}

#' Load a bundled BKBreed example dataset
#'
#' Convenience loader for the demonstration datasets shipped with the package
#' (stored as CSV in \code{inst/extdata}).
#'
#' @param name One of \code{"rbd"} (pearl-millet multi-trait RBD),
#'   \code{"frbd"} (nitrogen x variety factorial RBD),
#'   \code{"augmented"} (augmented alpha-lattice with 4 checks + 20 test entries),
#'   \code{"mlt"} (10 genotypes x 5 environments multi-location trial),
#'   \code{"lxt"} (5 lines x 3 testers line x tester crosses),
#'   \item \code{"diallel"} (6-parent half diallel: parents + F1s, 3 reps).
#' @return A \code{data.frame}.
#' @examples
#' head(bk_data("rbd"))
#' @export
bk_data <- function(name = c("rbd", "frbd", "augmented", "mlt", "lxt",
                             "diallel")) {
  name <- match.arg(name)
  file <- switch(name,
    rbd       = "pearlmillet_rbd.csv",
    frbd      = "factorial_rbd.csv",
    augmented = "augmented_alpha.csv",
    mlt       = "mlt_stability.csv",
    lxt       = "linextester.csv",
    diallel   = "diallel.csv")
  path <- system.file("extdata", file, package = "BKBreed")
  if (path == "") stop("Example data '", file, "' not found. Is BKBreed installed?")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Draw the signature figure for a BKBreed result
#'
#' A single generic entry point: \code{bk_plot(x)} returns the publication-ready
#' \code{ggplot2} figure appropriate to whichever analysis produced \code{x}.
#'
#' @param x A BKBreed result object.
#' @param ... Passed to the specific method.
#' @return A \code{ggplot2} object (or a \code{patchwork} composite).
#' @export
bk_plot <- function(x, ...) UseMethod("bk_plot")
