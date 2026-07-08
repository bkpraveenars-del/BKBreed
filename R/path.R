#' Path-coefficient analysis
#'
#' Partitions the correlation of each causal (independent) trait with a target
#' (dependent) trait into direct and indirect effects via Wright's path
#' analysis, and reports the residual effect.
#'
#' @param data A data frame in long format.
#' @param traits Character vector of causal traits (independent variables).
#' @param dependent Character; the target/effect trait (e.g. grain yield).
#' @param gen Character; genotype column.
#' @param rep Character; replication column.
#' @param type Correlation basis: \code{"genotypic"} (default) or
#'   \code{"phenotypic"}.
#' @return An object of class \code{bk_path}.
#' @examples
#' pa <- bk_path(bk_data("rbd"),
#'         traits = c("plant_height","tillers","panicle_len","test_weight"),
#'         dependent = "grain_yield", gen = "genotype", rep = "rep")
#' pa
#' bk_plot(pa)
#' @export
bk_path <- function(data, traits, dependent, gen, rep,
                    type = c("genotypic", "phenotypic")) {
  type <- match.arg(type)
  stopifnot(!(dependent %in% traits))
  cr  <- bk_correlation(data, c(traits, dependent), gen, rep)
  M   <- if (type == "genotypic") cr$rg else cr$rp
  R   <- M[traits, traits, drop = FALSE]
  rvec <- M[traits, dependent]
  if (any(is.na(R)) || any(is.na(rvec)))
    stop("Correlation matrix contains NA; check for zero-variance traits.")

  direct <- as.numeric(solve(R, rvec))
  names(direct) <- traits
  eff <- R * matrix(direct, nrow = length(traits), ncol = length(traits),
                    byrow = TRUE)          # eff[j,k] = P[k]*r[j,k]; diagonal=direct
  dimnames(eff) <- list(traits, traits)
  total <- rowSums(eff)                    # equals rvec (reconstruction check)
  resid_sq <- 1 - sum(direct * rvec)
  residual <- if (resid_sq >= 0) sqrt(resid_sq) else NA_real_
  if (is.na(residual))
    warning("Residual effect is imaginary (R^2 > 1); model may be over-determined.")

  structure(list(
    effects = eff, direct = direct, total = total, r_with_dep = rvec,
    residual = residual, type = type, traits = traits,
    dependent = dependent), class = "bk_path")
}

#' @export
print.bk_path <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Path analysis (", x$type, ")  ->  ", x$dependent, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  cat("Direct effects on the diagonal; indirect effects off-diagonal:\n")
  m <- round(x$effects, 3)
  print(m)
  cat("\nDirect effects:\n")
  print(round(x$direct, 3))
  cat("\nCorrelation with ", x$dependent, " (reconstructed):\n", sep = "")
  print(round(x$total, 3))
  cat("\nResidual effect:", round(x$residual, 3), "\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_path <- function(x, ...) {
  tr <- x$traits
  long <- expand.grid(via = tr, cause = tr, stringsAsFactors = FALSE)
  long$value <- as.vector(x$effects)   # column-major: value = eff[via, cause]
  long$type  <- ifelse(long$via == long$cause, "direct", "indirect")
  long$via   <- factor(long$via, levels = rev(tr))
  long$cause <- factor(long$cause, levels = tr)
  ggplot2::ggplot(long, ggplot2::aes(.data$cause, .data$via,
                                     fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_tile(data = long[long$type == "direct", , drop = FALSE],
                       colour = "#1A1A1A", linewidth = 1.1, fill = NA) +
    ggplot2::geom_text(ggplot2::aes(label = formatC(.data$value, format = "f",
                                                    digits = 2)),
                       size = 3, colour = "#1A1A1A") +
    ggplot2::scale_fill_gradientn(colours = bk_palette("spectrum"),
                                  name = "effect") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("Path analysis on ", x$dependent, " (", x$type, ")"),
      subtitle = paste0("boxed diagonal = direct effect   |   residual = ",
                        round(x$residual, 3)),
      x = "causal trait (direct effect)",
      y = "acting via", caption = "BKBreed") +
    theme_bk() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}
