#' Genotypic and phenotypic correlation
#'
#' Estimates genotypic and phenotypic correlation coefficients among a set of
#' traits from a replicated (RBD) trial, using analysis of covariance to obtain
#' genotypic and error mean cross-products.
#'
#' @param data A data frame in long format.
#' @param traits Character vector of trait columns (>= 2).
#' @param gen Character; genotype column.
#' @param rep Character; replication column.
#' @return An object of class \code{bk_correlation} holding the genotypic
#'   (\code{rg}) and phenotypic (\code{rp}) correlation matrices.
#' @examples
#' cr <- bk_correlation(bk_data("rbd"),
#'         traits = c("grain_yield","plant_height","tillers",
#'                    "panicle_len","test_weight"),
#'         gen = "genotype", rep = "rep")
#' cr
#' bk_plot(cr)
#' @export
bk_correlation <- function(data, traits, gen, rep) {
  stopifnot(length(traits) >= 2, all(c(traits, gen, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(traits, gen, rep)]), ]
  p <- length(traits)

  # per-trait genotypic and error variances
  st <- lapply(traits, function(tr) .gen_err_stats(d, tr, gen, rep))
  names(st) <- traits
  Vg <- vapply(st, function(s) max((s$Mg - s$Me) / s$r, 0), numeric(1))
  Ve <- vapply(st, function(s) s$Me, numeric(1))
  Vp <- Vg + Ve
  r  <- st[[1]]$r

  rg <- diag(p); rp <- diag(p)
  dimnames(rg) <- dimnames(rp) <- list(traits, traits)
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    mcg <- .mcp(d, traits[i], traits[j], gen, rep, "gen")
    mce <- .mcp(d, traits[i], traits[j], gen, rep, "err")
    covg <- (mcg - mce) / r
    covp <- covg + mce
    dg <- sqrt(Vg[i] * Vg[j]); dp <- sqrt(Vp[i] * Vp[j])
    rg[i, j] <- rg[j, i] <- if (dg > 0) max(min(covg / dg, 1), -1) else NA
    rp[i, j] <- rp[j, i] <- if (dp > 0) max(min(covp / dp, 1), -1) else NA
  }
  structure(list(rg = rg, rp = rp, traits = traits, Vg = Vg, Vp = Vp),
            class = "bk_correlation")
}

#' @export
print.bk_correlation <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Correlation coefficients\n")
  cat(strrep("=", 62), "\n\n", sep = "")
  cat("Genotypic correlation (rg):\n")
  print(round(x$rg, 3))
  cat("\nPhenotypic correlation (rp):\n")
  print(round(x$rp, 3))
  cat("\n")
  invisible(x)
}

#' @export
bk_plot.bk_correlation <- function(x, type = c("genotypic", "phenotypic"), ...) {
  type <- match.arg(type)
  m <- if (type == "genotypic") x$rg else x$rp
  tr <- x$traits
  long <- expand.grid(row = tr, col = tr, stringsAsFactors = FALSE)
  long$value <- as.vector(m)
  long$row <- factor(long$row, levels = rev(tr))
  long$col <- factor(long$col, levels = tr)
  # show numbers only on lower triangle + diagonal for clarity
  idx <- as.integer(long$col) <= (length(tr) - as.integer(long$row) + 1)
  long$lab <- ifelse(idx, formatC(long$value, format = "f", digits = 2), "")
  ggplot2::ggplot(long, ggplot2::aes(.data$col, .data$row, fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 3,
                       colour = "#1A1A1A") +
    ggplot2::scale_fill_gradientn(
      colours = bk_palette("spectrum"), limits = c(-1, 1),
      name = "r") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(toupper(substring(type, 1, 1)), substring(type, 2),
                     " correlation matrix"),
      subtitle = "blue = negative, red = positive association",
      x = NULL, y = NULL, caption = "BKBreed") +
    theme_bk() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   panel.grid = ggplot2::element_blank())
}
