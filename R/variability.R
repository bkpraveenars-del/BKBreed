#' Genetic variability parameters (GCV, PCV, heritability, genetic advance)
#'
#' Estimates the standard biometrical-genetics variability parameters from a
#' replicated (RBD) trial: genotypic, phenotypic and environmental variances and
#' coefficients of variation, broad-sense heritability, expected genetic advance
#' and genetic advance as per cent of mean, for one or several traits.
#'
#' Uses the conventional plant-breeding formulation (Singh & Chaudhary):
#' \eqn{\sigma^2_g = (M_g - M_e)/r}, \eqn{\sigma^2_e = M_e},
#' \eqn{\sigma^2_p = \sigma^2_g + \sigma^2_e}.
#'
#' @param data A data frame in long format.
#' @param traits Character vector of trait columns.
#' @param gen Character; genotype column.
#' @param rep Character; replication column.
#' @param k Selection differential (default 2.063 for 5\% selection intensity).
#' @return An object of class \code{bk_variability} (a data frame of parameters
#'   with metadata attributes).
#' @examples
#' v <- bk_variability(bk_data("rbd"),
#'        traits = c("grain_yield","plant_height","tillers","test_weight"),
#'        gen = "genotype", rep = "rep")
#' v
#' bk_plot(v)
#' @export
bk_variability <- function(data, traits, gen, rep, k = 2.063) {
  stopifnot(all(c(traits, gen, rep) %in% names(data)))
  rows <- lapply(traits, function(tr) {
    s  <- .gen_err_stats(data, tr, gen, rep)
    Vg <- max((s$Mg - s$Me) / s$r, 0)
    Ve <- s$Me
    Vp <- Vg + Ve
    mu <- s$mean
    h2 <- if (Vp > 0) Vg / Vp else NA_real_
    GA <- k * h2 * sqrt(Vp)
    data.frame(
      trait = tr, mean = mu,
      Vg = Vg, Ve = Ve, Vp = Vp,
      GCV = 100 * sqrt(Vg) / mu,
      PCV = 100 * sqrt(Vp) / mu,
      ECV = 100 * sqrt(Ve) / mu,
      h2_bs = 100 * h2,
      GA = GA, GAM = 100 * GA / mu,
      gen_p = s$p,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$GCV_class <- .band(out$GCV, c(10, 20))
  out$PCV_class <- .band(out$PCV, c(10, 20))
  out$h2_class  <- .band(out$h2_bs, c(30, 60))
  out$GAM_class <- .band(out$GAM, c(10, 20))
  class(out) <- c("bk_variability", "data.frame")
  attr(out, "gen") <- gen
  out
}

.band <- function(x, br) {
  ifelse(is.na(x), NA_character_,
  ifelse(x < br[1], "low",
  ifelse(x <= br[2], "moderate", "high")))
}

#' @export
print.bk_variability <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Genetic variability parameters\n")
  cat(strrep("=", 62), "\n\n", sep = "")
  d <- as.data.frame(x)
  num <- c("mean","Vg","Ve","Vp","GCV","PCV","ECV","h2_bs","GA","GAM")
  d[num] <- lapply(d[num], function(z) round(z, 2))
  print(d[c("trait","mean","GCV","PCV","h2_bs","GA","GAM",
            "GCV_class","h2_class","GAM_class")], row.names = FALSE)
  cat("\nGCV/PCV: low <10, moderate 10-20, high >20 (%)",
      "\nh2(bs) : low <30, moderate 30-60, high >60 (%)",
      "\nGAM    : low <10, moderate 10-20, high >20 (%)\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_variability <- function(x, ...) {
  d <- as.data.frame(x)
  long <- data.frame(
    trait = rep(d$trait, 2),
    kind  = rep(c("GCV", "PCV"), each = nrow(d)),
    value = c(d$GCV, d$PCV), stringsAsFactors = FALSE)
  long$trait <- factor(long$trait, levels = d$trait[order(d$PCV)])
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$trait,
                                     fill = .data$kind)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.65) +
    ggplot2::geom_text(
      data = d,
      ggplot2::aes(x = .data$PCV, y = .data$trait,
                   label = paste0("h2=", round(.data$h2_bs), "%")),
      inherit.aes = FALSE, hjust = -0.1, size = 2.8, colour = "#0B3D2E") +
    ggplot2::scale_fill_manual(values = c(GCV = "#1B7340", PCV = "#F2C14E"),
                               name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(
      title = "Genotypic and phenotypic variability",
      subtitle = "Bars: GCV vs PCV (%)   label: broad-sense heritability",
      x = "coefficient of variation (%)", y = NULL, caption = "BKBreed") +
    theme_bk()
}
