#' Two-factor Factorial RBD analysis
#'
#' ANOVA for a two-factor factorial laid out in a randomised block design,
#' with main-effect and interaction-cell means, factor-specific standard errors
#' and critical differences, and coefficient of variation.
#'
#' @param data A data frame in long format.
#' @param trait Character; response column.
#' @param factorA,factorB Character; the two treatment factor columns.
#' @param rep Character; replication/block column.
#' @param alpha Significance level for the critical differences (default 0.05).
#' @return An object of class \code{bk_frbd}.
#' @examples
#' res <- bk_frbd(bk_data("frbd"), trait = "grain_yield",
#'                factorA = "nitrogen", factorB = "variety", rep = "rep")
#' res
#' bk_plot(res)
#' @export
bk_frbd <- function(data, trait, factorA, factorB, rep, alpha = 0.05) {
  stopifnot(all(c(trait, factorA, factorB, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(trait, factorA, factorB, rep)]), ]
  d[[factorA]] <- factor(d[[factorA]])
  d[[factorB]] <- factor(d[[factorB]])
  d[[rep]]     <- factor(d[[rep]])
  a <- nlevels(d[[factorA]]); b <- nlevels(d[[factorB]]); r <- nlevels(d[[rep]])

  f <- stats::as.formula(paste0("`", trait, "` ~ ", rep, " + ",
                                factorA, " * ", factorB))
  fit <- stats::aov(f, data = d)
  av  <- .tidy_aov(fit)

  Me  <- av$MS[nrow(av)]; edf <- av$Df[nrow(av)]
  gm  <- mean(d[[trait]], na.rm = TRUE)
  tval <- stats::qt(1 - alpha / 2, edf)
  cv  <- 100 * sqrt(Me) / gm

  seA <- sqrt(2 * Me / (r * b)); cdA <- tval * seA
  seB <- sqrt(2 * Me / (r * a)); cdB <- tval * seB
  seAB <- sqrt(2 * Me / r);      cdAB <- tval * seAB

  mA <- tapply(d[[trait]], d[[factorA]], mean, na.rm = TRUE)
  mB <- tapply(d[[trait]], d[[factorB]], mean, na.rm = TRUE)
  meansA <- data.frame(level = names(mA), mean = as.numeric(mA),
                       group = .cld(mA, cdA)[names(mA)],
                       row.names = NULL, stringsAsFactors = FALSE)
  meansB <- data.frame(level = names(mB), mean = as.numeric(mB),
                       group = .cld(mB, cdB)[names(mB)],
                       row.names = NULL, stringsAsFactors = FALSE)
  cell <- stats::aggregate(stats::as.formula(paste0("`", trait, "` ~ ",
                    factorA, " + ", factorB)), data = d, FUN = mean)
  names(cell) <- c(factorA, factorB, "mean")

  structure(list(
    anova = av, meansA = meansA, meansB = meansB, cell = cell,
    grand_mean = gm, cv = cv,
    se = list(A = seA, B = seB, AB = seAB),
    cd = list(A = cdA, B = cdB, AB = cdAB),
    alpha = alpha, trait = trait, factorA = factorA, factorB = factorB,
    a = a, b = b, r = r
  ), class = "bk_frbd")
}

#' @export
print.bk_frbd <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Factorial RBD (", x$factorA, " x ", x$factorB,
      ")  |  trait: ", x$trait, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  a <- x$anova
  a$SS <- round(a$SS, 3); a$MS <- round(a$MS, 3)
  a$F  <- ifelse(is.na(a$F), "", formatC(a$F, format = "f", digits = 2))
  a$p  <- ifelse(is.na(a$p), "", formatC(a$p, format = "f", digits = 4))
  print(a, row.names = FALSE)
  cat("\nGrand mean:", round(x$grand_mean, 3), "  CV%:", round(x$cv, 2), "\n")
  cat("CD(", x$alpha, ")  ", x$factorA, ": ", round(x$cd$A, 3),
      "  |  ", x$factorB, ": ", round(x$cd$B, 3),
      "  |  interaction: ", round(x$cd$AB, 3), "\n", sep = "")
  cat("\n", x$factorA, " means:\n", sep = "")
  mAd <- x$meansA; mAd$mean <- round(mAd$mean, 3); print(mAd, row.names = FALSE)
  cat("\n", x$factorB, " means:\n", sep = "")
  mBd <- x$meansB; mBd$mean <- round(mBd$mean, 3); print(mBd, row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
bk_plot.bk_frbd <- function(x, ...) {
  cell <- x$cell
  ggplot2::ggplot(cell, ggplot2::aes(x = .data[[x$factorA]], y = .data$mean,
                                     colour = .data[[x$factorB]],
                                     group = .data[[x$factorB]])) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3.2) +
    scale_colour_bk("sunrise", name = x$factorB) +
    ggplot2::labs(
      title = paste0("Factorial RBD interaction — ", x$trait),
      subtitle = paste0(x$factorA, " x ", x$factorB,
                        "   CD(int) = ", round(x$cd$AB, 2),
                        "   CV% = ", round(x$cv, 1)),
      x = x$factorA, y = paste0("mean ", x$trait),
      caption = "BKBreed") +
    theme_bk()
}
