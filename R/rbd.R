#' Randomised Block Design (RBD) analysis
#'
#' One-call ANOVA for a randomised complete block design, returning the ANOVA
#' table (with SS, MS, F, p and significance stars), adjusted genotype means,
#' standard errors, critical difference (CD/LSD), coefficient of variation and
#' compact-letter groupings.
#'
#' @param data A data frame in long format.
#' @param trait Character; name of the response (yield, height, ...).
#' @param gen Character; name of the treatment/genotype column.
#' @param rep Character; name of the replication/block column.
#' @param alpha Significance level for the critical difference (default 0.05).
#' @return An object of class \code{bk_rbd}: a list with \code{anova},
#'   \code{means}, \code{cv}, \code{sem}, \code{sed}, \code{cd} and metadata.
#' @examples
#' res <- bk_rbd(bk_data("rbd"), trait = "grain_yield",
#'               gen = "genotype", rep = "rep")
#' res
#' bk_plot(res)
#' @export
bk_rbd <- function(data, trait, gen, rep, alpha = 0.05) {
  stopifnot(all(c(trait, gen, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(trait, gen, rep)]), ]
  d[[gen]] <- factor(d[[gen]])
  d[[rep]] <- factor(d[[rep]])
  r <- nlevels(d[[rep]])

  f   <- stats::as.formula(paste0("`", trait, "` ~ ", rep, " + ", gen))
  fit <- stats::aov(f, data = d)
  av  <- .tidy_aov(fit)

  Me   <- av$MS[nrow(av)]
  edf  <- av$Df[nrow(av)]
  gm   <- mean(d[[trait]], na.rm = TRUE)
  sem  <- sqrt(Me / r)
  sed  <- sqrt(2 * Me / r)
  cd   <- stats::qt(1 - alpha / 2, edf) * sed
  cv   <- 100 * sqrt(Me) / gm

  mns  <- tapply(d[[trait]], d[[gen]], mean, na.rm = TRUE)
  grp  <- .cld(mns, cd)
  means <- data.frame(
    genotype = names(mns),
    mean = as.numeric(mns),
    group = grp[names(mns)],
    row.names = NULL, stringsAsFactors = FALSE
  )
  means <- means[order(-means$mean), ]
  means$rank <- seq_len(nrow(means))
  rownames(means) <- NULL

  structure(list(
    anova = av, means = means, grand_mean = gm,
    cv = cv, sem = sem, sed = sed, cd = cd, alpha = alpha,
    trait = trait, gen = gen, rep = rep, r = r, n_gen = nlevels(d[[gen]])
  ), class = "bk_rbd")
}

#' @export
print.bk_rbd <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Randomised Block Design  |  trait: ", x$trait, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  a <- x$anova
  a$SS <- round(a$SS, 3); a$MS <- round(a$MS, 3)
  a$F  <- ifelse(is.na(a$F), "", formatC(a$F, format = "f", digits = 2))
  a$p  <- ifelse(is.na(a$p), "", formatC(a$p, format = "f", digits = 4))
  print(a, row.names = FALSE)
  cat("\nGrand mean : ", round(x$grand_mean, 3),
      "   CV% : ", round(x$cv, 2),
      "\nSE(m)     : ", round(x$sem, 3),
      "   SE(d): ", round(x$sed, 3),
      "   CD(", x$alpha, "): ", round(x$cd, 3), "\n", sep = "")
  cat("\nGenotype means (descending; shared letters = not significant):\n")
  m <- x$means; m$mean <- round(m$mean, 3)
  print(m[c("rank", "genotype", "mean", "group")], row.names = FALSE)
  cat("\nSignif: *** p<0.001  ** p<0.01  * p<0.05  . p<0.10  ns >=0.10\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_rbd <- function(x, ...) {
  m <- x$means
  m$genotype <- factor(m$genotype, levels = rev(m$genotype))
  ggplot2::ggplot(m, ggplot2::aes(x = mean, y = .data$genotype,
                                  fill = .data$mean)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = mean - x$sem, xmax = mean + x$sem),
      width = 0.28, colour = "#333333", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = .data$group),
                       hjust = -0.35, size = 3.1, colour = "#0B3D2E",
                       fontface = "bold") +
    ggplot2::scale_fill_gradientn(colours = bk_palette("field"),
                                  name = x$trait) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      title = paste0("RBD genotype means — ", x$trait),
      subtitle = paste0("CD(", x$alpha, ") = ", round(x$cd, 2),
                        "   CV% = ", round(x$cv, 1),
                        "   letters share = non-significant"),
      x = x$trait, y = NULL,
      caption = "BKBreed") +
    theme_bk() +
    ggplot2::theme(legend.position = "none")
}
