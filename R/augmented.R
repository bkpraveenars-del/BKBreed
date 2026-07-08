#' Augmented Alpha-Lattice (and Augmented RCBD) analysis
#'
#' Analysis of augmented designs in which a set of replicated \emph{checks} is
#' used to estimate block effects, and a large number of unreplicated
#' \emph{test entries} is adjusted for those block effects. Block effects are
#' estimated by intra-block least squares from the checks only; every genotype
#' mean is then adjusted and the classical Federer standard errors are reported
#' for each type of comparison.
#'
#' If \code{rep} is \code{NULL} the layout is treated as an augmented RCBD with
#' \code{block} as complete blocks. If \code{rep} is supplied, blocks are treated
#' as incomplete blocks nested within replications (an alpha-lattice layout) and
#' the blocking identifier becomes \code{rep:block}.
#'
#' @param data A data frame in long format.
#' @param trait Character; response column.
#' @param gen Character; genotype column (checks and test entries together).
#' @param block Character; (incomplete) block column.
#' @param rep Character or \code{NULL}; replication column for alpha-lattice.
#' @param checks Character vector of check genotype names. If \code{NULL}, any
#'   genotype occurring more than once is treated as a check.
#' @param alpha Significance level for critical differences (default 0.05).
#' @return An object of class \code{bk_augmented}.
#' @examples
#' d <- bk_data("augmented")
#' res <- bk_augmented(d, trait = "grain_yield", gen = "genotype",
#'                     block = "block", rep = "rep",
#'                     checks = c("CHK-1","CHK-2","CHK-3","CHK-4"))
#' res
#' bk_plot(res)
#' @export
bk_augmented <- function(data, trait, gen, block, rep = NULL,
                         checks = NULL, alpha = 0.05) {
  stopifnot(all(c(trait, gen, block) %in% names(data)))
  need <- c(trait, gen, block, if (!is.null(rep)) rep)
  d <- data[stats::complete.cases(data[need]), ]
  d[[gen]] <- as.character(d[[gen]])

  # blocking identifier
  blockid <- if (is.null(rep)) as.character(d[[block]]) else
             paste(d[[rep]], d[[block]], sep = ":")
  d$.bl <- factor(blockid)

  # infer checks if not given
  if (is.null(checks)) {
    tab <- table(d[[gen]])
    checks <- names(tab)[tab > 1]
  }
  if (length(checks) < 2)
    stop("At least two replicated checks are required for an augmented analysis.")
  is_chk <- d[[gen]] %in% checks

  # ---- estimate block effects from checks by intra-block least squares ----
  chk <- d[is_chk, ]
  chk$.g  <- droplevels(factor(chk[[gen]]))
  chk$.bl <- droplevels(chk$.bl)
  if (nlevels(chk$.bl) < 2)
    stop("Checks must appear in at least two blocks to estimate block effects.")
  fit <- stats::lm(stats::as.formula(paste0("`", trait, "` ~ .g + .bl")),
                   data = chk)
  edf <- fit$df.residual
  if (edf < 1)
    stop("No residual degrees of freedom: not enough replicated check data.")
  MSE <- sum(stats::residuals(fit)^2) / edf

  # centred block effects via prediction at a fixed reference genotype
  bls <- levels(chk$.bl)
  ref <- data.frame(.g = factor(levels(chk$.g)[1], levels = levels(chk$.g)),
                    .bl = factor(bls, levels = levels(chk$.bl)))
  be  <- stats::predict(fit, ref)
  be  <- be - mean(be)
  names(be) <- bls

  # ---- adjust every observation, then average per genotype ----
  adj <- be[as.character(d$.bl)]
  n_missing <- sum(is.na(adj))
  if (n_missing > 0)
    warning(n_missing, " plot(s) sit in a block with no check; ",
            "their block adjustment is set to 0.")
  adj[is.na(adj)] <- 0
  d$.yadj <- d[[trait]] - adj

  amean <- tapply(d$.yadj, d[[gen]], mean)
  cnt   <- tapply(d[[gen]], d[[gen]], length)
  means <- data.frame(
    genotype = names(amean),
    adj_mean = as.numeric(amean),
    n = as.integer(cnt[names(amean)]),
    type = ifelse(names(amean) %in% checks, "check", "test"),
    row.names = NULL, stringsAsFactors = FALSE)
  means <- means[order(-means$adj_mean), ]
  means$rank <- seq_len(nrow(means))
  rownames(means) <- NULL

  b  <- nlevels(chk$.bl)                # number of blocks with checks
  cc <- length(checks)                  # number of checks
  rc <- stats::median(cnt[checks])      # replication of checks
  tval <- stats::qt(1 - alpha / 2, edf)
  gm  <- mean(d$.yadj, na.rm = TRUE)

  se <- list(
    check_check    = sqrt(2 * MSE / rc),
    test_test_same = sqrt(2 * MSE),
    test_test_diff = sqrt(2 * MSE * (b + 1) / b),
    test_check     = sqrt(MSE * (1 + 1 / b) * (1 + 1 / cc))
  )
  cd <- lapply(se, function(s) tval * s)

  structure(list(
    means = means, block_effects = be, mse = MSE, edf = edf,
    cv = 100 * sqrt(MSE) / gm, grand_mean = gm,
    se = se, cd = cd, alpha = alpha,
    n_block = b, n_check = cc, rc = rc,
    checks = checks, trait = trait, gen = gen,
    design = if (is.null(rep)) "Augmented RCBD" else "Augmented Alpha-Lattice"
  ), class = "bk_augmented")
}

#' @export
print.bk_augmented <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | ", x$design, "  |  trait: ", x$trait, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  cat("Blocks:", x$n_block, "  Checks:", x$n_check,
      "  Test entries:", sum(x$means$type == "test"),
      "\nError MS:", round(x$mse, 3), " (df=", x$edf, ")",
      "   CV%:", round(x$cv, 2), "\n\n", sep = "")
  cat("Critical differences (CD, alpha=", x$alpha, "):\n", sep = "")
  cat("  two checks               :", round(x$cd$check_check, 3), "\n")
  cat("  two tests (same block)   :", round(x$cd$test_test_same, 3), "\n")
  cat("  two tests (diff. blocks) :", round(x$cd$test_test_diff, 3), "\n")
  cat("  test vs check            :", round(x$cd$test_check, 3), "\n\n")
  cat("Adjusted means (top 15 by rank):\n")
  m <- x$means; m$adj_mean <- round(m$adj_mean, 3)
  print(utils::head(m[c("rank", "genotype", "adj_mean", "type")], 15),
        row.names = FALSE)
  cat("\n")
  invisible(x)
}

#' @export
bk_plot.bk_augmented <- function(x, ...) {
  m <- x$means
  m$genotype <- factor(m$genotype, levels = rev(m$genotype))
  chk_mean <- mean(m$adj_mean[m$type == "check"])
  ggplot2::ggplot(m, ggplot2::aes(x = .data$adj_mean, y = .data$genotype,
                                  colour = .data$type)) +
    ggplot2::geom_vline(xintercept = chk_mean, linetype = "dashed",
                        colour = "#B5341C", linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = x$grand_mean, xend = .data$adj_mean,
                                       yend = .data$genotype),
                          colour = "#CFCFC4", linewidth = 0.5) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = c(check = "#B5341C", test = "#1B7340"),
                                 name = NULL) +
    ggplot2::annotate("text", x = chk_mean, y = nrow(m),
                      label = "check mean", hjust = -0.05, vjust = 1,
                      colour = "#B5341C", size = 3) +
    ggplot2::labs(
      title = paste0(x$design, " - adjusted means"),
      subtitle = paste0(x$trait, "   CV% = ", round(x$cv, 1),
                        "   error df = ", x$edf),
      x = paste0("block-adjusted ", x$trait), y = NULL,
      caption = "BKBreed  |  red dashed = mean of checks") +
    theme_bk()
}
