#' Genotype x Environment stability analysis (Eberhart-Russell + AMMI)
#'
#' Analyses multi-location / multi-environment trial (MLT) data. Returns the
#' combined ANOVA, the Eberhart & Russell (1966) stability parameters
#' (regression coefficient \eqn{b_i} and deviation from regression
#' \eqn{S^2_{di}}), and an AMMI (Additive Main effects and Multiplicative
#' Interaction) decomposition with interaction principal-component scores for
#' biplots.
#'
#' @param data A data frame in long format.
#' @param trait Character; response column.
#' @param gen Character; genotype column.
#' @param env Character; environment / location column.
#' @param rep Character; replication column (within environment).
#' @return An object of class \code{bk_stability}.
#' @examples
#' st <- bk_stability(bk_data("mlt"), trait = "grain_yield",
#'                    gen = "genotype", env = "environment", rep = "rep")
#' st
#' bk_plot(st)                 # AMMI-2 biplot
#' bk_plot(st, type = "er")    # Eberhart-Russell
#' @export
bk_stability <- function(data, trait, gen, env, rep) {
  stopifnot(all(c(trait, gen, env, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(trait, gen, env, rep)]), ]
  d[[gen]] <- factor(d[[gen]]); d[[env]] <- factor(d[[env]])
  d[[rep]] <- factor(d[[rep]])
  gens <- levels(d[[gen]]); envs <- levels(d[[env]])
  g <- length(gens); e <- length(envs); r <- nlevels(d[[rep]])

  # genotype x environment mean matrix
  M <- tapply(d[[trait]], list(d[[gen]], d[[env]]), mean, na.rm = TRUE)
  M <- matrix(as.numeric(M), nrow = g, dimnames = list(gens, envs))
  grand <- mean(M)

  # combined ANOVA: Env, Rep(Env), Gen, Gen:Env, Error
  f <- stats::as.formula(paste0("`", trait, "` ~ ", env, " + ",
        env, ":", rep, " + ", gen, " + ", gen, ":", env))
  av  <- .tidy_aov(stats::aov(f, data = d))
  Me  <- av$MS[nrow(av)]; edf <- av$Df[nrow(av)]

  # Eberhart-Russell parameters
  gmean <- rowMeans(M)
  Ij    <- colMeans(M) - grand              # environmental index (centred)
  sumI2 <- sum(Ij^2)
  bi    <- as.numeric(M %*% Ij) / sumI2
  names(bi) <- gens
  pred  <- outer(gmean, rep(1, e)) + outer(bi, Ij)
  dev2  <- rowSums((M - pred)^2)
  S2di  <- dev2 / (e - 2) - Me / r
  er <- data.frame(genotype = gens, mean = gmean, bi = bi, S2di = S2di,
                   row.names = NULL, stringsAsFactors = FALSE)
  er <- er[order(-er$mean), ]; rownames(er) <- NULL

  # AMMI decomposition of the interaction matrix
  Rmat <- M - outer(gmean, rep(1, e)) -
          outer(rep(1, g), colMeans(M)) + grand
  sv   <- svd(Rmat)
  m    <- length(sv$d)
  gsc  <- sv$u %*% diag(sqrt(sv$d), m, m)
  esc  <- sv$v %*% diag(sqrt(sv$d), m, m)
  colnames(gsc) <- colnames(esc) <- paste0("IPCA", seq_len(m))
  rownames(gsc) <- gens; rownames(esc) <- envs
  ipca_pct <- 100 * sv$d^2 / sum(sv$d^2)

  structure(list(
    anova = av, eberhart = er, M = M, grand = grand,
    gen_scores = gsc, env_scores = esc, ipca_pct = ipca_pct,
    gen_mean = gmean, env_mean = colMeans(M),
    trait = trait, gens = gens, envs = envs, Me = Me, edf = edf,
    r = r), class = "bk_stability")
}

#' @export
print.bk_stability <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | GxE stability  |  trait: ", x$trait, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  a <- x$anova
  a$SS <- round(a$SS, 2); a$MS <- round(a$MS, 2)
  a$F  <- ifelse(is.na(a$F), "", formatC(a$F, format = "f", digits = 2))
  a$p  <- ifelse(is.na(a$p), "", formatC(a$p, format = "f", digits = 4))
  print(a, row.names = FALSE)
  cat("\nEberhart-Russell stability (ideal: bi ~ 1, S2di ~ 0):\n")
  er <- x$eberhart
  er$mean <- round(er$mean, 2); er$bi <- round(er$bi, 3)
  er$S2di <- round(er$S2di, 3)
  print(er, row.names = FALSE)
  cat("\nAMMI interaction captured:  IPCA1 =", round(x$ipca_pct[1], 1),
      "%   IPCA2 =", round(x$ipca_pct[2], 1), "%\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_stability <- function(x, type = c("ammi2", "ammi1", "er"), ...) {
  type <- match.arg(type)
  has_repel <- requireNamespace("ggrepel", quietly = TRUE)

  if (type == "er") {
    df <- x$eberhart
    ggplot2::ggplot(df, ggplot2::aes(.data$mean, .data$bi)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                          colour = "#B5341C") +
      ggplot2::geom_vline(xintercept = mean(df$mean), linetype = "dotted",
                          colour = "#6E6E6E") +
      ggplot2::geom_point(ggplot2::aes(colour = .data$S2di), size = 3.6) +
      (if (has_repel)
        ggrepel::geom_text_repel(ggplot2::aes(label = .data$genotype),
                                 size = 2.8)
       else ggplot2::geom_text(ggplot2::aes(label = .data$genotype),
                               size = 2.8, vjust = -0.9)) +
      ggplot2::scale_colour_gradientn(colours = bk_palette("field"),
                                      name = "S2di") +
      ggplot2::labs(
        title = "Eberhart-Russell stability",
        subtitle = "high mean + bi near 1 + low S2di = widely adapted & stable",
        x = paste0("mean ", x$trait), y = "regression coefficient (bi)",
        caption = "BKBreed") +
      theme_bk()

  } else if (type == "ammi1") {
    gd <- data.frame(label = x$gens, main = x$gen_mean,
                     ipca1 = x$gen_scores[, 1], kind = "genotype")
    ed <- data.frame(label = x$envs, main = x$env_mean,
                     ipca1 = x$env_scores[, 1], kind = "environment")
    df <- rbind(gd, ed)
    ggplot2::ggplot(df, ggplot2::aes(.data$main, .data$ipca1,
                                     colour = .data$kind)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#B5341C",
                          linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(shape = .data$kind), size = 3) +
      (if (has_repel)
        ggrepel::geom_text_repel(ggplot2::aes(label = .data$label), size = 2.7)
       else ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.7,
                               vjust = -0.9)) +
      ggplot2::scale_colour_manual(
        values = c(genotype = "#1B7340", environment = "#B5341C"), name = NULL) +
      ggplot2::labs(
        title = "AMMI-1 biplot",
        subtitle = paste0("main effect vs IPCA1 (",
                          round(x$ipca_pct[1], 1), "% of interaction)"),
        x = paste0("mean ", x$trait), y = "IPCA1", caption = "BKBreed") +
      theme_bk()

  } else { # ammi2
    gd <- data.frame(label = x$gens, ipca1 = x$gen_scores[, 1],
                     ipca2 = x$gen_scores[, 2], kind = "genotype")
    ed <- data.frame(label = x$envs, ipca1 = x$env_scores[, 1],
                     ipca2 = x$env_scores[, 2], kind = "environment")
    ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, colour = "#E4E4DC") +
      ggplot2::geom_vline(xintercept = 0, colour = "#E4E4DC") +
      ggplot2::geom_segment(data = ed,
        ggplot2::aes(x = 0, y = 0, xend = .data$ipca1, yend = .data$ipca2),
        colour = "#B5341C", linewidth = 0.5,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.16, "cm"))) +
      ggplot2::geom_point(data = gd,
        ggplot2::aes(.data$ipca1, .data$ipca2), colour = "#1B7340",
        size = 3.2) +
      (if (has_repel) list(
        ggrepel::geom_text_repel(data = gd,
          ggplot2::aes(.data$ipca1, .data$ipca2, label = .data$label),
          colour = "#0B3D2E", size = 2.7),
        ggrepel::geom_text_repel(data = ed,
          ggplot2::aes(.data$ipca1, .data$ipca2, label = .data$label),
          colour = "#B5341C", size = 2.7))
       else list(
        ggplot2::geom_text(data = gd,
          ggplot2::aes(.data$ipca1, .data$ipca2, label = .data$label),
          colour = "#0B3D2E", size = 2.7, vjust = -0.9),
        ggplot2::geom_text(data = ed,
          ggplot2::aes(.data$ipca1, .data$ipca2, label = .data$label),
          colour = "#B5341C", size = 2.7, vjust = -0.9))) +
      ggplot2::labs(
        title = "AMMI-2 biplot (IPCA1 vs IPCA2)",
        subtitle = paste0("green = genotypes, red = environments   |   IPCA1 ",
                          round(x$ipca_pct[1], 1), "% , IPCA2 ",
                          round(x$ipca_pct[2], 1), "%"),
        x = "IPCA1", y = "IPCA2", caption = "BKBreed") +
      theme_bk()
  }
}
