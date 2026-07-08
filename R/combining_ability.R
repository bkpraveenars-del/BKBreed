#' Line x Tester combining-ability analysis
#'
#' Kempthorne's Line x Tester analysis for a set of \eqn{l \times t} crosses
#' evaluated in a randomised block design. Returns the combining-ability ANOVA
#' (crosses partitioned into lines, testers and line x tester), general combining
#' ability (GCA) effects for lines and testers, specific combining ability (SCA)
#' effects for each cross, the proportional contribution of lines, testers and
#' their interaction, and the GCA/SCA variance components with derived additive
#' and dominance variances.
#'
#' Variance components assume homozygous (inbred) parents (F = 1), so
#' \eqn{\sigma^2_A = 2\sigma^2_{gca}} and \eqn{\sigma^2_D = \sigma^2_{sca}}.
#'
#' @param data A data frame of the crosses in long format.
#' @param trait Character; response column.
#' @param line Character; line (female) column.
#' @param tester Character; tester (male) column.
#' @param rep Character; replication column.
#' @param alpha Significance level for testing GCA/SCA effects (default 0.05).
#' @return An object of class \code{bk_lxt}.
#' @examples
#' res <- bk_lxt(bk_data("lxt"), trait = "grain_yield",
#'               line = "line", tester = "tester", rep = "rep")
#' res
#' bk_plot(res)              # GCA effects
#' bk_plot(res, type = "sca")  # SCA heatmap
#' @export
bk_lxt <- function(data, trait, line, tester, rep, alpha = 0.05) {
  stopifnot(all(c(trait, line, tester, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(trait, line, tester, rep)]), ]
  d[[line]]   <- factor(d[[line]])
  d[[tester]] <- factor(d[[tester]])
  d[[rep]]    <- factor(d[[rep]])
  l <- nlevels(d[[line]]); t <- nlevels(d[[tester]]); r <- nlevels(d[[rep]])
  grand <- mean(d[[trait]], na.rm = TRUE)

  fit <- stats::aov(stats::as.formula(paste0("`", trait, "` ~ ", rep,
           " + ", line, " * ", tester)), data = d)
  av  <- .tidy_aov(fit)
  pick <- function(nm) which(trimws(av$Source) == nm)
  i_rep <- pick(rep); i_l <- pick(line); i_t <- pick(tester)
  i_lt <- pick(paste0(line, ":", tester)); i_e <- nrow(av)

  Me  <- av$MS[i_e]; dfe <- av$Df[i_e]
  ss_line <- av$SS[i_l]; ss_test <- av$SS[i_t]; ss_lt <- av$SS[i_lt]
  ss_cross <- ss_line + ss_test + ss_lt; df_cross <- l * t - 1
  Ml <- av$MS[i_l]; Mt <- av$MS[i_t]; Mlt <- av$MS[i_lt]

  mkrow <- function(src, ss, df) {
    ms <- ss / df; F <- ms / Me
    data.frame(Source = src, Df = df, SS = ss, MS = ms, F = F,
               p = stats::pf(F, df, dfe, lower.tail = FALSE),
               Sig = .stars(stats::pf(F, df, dfe, lower.tail = FALSE)),
               stringsAsFactors = FALSE)
  }
  anova_tab <- rbind(
    mkrow("Replication", av$SS[i_rep], av$Df[i_rep]),
    mkrow("Crosses",     ss_cross, df_cross),
    mkrow("  Lines",     ss_line, l - 1),
    mkrow("  Testers",   ss_test, t - 1),
    mkrow("  Line x Tester", ss_lt, (l - 1) * (t - 1)),
    data.frame(Source = "Error", Df = dfe, SS = Me * dfe, MS = Me,
               F = NA_real_, p = NA_real_, Sig = "", stringsAsFactors = FALSE))

  # GCA / SCA from cross totals
  Yij <- tapply(d[[trait]], list(d[[line]], d[[tester]]), sum)
  Yij[is.na(Yij)] <- 0
  Yi <- rowSums(Yij); Yj <- colSums(Yij)
  gi <- Yi / (t * r) - grand
  gj <- Yj / (l * r) - grand
  Mij <- Yij / r
  sij <- sweep(sweep(Mij, 1, Yi / (t * r)), 2, Yj / (l * r)) + grand

  se_gl <- sqrt(Me / (r * t)); se_gt <- sqrt(Me / (r * l)); se_s <- sqrt(Me / r)
  tcrit <- stats::qt(1 - alpha / 2, dfe)
  gca_lines <- data.frame(parent = names(gi), type = "line",
    gca = as.numeric(gi), se = se_gl,
    sig = .stars(2 * stats::pt(abs(as.numeric(gi) / se_gl), dfe,
                               lower.tail = FALSE)),
    stringsAsFactors = FALSE)
  gca_testers <- data.frame(parent = names(gj), type = "tester",
    gca = as.numeric(gj), se = se_gt,
    sig = .stars(2 * stats::pt(abs(as.numeric(gj) / se_gt), dfe,
                               lower.tail = FALSE)),
    stringsAsFactors = FALSE)

  contrib <- c(lines = 100 * ss_line / ss_cross,
               testers = 100 * ss_test / ss_cross,
               line_tester = 100 * ss_lt / ss_cross)

  s2gca_l <- (Ml - Mlt) / (r * t)
  s2gca_t <- (Mt - Mlt) / (r * l)
  s2gca   <- ((Ml - Mlt) + (Mt - Mlt)) / (r * (l + t))
  s2sca   <- (Mlt - Me) / r
  s2A <- 2 * s2gca; s2D <- s2sca
  avg_dom <- if (s2A > 0) sqrt(2 * s2sca / s2A) else NA_real_
  predict <- if ((2 * s2gca + s2sca) != 0) (2 * s2gca) / (2 * s2gca + s2sca)
             else NA_real_

  structure(list(
    anova = anova_tab, gca_lines = gca_lines, gca_testers = gca_testers,
    sca = sij, contribution = contrib,
    se = list(gca_line = se_gl, gca_tester = se_gt, sca = se_s,
              gi_gj_line = sqrt(2 * Me / (r * t)),
              gi_gj_tester = sqrt(2 * Me / (r * l)),
              sij_skl = sqrt(2 * Me / r)),
    var_comp = c(gca_line = s2gca_l, gca_tester = s2gca_t, gca = s2gca,
                 sca = s2sca, sigma2A = s2A, sigma2D = s2D,
                 avg_dominance = avg_dom, predictability = predict),
    Me = Me, dfe = dfe, grand_mean = grand, alpha = alpha,
    l = l, t = t, r = r, trait = trait,
    line = line, tester = tester), class = "bk_lxt")
}

#' @export
print.bk_lxt <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Line x Tester combining ability  |  trait: ", x$trait,
      "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  a <- x$anova
  a$SS <- round(a$SS, 2); a$MS <- round(a$MS, 2)
  a$F  <- ifelse(is.na(a$F), "", formatC(a$F, format = "f", digits = 2))
  a$p  <- ifelse(is.na(a$p), "", formatC(a$p, format = "f", digits = 4))
  print(a, row.names = FALSE)
  cat("\nError MS:", round(x$Me, 3), " (df=", x$dfe, ")   Grand mean:",
      round(x$grand_mean, 2), "\n")

  cat("\nGCA effects - LINES  (SE=", round(x$se$gca_line, 3), "):\n", sep = "")
  gl <- x$gca_lines; gl$gca <- round(gl$gca, 3)
  print(gl[c("parent", "gca", "sig")], row.names = FALSE)
  cat("\nGCA effects - TESTERS  (SE=", round(x$se$gca_tester, 3), "):\n", sep = "")
  gt <- x$gca_testers; gt$gca <- round(gt$gca, 3)
  print(gt[c("parent", "gca", "sig")], row.names = FALSE)

  cat("\nSCA effects (line x tester;  SE=", round(x$se$sca, 3), "):\n", sep = "")
  print(round(x$sca, 3))

  cat("\nProportional contribution to crosses SS (%):\n")
  print(round(x$contribution, 2))

  cat("\nVariance components (inbred parents, F=1):\n")
  vc <- x$var_comp
  cat("  sigma2 GCA =", round(vc["gca"], 3),
      "   sigma2 SCA =", round(vc["sca"], 3), "\n")
  cat("  sigma2 A   =", round(vc["sigma2A"], 3),
      "   sigma2 D   =", round(vc["sigma2D"], 3), "\n")
  cat("  avg. degree of dominance =", round(vc["avg_dominance"], 3),
      "   predictability ratio =", round(vc["predictability"], 3), "\n")
  cat("\nSignif: *** p<0.001  ** p<0.01  * p<0.05  . p<0.10  ns >=0.10\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_lxt <- function(x, type = c("gca", "sca"), ...) {
  type <- match.arg(type)
  if (type == "gca") {
    df <- rbind(x$gca_lines, x$gca_testers)
    df$sign <- ifelse(df$gca >= 0, "positive", "negative")
    df$type <- factor(df$type, levels = c("line", "tester"),
                      labels = c("Lines", "Testers"))
    df <- df[order(df$type, df$gca), ]
    df$parent <- factor(df$parent, levels = df$parent)
    ggplot2::ggplot(df, ggplot2::aes(x = .data$gca, y = .data$parent,
                                     fill = .data$sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$gca - .data$se, xmax = .data$gca + .data$se),
        width = 0.25, colour = "#333333", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, colour = "#1A1A1A", linewidth = 0.4) +
      ggplot2::facet_grid(rows = ggplot2::vars(.data$type), scales = "free_y",
                          space = "free_y") +
      ggplot2::scale_fill_manual(
        values = c(positive = "#1B7340", negative = "#B5341C"), name = NULL) +
      ggplot2::labs(
        title = paste0("General combining ability - ", x$trait),
        subtitle = "positive GCA = good general combiner (bars = +/- 1 SE)",
        x = "GCA effect", y = NULL, caption = "BKBreed") +
      theme_bk()
  } else {
    m <- x$sca
    long <- expand.grid(line = rownames(m), tester = colnames(m),
                        stringsAsFactors = FALSE)
    long$value <- as.vector(m)
    long$line   <- factor(long$line, levels = rev(rownames(m)))
    long$tester <- factor(long$tester, levels = colnames(m))
    ggplot2::ggplot(long, ggplot2::aes(.data$tester, .data$line,
                                       fill = .data$value)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = formatC(.data$value,
                         format = "f", digits = 2)), size = 3,
                         colour = "#1A1A1A") +
      ggplot2::scale_fill_gradientn(colours = bk_palette("spectrum"),
                                    name = "SCA") +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = paste0("Specific combining ability - ", x$trait),
        subtitle = "red = superior specific combination for that cross",
        x = "tester", y = "line", caption = "BKBreed") +
      theme_bk() +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
}
