#' Griffing diallel combining-ability analysis (Method 2, Model I)
#'
#' Griffing's (1956) diallel analysis for a half diallel consisting of the
#' parents (selfs) plus one set of F1 crosses without reciprocals
#' (Method 2), treating genotypes as fixed effects (Model I). Returns the
#' combining-ability ANOVA (GCA and SCA), general combining-ability effects for
#' the parents, the specific combining-ability matrix, standard errors, Baker's
#' predictability ratio, and Model-II variance components.
#'
#' Standard errors are computed exactly, by propagating the per-entry error
#' variance through the linear GCA/SCA estimators, so they do not rely on
#' tabulated constants.
#'
#' @param data A data frame in long format with one row per plot.
#' @param trait Character; response column.
#' @param parent1,parent2 Character; the two parent columns of each entry.
#'   Diagonal entries (parent1 == parent2) are the parents/selfs.
#' @param rep Character; replication column.
#' @param alpha Significance level for testing effects (default 0.05).
#' @return An object of class \code{bk_diallel}.
#' @examples
#' res <- bk_diallel(bk_data("diallel"), trait = "grain_yield",
#'                   parent1 = "parent1", parent2 = "parent2", rep = "rep")
#' res
#' bk_plot(res)              # GCA effects
#' bk_plot(res, type = "sca")  # SCA heatmap
#' @export
bk_diallel <- function(data, trait, parent1, parent2, rep, alpha = 0.05) {
  stopifnot(all(c(trait, parent1, parent2, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(trait, parent1, parent2, rep)]), ]
  d[[parent1]] <- as.character(d[[parent1]])
  d[[parent2]] <- as.character(d[[parent2]])
  parents <- sort(unique(c(d[[parent1]], d[[parent2]])))
  n <- length(parents)
  if (n < 3) stop("A diallel needs at least 3 parents.")

  # unordered entry key so reciprocals collapse
  key <- apply(d[c(parent1, parent2)], 1,
               function(z) paste(sort(z), collapse = " x "))
  d$.entry <- key
  d[[rep]] <- factor(d[[rep]])
  r <- nlevels(d[[rep]])

  # entry-mean symmetric matrix
  em <- tapply(d[[trait]], d$.entry, mean, na.rm = TRUE)
  idx <- setNames(seq_len(n), parents)
  M <- matrix(NA_real_, n, n, dimnames = list(parents, parents))
  for (i in seq_len(n)) for (j in i:n) {
    k1 <- paste(sort(c(parents[i], parents[j])), collapse = " x ")
    if (!is.na(idx[parents[i]]) && k1 %in% names(em)) {
      M[i, j] <- M[j, i] <- em[[k1]]
    }
  }
  if (any(is.na(diag(M))))
    stop("Method 2 requires every parent (self) on the diagonal; ",
         "some selfs are missing. Crosses-only Method 4 is not yet supported.")
  if (any(is.na(M)))
    stop("The diallel table is incomplete (some crosses are missing).")

  k <- n * (n + 1) / 2
  uniq <- do.call(rbind, lapply(seq_len(n), function(i)
    cbind(i, i:n)))                         # k rows of (i,j), i<=j
  pos <- function(i, j) { a <- min(i, j); b <- max(i, j)
    which(uniq[, 1] == a & uniq[, 2] == b) }
  m_uniq <- apply(uniq, 1, function(z) M[z[1], z[2]])
  mdot <- sum(m_uniq); mu <- mdot / k
  mi <- rowSums(M); mii <- diag(M)

  # GCA and SCA (Griffing Method 2, Model I)
  g <- (mi + mii - (2 / n) * mdot) / (n + 2)
  S <- matrix(0, n, n, dimnames = list(parents, parents))
  for (i in seq_len(n)) for (j in i:n) {
    S[i, j] <- S[j, i] <- M[i, j] -
      (mi[i] + mii[i] + mi[j] + mii[j]) / (n + 2) +
      2 / ((n + 1) * (n + 2)) * mdot
  }

  # ANOVA: error MS from RBD of entries
  fit <- stats::aov(stats::as.formula(paste0("`", trait, "` ~ ", rep,
                    " + .entry")), data = d)
  av  <- .tidy_aov(fit)
  Me  <- av$MS[nrow(av)]; dfe <- av$Df[nrow(av)]

  SS_gca <- r * (1 / (n + 2)) * (sum((mi + mii)^2) - (4 / n) * mdot^2)
  SS_sca <- r * (sum(m_uniq^2) - (1 / (n + 2)) * sum((mi + mii)^2) +
                 2 / ((n + 1) * (n + 2)) * mdot^2)
  df_gca <- n - 1; df_sca <- n * (n - 1) / 2
  Mg <- SS_gca / df_gca; Ms <- SS_sca / df_sca
  mkrow <- function(src, ss, df) {
    ms <- ss / df; F <- ms / Me
    data.frame(Source = src, Df = df, SS = ss, MS = ms, F = F,
               p = stats::pf(F, df, dfe, lower.tail = FALSE),
               Sig = .stars(stats::pf(F, df, dfe, lower.tail = FALSE)),
               stringsAsFactors = FALSE)
  }
  anova_tab <- rbind(
    mkrow("GCA", SS_gca, df_gca),
    mkrow("SCA", SS_sca, df_sca),
    data.frame(Source = "Error", Df = dfe, SS = Me * dfe, MS = Me,
               F = NA_real_, p = NA_real_, Sig = "", stringsAsFactors = FALSE))

  # exact SEs via linear-combination coefficients over the k unique entries
  s2e <- Me / r
  coef_g <- function(i) {
    c <- numeric(k)
    for (j in seq_len(n)) c[pos(i, j)] <- c[pos(i, j)] + 1 / (n + 2)
    c[pos(i, i)] <- c[pos(i, i)] + 1 / (n + 2)
    c <- c - 2 / (n * (n + 2))
    c
  }
  coef_s <- function(i, j) {
    c <- numeric(k)
    c[pos(i, j)] <- c[pos(i, j)] + 1
    for (m in seq_len(n)) c[pos(i, m)] <- c[pos(i, m)] - 1 / (n + 2)
    c[pos(i, i)] <- c[pos(i, i)] - 1 / (n + 2)
    for (m in seq_len(n)) c[pos(j, m)] <- c[pos(j, m)] - 1 / (n + 2)
    c[pos(j, j)] <- c[pos(j, j)] - 1 / (n + 2)
    c <- c + 2 / ((n + 1) * (n + 2))
    c
  }
  se_g <- sqrt(s2e * sum(coef_g(1)^2))
  se_s_off <- sqrt(s2e * sum(coef_s(1, 2)^2))
  se_s_diag <- sqrt(s2e * sum(coef_s(1, 1)^2))

  tcrit <- stats::qt(1 - alpha / 2, dfe)
  gca <- data.frame(parent = parents, gca = as.numeric(g), se = se_g,
    sig = .stars(2 * stats::pt(abs(as.numeric(g) / se_g), dfe,
                               lower.tail = FALSE)),
    stringsAsFactors = FALSE)
  gca <- gca[order(-gca$gca), ]; rownames(gca) <- NULL

  baker <- 2 * Mg / (2 * Mg + Ms)
  s2gca <- (Mg - Me) / (n + 2)     # Model-II interpretation
  s2sca <- Ms - Me

  structure(list(
    anova = anova_tab, gca = gca, sca = S, mu = mu,
    se = list(gca = se_g, sca_offdiag = se_s_off, sca_diag = se_s_diag,
              gca_diff = sqrt(s2e * sum((coef_g(1) - coef_g(2))^2))),
    baker = baker, var_comp = c(gca = s2gca, sca = s2sca),
    Me = Me, dfe = dfe, n = n, r = r, alpha = alpha,
    trait = trait, parents = parents), class = "bk_diallel")
}

#' @export
print.bk_diallel <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Griffing diallel (Method 2, Model I)  |  trait: ",
      x$trait, "\n", sep = "")
  cat(strrep("=", 62), "\n\n", sep = "")
  cat("Parents:", x$n, "  Replications:", x$r, "  Grand mean:",
      round(x$mu, 2), "\n\n")
  a <- x$anova
  a$SS <- round(a$SS, 2); a$MS <- round(a$MS, 2)
  a$F  <- ifelse(is.na(a$F), "", formatC(a$F, format = "f", digits = 2))
  a$p  <- ifelse(is.na(a$p), "", formatC(a$p, format = "f", digits = 4))
  print(a, row.names = FALSE)
  cat("\nGCA effects  (SE=", round(x$se$gca, 3),
      ", SE of difference=", round(x$se$gca_diff, 3), "):\n", sep = "")
  gg <- x$gca; gg$gca <- round(gg$gca, 3)
  print(gg[c("parent", "gca", "sig")], row.names = FALSE)
  cat("\nSCA effects (SE off-diagonal=", round(x$se$sca_offdiag, 3), "):\n",
      sep = "")
  print(round(x$sca, 3))
  cat("\nBaker's predictability ratio:", round(x$baker, 3),
      "\nVariance components (Model II):  sigma2 GCA =",
      round(x$var_comp["gca"], 3), "  sigma2 SCA =",
      round(x$var_comp["sca"], 3), "\n")
  cat("\nSignif: *** p<0.001  ** p<0.01  * p<0.05  . p<0.10  ns >=0.10\n\n")
  invisible(x)
}

#' @export
bk_plot.bk_diallel <- function(x, type = c("gca", "sca"), ...) {
  type <- match.arg(type)
  if (type == "gca") {
    df <- x$gca
    df$sign <- ifelse(df$gca >= 0, "positive", "negative")
    df$parent <- factor(df$parent, levels = rev(df$parent))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$gca, y = .data$parent,
                                     fill = .data$sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data$gca - .data$se, xmax = .data$gca + .data$se),
        width = 0.25, colour = "#333333", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = 0, colour = "#1A1A1A", linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = c(positive = "#1B7340", negative = "#B5341C"), name = NULL) +
      ggplot2::labs(
        title = paste0("Diallel GCA effects - ", x$trait),
        subtitle = paste0("positive = good general combiner (bars = +/- 1 SE)",
                          "   Baker's ratio = ", round(x$baker, 2)),
        x = "GCA effect", y = NULL, caption = "BKBreed") +
      theme_bk() +
      ggplot2::theme(legend.position = "none")
  } else {
    m <- x$sca
    long <- expand.grid(p1 = rownames(m), p2 = colnames(m),
                        stringsAsFactors = FALSE)
    long$value <- as.vector(m)
    long$p1 <- factor(long$p1, levels = rev(rownames(m)))
    long$p2 <- factor(long$p2, levels = colnames(m))
     up <- as.integer(long$p2) >= (nrow(m) - as.integer(long$p1) + 1)
    long$lab <- ifelse(up, formatC(long$value, format = "f", digits = 2), "")
    ggplot2::ggplot(long, ggplot2::aes(.data$p2, .data$p1, fill = .data$value)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = .data$lab), size = 2.9,
                         colour = "#1A1A1A") +
      ggplot2::scale_fill_gradientn(colours = bk_palette("spectrum"),
                                    name = "SCA") +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = paste0("Diallel SCA effects - ", x$trait),
        subtitle = "diagonal = parental SCA; off-diagonal = cross SCA",
        x = NULL, y = NULL, caption = "BKBreed") +
      theme_bk() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  }
}
