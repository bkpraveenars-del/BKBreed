#' Mahalanobis D-square genetic-divergence analysis
#'
#' Computes Mahalanobis \eqn{D^2} distances among genotypes from a replicated
#' multi-trait trial using the pooled error variance-covariance matrix, then
#' groups genotypes by both Tocher's method and hierarchical (Ward) clustering.
#'
#' @param data A data frame in long format.
#' @param traits Character vector of trait columns.
#' @param gen Character; genotype column.
#' @param rep Character; replication column.
#' @param method Clustering used for the default figure and summary:
#'   \code{"tocher"} (default) or \code{"hierarchical"}.
#' @param clusters Number of hierarchical clusters. If \code{NULL}, matches the
#'   number of clusters found by Tocher's method.
#' @return An object of class \code{bk_diversity}.
#' @examples
#' dv <- bk_diversity(bk_data("rbd"),
#'         traits = c("grain_yield","plant_height","tillers",
#'                    "panicle_len","test_weight"),
#'         gen = "genotype", rep = "rep")
#' dv
#' bk_plot(dv)
#' @export
bk_diversity <- function(data, traits, gen, rep,
                         method = c("tocher", "hierarchical"),
                         clusters = NULL) {
  method <- match.arg(method)
  stopifnot(length(traits) >= 2, all(c(traits, gen, rep) %in% names(data)))
  d <- data[stats::complete.cases(data[c(traits, gen, rep)]), ]
  d[[gen]] <- factor(d[[gen]]); d[[rep]] <- factor(d[[rep]])
  gens <- levels(d[[gen]]); g <- length(gens); p <- length(traits)

  # genotype-mean matrix
  G <- sapply(traits, function(tr) tapply(d[[tr]], d[[gen]], mean, na.rm = TRUE))
  G <- matrix(as.numeric(G), nrow = g, dimnames = list(gens, traits))

  # pooled error variance-covariance from residuals of y ~ rep + gen
  res <- sapply(traits, function(tr) {
    f <- stats::as.formula(paste0("`", tr, "` ~ ", rep, " + ", gen))
    stats::residuals(stats::aov(f, data = d))
  })
  edf <- nrow(d) - g - nlevels(d[[rep]]) + 1
  Se  <- crossprod(res) / edf
  Sinv <- .safe_inv(Se)

  # D^2 matrix
  D2 <- matrix(0, g, g, dimnames = list(gens, gens))
  for (i in 1:(g - 1)) for (j in (i + 1):g) {
    dv <- G[i, ] - G[j, ]
    v  <- as.numeric(t(dv) %*% Sinv %*% dv)
    D2[i, j] <- D2[j, i] <- v
  }

  toch <- .tocher(D2)
  k    <- clusters %||% length(unique(toch))
  hc   <- stats::hclust(stats::as.dist(sqrt(D2)), method = "ward.D2")
  hcl  <- stats::cutree(hc, k = k)

  memb <- if (method == "tocher") toch else hcl
  # 2-D ordination of the D matrix for plotting
  ord <- tryCatch(stats::cmdscale(stats::as.dist(sqrt(D2)), k = 2),
                  error = function(e) matrix(0, g, 2))
  colnames(ord) <- c("Dim1", "Dim2")

  clustab <- data.frame(genotype = gens,
                        tocher = toch, hierarchical = hcl,
                        Dim1 = ord[, 1], Dim2 = ord[, 2],
                        row.names = NULL, stringsAsFactors = FALSE)

  # cluster mean vectors (using selected membership)
  cm <- t(sapply(sort(unique(memb)), function(cl)
    colMeans(G[memb == cl, , drop = FALSE])))
  rownames(cm) <- paste0("C", sort(unique(memb)))

  structure(list(
    D2 = D2, clusters = clustab, membership = memb, method = method,
    n_cluster = length(unique(memb)), cluster_means = cm,
    hc = hc, traits = traits, gens = gens,
    avg_D2 = mean(D2[upper.tri(D2)])), class = "bk_diversity")
}

# safe symmetric inverse (ridge fallback if near-singular)
.safe_inv <- function(S) {
  out <- try(solve(S), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  solve(S + diag(1e-6 * mean(diag(S)), nrow(S)))
}

# Tocher's method of cluster formation from a D^2 matrix
.tocher <- function(D2) {
  g <- nrow(D2); nm <- rownames(D2)
  theta <- max(vapply(seq_len(g), function(i) min(D2[i, -i]), numeric(1)))
  remaining <- seq_len(g)
  memb <- integer(g); cl <- 0L
  avg_intra <- function(S) if (length(S) < 2) 0 else
    mean(D2[t(utils::combn(S, 2))])
  while (length(remaining) > 0) {
    cl <- cl + 1L
    if (length(remaining) == 1) { memb[remaining] <- cl; break }
    sub <- D2[remaining, remaining, drop = FALSE]
    diag(sub) <- Inf
    wp <- which(sub == min(sub), arr.ind = TRUE)[1, ]
    C  <- remaining[c(wp[1], wp[2])]
    repeat {
      pool <- setdiff(remaining, C)
      if (length(pool) == 0) break
      old  <- avg_intra(C)
      incs <- vapply(pool, function(k) avg_intra(c(C, k)) - old, numeric(1))
      best <- which.min(incs)
      if (incs[best] <= theta) C <- c(C, pool[best]) else break
    }
    memb[C] <- cl
    remaining <- setdiff(remaining, C)
  }
  names(memb) <- nm
  memb
}

#' @export
print.bk_diversity <- function(x, ...) {
  cat("\n", strrep("=", 62), "\n", sep = "")
  cat("  BKBreed | Mahalanobis D-square divergence analysis\n")
  cat(strrep("=", 62), "\n\n", sep = "")
  cat("Genotypes:", length(x$gens), "  Traits:", length(x$traits),
      "  Avg D2:", round(x$avg_D2, 2), "\n")
  cat("Clustering shown:", x$method, " -> ", x$n_cluster, "clusters\n\n")
  tb <- table(x$membership)
  for (cl in names(tb)) {
    mem <- x$gens[x$membership == as.integer(cl)]
    cat("Cluster ", cl, " (", tb[cl], "): ",
        paste(mem, collapse = ", "), "\n", sep = "")
  }
  cat("\nCluster mean vectors:\n")
  print(round(x$cluster_means, 2))
  cat("\n")
  invisible(x)
}

#' @export
bk_plot.bk_diversity <- function(x, ...) {
  df <- x$clusters
  df$cluster <- factor(x$membership)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$Dim1, .data$Dim2,
                                        colour = .data$cluster))
  if (all(table(df$cluster) >= 3))
    p <- p + ggplot2::stat_ellipse(ggplot2::aes(fill = .data$cluster),
                                   geom = "polygon", alpha = 0.12,
                                   colour = NA, type = "norm", level = 0.68)
  p <- p +
    ggplot2::geom_hline(yintercept = 0, colour = "#E4E4DC") +
    ggplot2::geom_vline(xintercept = 0, colour = "#E4E4DC") +
    ggplot2::geom_point(size = 3.4)
  lab <- ggplot2::geom_text(ggplot2::aes(label = .data$genotype),
                            size = 2.7, vjust = -0.9, show.legend = FALSE)
  if (requireNamespace("ggrepel", quietly = TRUE))
    lab <- ggrepel::geom_text_repel(ggplot2::aes(label = .data$genotype),
                                    size = 2.7, show.legend = FALSE,
                                    max.overlaps = Inf)
  p + lab +
    scale_colour_bk("sunrise", name = "cluster") +
    scale_fill_bk("sunrise") +
    ggplot2::labs(
      title = "Genetic divergence — principal coordinates of D-square",
      subtitle = paste0(x$method, " clustering, ", x$n_cluster,
                        " clusters   |   avg D2 = ", round(x$avg_D2, 1)),
      x = "PCo 1", y = "PCo 2", caption = "BKBreed") +
    theme_bk()
}
# end of diversity module
