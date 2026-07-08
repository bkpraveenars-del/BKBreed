\name{bk_diversity}
\alias{bk_diversity}
\title{Mahalanobis D-square genetic divergence}
\usage{
bk_diversity(data, traits, gen, rep,
             method = c("tocher", "hierarchical"), clusters = NULL)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{traits}{Character vector of trait columns.}
  \item{gen}{Character; genotype column.}
  \item{rep}{Character; replication column.}
  \item{method}{Clustering for the default figure: tocher (default) or hierarchical.}
  \item{clusters}{Number of hierarchical clusters; if NULL, matches the Tocher count.}
}
\value{
An object of class \code{bk_diversity} with the \code{D2} matrix, cluster
memberships and cluster means.
}
\description{
Computes Mahalanobis D-square distances among genotypes using the pooled error
covariance, then groups them by Tocher and Ward clustering.
}
\examples{
dv <- bk_diversity(bk_data("rbd"),
        c("grain_yield","plant_height","tillers","test_weight"),
        "genotype", "rep")
dv
bk_plot(dv)
}
