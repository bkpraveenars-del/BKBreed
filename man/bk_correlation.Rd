\name{bk_correlation}
\alias{bk_correlation}
\alias{bk_correlation}
\title{Genotypic and phenotypic correlation}
\usage{
bk_correlation(data, traits, gen, rep)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{traits}{Character vector of trait columns (>= 2).}
  \item{gen}{Character; genotype column.}
  \item{rep}{Character; replication column.}
}
\value{
An object of class \code{bk_correlation} with matrices \code{rg} and \code{rp}.
}
\description{
Genotypic and phenotypic correlation coefficients among traits, via analysis of covariance.
}
\examples{
cr <- bk_correlation(bk_data("rbd"),
        c("grain_yield","plant_height","tillers","test_weight"),
        "genotype", "rep")
cr
bk_plot(cr)
}
