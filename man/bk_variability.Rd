\name{bk_variability}
\alias{bk_variability}
\alias{bk_variability}
\title{Genetic variability parameters}
\usage{
bk_variability(data, traits, gen, rep, k = 2.063)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{traits}{Character vector of trait columns.}
  \item{gen}{Character; genotype column.}
  \item{rep}{Character; replication column.}
  \item{k}{Selection differential (default 2.063 for 5 percent intensity).}
}
\value{
An object of class \code{bk_variability} (a data frame with metadata).
}
\description{
Estimates GCV, PCV, ECV, broad-sense heritability, genetic advance and genetic advance as percent of mean for one or more traits.
}
\examples{
v <- bk_variability(bk_data("rbd"),
       c("grain_yield","plant_height","tillers","test_weight"),
       "genotype", "rep")
v
bk_plot(v)
}
