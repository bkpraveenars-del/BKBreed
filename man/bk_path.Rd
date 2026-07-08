\name{bk_path}
\alias{bk_path}
\alias{bk_path}
\title{Path-coefficient analysis}
\usage{
bk_path(data, traits, dependent, gen, rep, type = c("genotypic", "phenotypic"))
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{traits}{Character vector of causal (independent) traits.}
  \item{dependent}{Character; the target/effect trait.}
  \item{gen}{Character; genotype column.}
  \item{rep}{Character; replication column.}
  \item{type}{Correlation basis: genotypic (default) or phenotypic.}
}
\value{
An object of class \code{bk_path} with direct/indirect effects and the residual.
}
\description{
Partitions the correlation of each causal trait with a target trait into direct and indirect effects (Wright's path analysis).
}
\examples{
pa <- bk_path(bk_data("rbd"),
        c("plant_height","tillers","panicle_len","test_weight"),
        "grain_yield", "genotype", "rep")
pa
bk_plot(pa)
}
