\name{bk_augmented}
\alias{bk_augmented}
\alias{bk_augmented}
\title{Augmented Alpha-Lattice / Augmented RCBD analysis}
\usage{
bk_augmented(data, trait, gen, block, rep = NULL, checks = NULL, alpha = 0.05)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{trait}{Character; response column.}
  \item{gen}{Character; genotype column (checks and test entries).}
  \item{block}{Character; (incomplete) block column.}
  \item{rep}{Character or NULL; replication column for the alpha-lattice layout.}
  \item{checks}{Character vector of check names; if NULL, genotypes occurring more than once are treated as checks.}
  \item{alpha}{Significance level (default 0.05).}
}
\value{
An object of class \code{bk_augmented} with adjusted means, block effects, error MS and Federer critical differences.
}
\description{
Analysis of augmented designs. Block effects are estimated by intra-block least squares from replicated checks; every genotype mean is adjusted and classical Federer standard errors are reported.
}
\examples{
d <- bk_data("augmented")
res <- bk_augmented(d, "grain_yield", "genotype", "block", "rep",
                    checks = c("CHK-1","CHK-2","CHK-3","CHK-4"))
res
bk_plot(res)
}
