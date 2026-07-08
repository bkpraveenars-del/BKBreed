\name{bk_stability}
\alias{bk_stability}
\alias{bk_stability}
\title{GxE stability analysis (Eberhart-Russell + AMMI)}
\usage{
bk_stability(data, trait, gen, env, rep)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{trait}{Character; response column.}
  \item{gen}{Character; genotype column.}
  \item{env}{Character; environment/location column.}
  \item{rep}{Character; replication column.}
}
\value{
An object of class \code{bk_stability} with combined ANOVA, Eberhart-Russell parameters and AMMI scores.
}
\description{
Multi-location trial analysis: combined ANOVA, Eberhart-Russell regression (bi, S2di) and an AMMI decomposition with IPCA scores for biplots.
}
\examples{
st <- bk_stability(bk_data("mlt"), "grain_yield",
                   "genotype", "environment", "rep")
st
bk_plot(st)
bk_plot(st, type = "er")
}
