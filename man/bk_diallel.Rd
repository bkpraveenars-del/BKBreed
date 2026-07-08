\name{bk_diallel}
\alias{bk_diallel}
\title{Griffing diallel combining-ability analysis (Method 2, Model I)}
\usage{
bk_diallel(data, trait, parent1, parent2, rep, alpha = 0.05)
}
\arguments{
  \item{data}{A data frame in long format with one row per plot.}
  \item{trait}{Character; response column.}
  \item{parent1, parent2}{Character; the two parent columns of each entry.
    Diagonal entries (parent1 == parent2) are the parents/selfs.}
  \item{rep}{Character; replication column.}
  \item{alpha}{Significance level for testing effects (default 0.05).}
}
\value{
An object of class \code{bk_diallel} with the combining-ability ANOVA, parental
GCA effects, the SCA matrix, exact standard errors, Baker's predictability ratio
and variance components.
}
\description{
Griffing's (1956) diallel analysis for a half diallel of parents plus one set of
F1 crosses without reciprocals (Method 2), with genotypes as fixed effects
(Model I). Standard errors are computed exactly by propagating the per-entry
error variance through the linear GCA/SCA estimators.
}
\seealso{\code{\link{bk_lxt}}, \code{\link{bk_plot}}}
\examples{
res <- bk_diallel(bk_data("diallel"), trait = "grain_yield",
                  parent1 = "parent1", parent2 = "parent2", rep = "rep")
res
bk_plot(res)
bk_plot(res, type = "sca")
}
