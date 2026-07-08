\name{bk_lxt}
\alias{bk_lxt}
\title{Line x Tester combining-ability analysis}
\usage{
bk_lxt(data, trait, line, tester, rep, alpha = 0.05)
}
\arguments{
  \item{data}{A data frame of the crosses in long format.}
  \item{trait}{Character; response column.}
  \item{line}{Character; line (female) column.}
  \item{tester}{Character; tester (male) column.}
  \item{rep}{Character; replication column.}
  \item{alpha}{Significance level for testing GCA/SCA effects (default 0.05).}
}
\value{
An object of class \code{bk_lxt} with the combining-ability ANOVA, GCA effects
for lines and testers, SCA effects, proportional contributions and variance
components.
}
\description{
Kempthorne's Line x Tester analysis for \eqn{l \times t} crosses in a randomised
block design. Partitions crosses into lines, testers and line x tester;
estimates GCA and SCA effects; and derives additive and dominance variances
(assuming inbred parents, F = 1).
}
\seealso{\code{\link{bk_plot}}}
\examples{
res <- bk_lxt(bk_data("lxt"), trait = "grain_yield",
              line = "line", tester = "tester", rep = "rep")
res
bk_plot(res)
bk_plot(res, type = "sca")
}
