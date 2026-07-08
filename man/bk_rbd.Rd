\name{bk_rbd}
\alias{bk_rbd}
\alias{bk_rbd}
\title{Randomised Block Design (RBD) analysis}
\usage{
bk_rbd(data, trait, gen, rep, alpha = 0.05)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{trait}{Character; response column.}
  \item{gen}{Character; genotype/treatment column.}
  \item{rep}{Character; replication/block column.}
  \item{alpha}{Significance level for the critical difference (default 0.05).}
}
\value{
An object of class \code{bk_rbd} with elements \code{anova}, \code{means}, \code{cv}, \code{sem}, \code{sed} and \code{cd}.
}
\description{
One-call ANOVA for a randomised complete block design with SE, critical difference, CV and compact-letter groupings.
}
\seealso{\code{\link{bk_plot}}, \code{\link{bk_frbd}}}
\examples{
res <- bk_rbd(bk_data("rbd"), "grain_yield", "genotype", "rep")
res
bk_plot(res)
}
