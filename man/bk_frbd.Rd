\name{bk_frbd}
\alias{bk_frbd}
\alias{bk_frbd}
\title{Two-factor Factorial RBD analysis}
\usage{
bk_frbd(data, trait, factorA, factorB, rep, alpha = 0.05)
}
\arguments{
  \item{data}{A data frame in long format.}
  \item{trait}{Character; response column.}
  \item{factorA, factorB}{Character; the two treatment factor columns.}
  \item{rep}{Character; replication/block column.}
  \item{alpha}{Significance level (default 0.05).}
}
\value{
An object of class \code{bk_frbd}.
}
\description{
ANOVA for a two-factor factorial in a randomised block design, with main-effect and interaction means and factor-specific critical differences.
}
\examples{
res <- bk_frbd(bk_data("frbd"), "grain_yield", "nitrogen", "variety", "rep")
res
bk_plot(res)
}
