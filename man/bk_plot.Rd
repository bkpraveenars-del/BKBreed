\name{bk_plot}
\alias{bk_plot}
\alias{bk_plot}
\title{Draw the signature figure for a BKBreed result}
\usage{
bk_plot(x, ...)
}
\arguments{
  \item{x}{A BKBreed result object.}
  \item{...}{Passed to the specific method (e.g. \code{type} for stability and correlation).}
}
\value{
A \code{ggplot2} object.
}
\description{
Generic entry point returning the publication-ready figure for whichever analysis produced \code{x}.
}
\examples{
bk_plot(bk_rbd(bk_data("rbd"), "grain_yield", "genotype", "rep"))
}
