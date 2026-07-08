\name{bk_data}
\alias{bk_data}
\title{Load a bundled BKBreed example dataset}
\usage{
bk_data(name = c("rbd", "frbd", "augmented", "mlt", "lxt", "diallel"))
}
\arguments{
  \item{name}{One of \code{"rbd"}, \code{"frbd"}, \code{"augmented"},
    \code{"mlt"}, \code{"lxt"}, \code{"diallel"}.}
}
\value{
A data frame.
}
\description{
Convenience loader for the demonstration datasets shipped with the package.
}
\examples{
head(bk_data("rbd"))
}
