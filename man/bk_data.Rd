\name{bk_data}
\alias{bk_data}
\alias{bk_data}
\title{Load a bundled BKBreed example dataset}
\usage{
bk_data(name = c("rbd", "frbd", "augmented", "mlt"))
}
\arguments{
  \item{name}{One of rbd, frbd, augmented, mlt.}
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
