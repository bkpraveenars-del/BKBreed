\name{theme_bk}
\alias{theme_bk}
\alias{theme_bk}
\title{A clean, high-contrast ggplot2 theme}
\usage{
theme_bk(base_size = 12, base_family = "", grid = TRUE)
}
\arguments{
  \item{base_size}{Base font size in points.}
  \item{base_family}{Font family.}
  \item{grid}{Logical; draw light major grid lines.}
}
\value{
A \code{ggplot2} theme object.
}
\description{
The BKBreed figure theme.
}
\examples{
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_bk()
}
