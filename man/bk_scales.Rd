\name{scale_colour_bk}
\alias{scale_colour_bk}
\alias{scale_color_bk}
\alias{scale_fill_bk}
\title{Discrete BKBreed colour and fill scales}
\usage{
scale_colour_bk(pal_name = "sunrise", reverse = FALSE, ...)

scale_color_bk(pal_name = "sunrise", reverse = FALSE, ...)

scale_fill_bk(pal_name = "sunrise", reverse = FALSE, ...)
}
\arguments{
  \item{pal_name}{BKBreed palette name (see \code{\link{bk_palette}}).}
  \item{reverse}{Logical; reverse the palette.}
  \item{...}{Passed to \code{ggplot2::discrete_scale} (e.g. \code{name} for the
    legend title).}
}
\value{
A ggplot2 scale.
}
\description{
Discrete colour and fill scales built on the BKBreed palettes.
}
\examples{
library(ggplot2)
ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
  geom_point() + scale_colour_bk("sunrise")
}
