\name{bk_palette}
\alias{bk_palette}
\alias{bk_palette}
\title{BKBreed colour palettes}
\usage{
bk_palette(name = c("field", "spectrum", "canopy", "sunrise", "earth"),
          n = NULL, reverse = FALSE)
}
\arguments{
  \item{name}{Palette name.}
  \item{n}{Number of colours to return (interpolated if needed).}
  \item{reverse}{Logical; reverse the palette order.}
}
\value{
A character vector of hex colours.
}
\description{
A curated set of high-contrast palettes used across all BKBreed figures.
}
\examples{
bk_palette("sunrise", 5)
bk_palette("spectrum", 11)
}
