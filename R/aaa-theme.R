#' BKBreed colour palettes
#'
#' A curated set of high-contrast, colour-blind-considerate palettes used across
#' all BKBreed figures. These give the package its distinctive look.
#'
#' @param name Palette name. One of \code{"field"} (default, warm agronomic),
#'   \code{"spectrum"} (diverging blue-white-red, for correlations),
#'   \code{"canopy"} (green sequential), \code{"sunrise"} (categorical vivid),
#'   \code{"earth"} (muted categorical).
#' @param n Number of colours to return. If \code{NULL} the full anchor palette
#'   is returned; otherwise colours are interpolated to length \code{n}.
#' @param reverse Logical; reverse the palette order.
#' @return A character vector of hex colours.
#' @examples
#' bk_palette("sunrise", 5)
#' bk_palette("spectrum", 11)
#' @export
bk_palette <- function(name = c("field", "spectrum", "canopy",
                                "sunrise", "earth"), n = NULL,
                       reverse = FALSE) {
  name <- match.arg(name)
  anchors <- list(
    field    = c("#0B3D2E", "#1B7340", "#5BB35A", "#C7E29A",
                 "#F2C14E", "#E8792B", "#B5341C"),
    spectrum = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
                 "#FDDBC7", "#EF8A62", "#B2182B"),
    canopy   = c("#EDF8E9", "#BAE4B3", "#74C476", "#31A354", "#006D2C"),
    sunrise  = c("#264653", "#2A9D8F", "#E9C46A", "#F4A261",
                 "#E76F51", "#9B5DE5", "#00BBF9"),
    earth    = c("#8C6A5D", "#C5A880", "#6B8E23", "#4E6E58",
                 "#A63A50", "#3D5A80", "#D9A21B")
  )
  pal <- anchors[[name]]
  if (reverse) pal <- rev(pal)
  if (is.null(n)) return(pal)
  if (n <= length(pal) && name %in% c("sunrise", "earth")) {
    return(pal[seq_len(n)])
  }
  grDevices::colorRampPalette(pal)(n)
}

#' A clean, high-contrast ggplot2 theme for BKBreed
#'
#' @param base_size Base font size in points.
#' @param base_family Font family.
#' @param grid Logical; draw light major grid lines.
#' @return A \code{ggplot2} theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_bk()
#' @export
theme_bk <- function(base_size = 12, base_family = "", grid = TRUE) {
  ink    <- "#1A1A1A"
  panel  <- "#FCFCF9"
  subtle <- "#6E6E6E"
  gl     <- if (grid) ggplot2::element_line(colour = "#E4E4DC", linewidth = 0.3)
            else ggplot2::element_blank()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background  = ggplot2::element_rect(fill = panel, colour = NA),
      panel.grid.major  = gl,
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(fill = NA, colour = "#D8D8CF",
                                                linewidth = 0.4),
      axis.title        = ggplot2::element_text(colour = ink, face = "bold"),
      axis.text         = ggplot2::element_text(colour = subtle),
      plot.title        = ggplot2::element_text(colour = ink, face = "bold",
                                                size = base_size * 1.35,
                                                margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(colour = subtle,
                                                size = base_size * 0.95,
                                                margin = ggplot2::margin(b = 8)),
      plot.caption      = ggplot2::element_text(colour = subtle,
                                                size = base_size * 0.75),
      legend.title      = ggplot2::element_text(colour = ink, face = "bold"),
      legend.text       = ggplot2::element_text(colour = subtle),
      legend.position   = "right",
      strip.background  = ggplot2::element_rect(fill = "#0B3D2E", colour = NA),
      strip.text        = ggplot2::element_text(colour = "white", face = "bold",
                                                margin = ggplot2::margin(3,3,3,3))
    )
}

#' Discrete BKBreed colour and fill scales
#'
#' @param name BKBreed palette name (see \code{\link{bk_palette}}).
#' @param reverse Logical; reverse the palette.
#' @param ... Passed to \code{ggplot2::discrete_scale}.
#' @return A ggplot2 scale.
#' @rdname bk_scales
#' @export
scale_colour_bk <- function(name = "sunrise", reverse = FALSE, ...) {
  pal <- function(n) bk_palette(name, n = n, reverse = reverse)
  ggplot2::discrete_scale("colour", palette = pal, ...)
}

#' @rdname bk_scales
#' @export
scale_color_bk <- scale_colour_bk

#' @rdname bk_scales
#' @export
scale_fill_bk <- function(name = "sunrise", reverse = FALSE, ...) {
  pal <- function(n) bk_palette(name, n = n, reverse = reverse)
  ggplot2::discrete_scale("fill", palette = pal, ...)
}
