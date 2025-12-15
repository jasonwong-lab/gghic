.colorHic <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c(
      "#FFFEF9", "#FCF9CE", "#FFF2A9", "#FDE188", "#FFCA67",
      "#FAAA4B", "#F78E40", "#F15C34", "#ED3024", "#D42027",
      "#B01F29", "#7A1128", "#1A0A10"
    ),
    na.value = "#FFFFFF"
  )
}

#' ggplot2 theme optimized for Hi-C contact maps
#'
#' @description
#' Creates a clean, publication-ready theme for Hi-C visualization with
#' optimized defaults: hidden y-axis, fixed coordinate ratio, and custom color
#' gradient. Provides fine control over axis formatting and limits.
#'
#' @param hide_y Logical. Hide y-axis elements (axis title, text, ticks)
#'   (default: TRUE). Y-axis typically uninformative for Hi-C plots.
#' @param coord_ratio Numeric. Aspect ratio for [ggplot2::coord_fixed()]. 1
#'   maintains equal x/y scaling (default: 1).
#' @param scale_fill_gradientn Logical. Apply custom Hi-C color gradient
#'   (white → yellow → orange → red → black) (default: TRUE). Set FALSE to use
#'   custom scales.
#' @param breaks Numeric vector or [ggplot2::waiver()]. X-axis break positions
#'   for [ggplot2::scale_x_continuous()]. Default uses automatic breaks.
#' @param labels Function or character vector. X-axis label formatter. Default
#'   displays genomic coordinates in megabases (e.g., "10M" for 10,000,000 bp).
#' @param xmin Numeric or NULL. Minimum x-axis limit. If NULL, uses data range
#'   (default: NULL).
#' @param xmax Numeric or NULL. Maximum x-axis limit. If NULL, uses data range
#'   (default: NULL).
#' @param expand_x Numeric vector (length 2). X-axis expansion in plot units:
#'   `c(left_expansion, right_expansion)` (default: `c(0, 0)` for no
#'   expansion).
#'
#' @details
#' ## Theme components
#' * Based on `theme_bw()` with minimal grid lines
#' * Hidden y-axis (configurable)
#' * Fixed 1:1 coordinate ratio for accurate distance representation
#' * Custom color gradient optimized for Hi-C data
#' * Genomic coordinate formatting (Mb units)
#'
#' ## Color gradient
#' The default gradient transitions through 13 colors from white (no contacts)
#' to black (maximum contacts), providing excellent dynamic range for typical
#' Hi-C data.
#'
#' ## Customization
#' Combine with other ggplot2 theme elements for additional customization:
#' ```r
#' gghic(cc) + theme_hic() + theme(text = element_text(size = 14))
#' ```
#'
#' @return List of ggplot2 theme elements and scales that can be added to a
#'   ggplot object.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()
#'
#' # Basic usage with default theme
#' gghic(cc) + theme_hic()
#'
#' # Custom color scale
#' gghic(cc) +
#'   theme_hic(scale_fill_gradientn = FALSE) +
#'   scale_fill_viridis_c()
#'
#' # Show y-axis
#' gghic(cc) + theme_hic(hide_y = FALSE)
#'
#' # Custom axis limits and breaks
#' gghic(cc) +
#'   theme_hic(
#'     xmin = 1e6,
#'     xmax = 5e6,
#'     breaks = seq(1e6, 5e6, by = 1e6)
#'   )
#'
#' # Custom axis labels (show kb instead of Mb)
#' gghic(cc) +
#'   theme_hic(labels = scales::unit_format(unit = "kb", scale = 1e-3))
#' }
#'
#' @seealso [gghic()], [geom_hic()], [ggplot2::theme()]
#' @export
#' @aliases theme_hic
theme_hic <- function(
  hide_y = TRUE, coord_ratio = 1, scale_fill_gradientn = TRUE,
  breaks = ggplot2::waiver(),
  labels = scales::unit_format(unit = "M", scale = 1e-6), xmin = NULL,
  xmax = NULL, expand_x = c(0, 0)
) {
  t <- ggplot2::theme_bw() %+replace%
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )

  if (hide_y) {
    t <- t %+replace%
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_line(color = "black")
      )
  }

  if (scale_fill_gradientn) {
    tt <- list(t, .colorHic())
  }

  tt <- c(
    tt,
    ggplot2::scale_x_continuous(
      expand = c(0, 0), breaks = breaks, labels = labels
    ),
    ggplot2::coord_fixed(ratio = coord_ratio)
  )

  if (!is.null(xmin) && !is.null(xmax)) {
    tt <- c(
      tt,
      ggplot2::expand_limits(x = c(xmin - expand_x[1], xmax + expand_x[2]))
    )
  }

  tt
}
