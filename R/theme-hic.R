colors_hic <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c(
      "#FFFEF9", "#FCF9CE", "#FFF2A9", "#FDE188", "#FFCA67", "#FAAA4B",
      "#F78E40", "#F15C34", "#ED3024", "#D42027", "#B01F29", "#7A1128",
      "#1A0A10"
    ),
    na.value = "#FFFFFF"
  )
}

#' theme_hic
#'
#' @description Generate a ggplot2 theme for Hi-C plots.
#' @param hide_y A logical value indicating whether to hide the y-axis.
#'   Default is `TRUE`.
#' @param coord_ratio The ratio of the [ggplot2::coord_fixed]. Default is `1`.
#' @param breaks This argument is passed to [ggplot2::scale_x_continuous].
#' @param labels This argument is passed to [ggplot2::scale_x_continuous].
#' @param xmin The minimum x-axis value. Default is `NULL`.
#' @param xmax The maximum x-axis value. Default is `NULL`.
#' @param expand_x A numeric vector of length 2 to expand the x-axis.
#'   Default is `c(0, 0)`.
#' @details
#' If either `xmin` or `xmax` is `NULL`, the x-axis will not be expanded.
#' @return A ggplot2 theme.
#' @inherit geom_hic examples
#' @export theme_hic
theme_hic <- function(
  hide_y = TRUE, coord_ratio = 1,
  breaks = ggplot2::waiver(),
  labels = scales::unit_format(unit = "M", scale = 1e-6),
  xmin = NULL, xmax = NULL, expand_x = c(0, 0)
) {
  `%+replace%` <- ggplot2::`%+replace%`

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

  tt <- list(
    t,
    colors_hic(),
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
