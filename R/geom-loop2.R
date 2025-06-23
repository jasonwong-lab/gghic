StatLoop2 <- ggplot2::ggproto(
  "StatLoop2",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(ggplot2::Stat$extra_params),
  compute_panel = function(data, scales) {
    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic > 0) {
      n_sn <- env$n_sn
      chroms_add <- env$chroms_add
      chroms_sub <- env$chroms_sub
    } else {
      stop("geom_loop2() requires a HiC plot to be drawn first.")
    }

    if (!is.null(chroms_add) && !is.null(chroms_sub)) {
      data <- data |>
        dplyr::rename(seqname = seqnames1) |>
        adjust_coordinates2(
          chroms_add, chroms_sub, c(start1 = "start1", end1 = "end1")
        ) |>
        dplyr::rename(seqname1 = seqname, seqname = seqnames2) |>
        adjust_coordinates2(
          chroms_add, chroms_sub, c(start2 = "start2", end2 = "end2")
        ) |>
        dplyr::rename(seqname2 = seqname)
    }

    dat <- data |>
      dplyr::mutate(
        x = (start1 + end2) / 2,
        y = (x - start1)
      )

    dat
  }
)

GeomLoop2 <- ggproto(
  "GeomLoop2",
  ggplot2::Geom,
  required_aes = c("x", "y"),
  extra_params = c(ggplot2::Geom$extra_params),
  draw_key = ggplot2::draw_key_point,
  default_aes = ggplot2::aes(
    colour = "black", shape = 21, fill = NA, size = 1, fontsize = 3, alpha = 1,
    stroke = 0.5
  ),
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)

    if (is.null(data$size)) {
      size <- grid::unit(1 / 80, "native")
    } else if (is.numeric(data$size)) {
      size <- grid::unit(data$size, "mm")
    } else if (grid::is.unit(data$size)) {
      size <- data$size
    }

    grid::pointsGrob(
      coords$x, coords$y,
      pch = coords$shape, size = size,
      default.units = "native",
      gp = grid::gpar(
        col = coords$colour, fill = coords$fill, lwd = coords$stroke
      )
    )
  }
)

#' geom_loop2
#'
#' @description A second version of [geom_loop()].
#' @inheritParams ggplot2::geom_point
#' @inheritParams geom_hic
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqnames1
#' * start1
#' * end1
#' * seqnames2
#' * start2
#' * end2
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # example usage
#' }
#' @export geom_loop2
#' @aliases geom_loop2
geom_loop2 <- function(
  mapping = NULL, data = NULL, stat = StatLoop2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, ...
) {
  ggplot2::layer(
    geom = GeomLoop2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = TRUE, params = list(na.rm = na.rm, ...)
  )
}
