calculate_hic_coordinates <- function(data) {
  n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))

  if (n_sn == 1 || (n_sn == 2 && all(data$seqnames1 != data$seqnames2))) {
    dat <- data
  } else {
    dat <- data |>
      adjust_coordinates(
        list(
          c(start1 = "start1", end1 = "end1"),
          c(start2 = "start2", end2 = "end2")
        )
      )
  }

  dat <- dat |>
    dplyr::mutate(
      x = (end1 + start2) / 2,
      xmin = (start1 + start2) / 2,
      xmax = (start1 + end2) / 2,
      xend = (end1 + end2) / 2,
      y = (x - end1),
      ymin = (xmin - start1),
      ymax = (xmax - start1),
      yend = (xend - end1)
    )

  dat
}

StatHic <- ggplot2::ggproto(
  "StatHic",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(data, scales) {
    # ======================================================================== #
    #   ^          /\ (xmax, ymax)                                             #
    #   |         /  \                                                         #
    #   |        /    \                                                        #
    #   |       /      \                                                       #
    #   |      /        \                                                      #
    #   |     /          \                                                     #
    #   |    /            \                                                    #
    #   |   /              \                                                   #
    #   |   \ (xmin, ymin) / (xend, yend)                                      #
    #   |    \            /                                                    #
    #   |     \          /                                                     #
    #   |      \        /                                                      #
    #   |       \      /                                                       #
    #   |        \    /                                                        #
    #   |         \  /                                                         #
    #   |          \/ (x, y)                                                   #
    # --+--------------------------------------------------------------------> #
    #   |                                                                      #
    # ======================================================================== #
    dat <- data |>
      calculate_hic_coordinates()

    env <- get(".env", envir = asNamespace(name_pkg))
    env$n_annotation <- 0
    env$n_track <- 0
    env$n_hic <- 1
    env$max_y <- max(dat$ymax, na.rm = TRUE)
    env$min_y <- min(dat$y, na.rm = TRUE)
    env$MIN_Y <- min(dat$y, na.rm = TRUE)
    env$max_x <- max(dat$xend, na.rm = TRUE)
    env$min_x <- min(dat$xmin, na.rm = TRUE)
    env$res <- data$end1[1] - data$start1[1] + 1
    env$n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
    env$maxs_x <- dat |>
      dplyr::group_by(seqnames1) |>
      dplyr::summarize(maxs_x = max(xend)) |>
      dplyr::pull(maxs_x) |>
      stats::setNames(unique(dat$seqnames1))

    dat
  }
)

GeomHic <- ggplot2::ggproto(
  "GeomHic",
  ggplot2::Geom,
  required_aes = c(
    "x", "xmin", "xmax", "xend", "y", "ymin", "ymax", "yend", "fill"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params, "draw_boundary", "boundary_colour", "linetype"
  ),
  draw_key = ggplot2::draw_key_polygon,
  draw_panel = function(
    data, panel_params, coord,
    draw_boundary, boundary_colour, linetype
  ) {
    coords <- coord$transform(data, panel_params)
    env <- get(".env", envir = asNamespace(name_pkg))
    n_sn <- env$n_sn
    grob_boundary_left <- grob_boundary_right <- grid::nullGrob()
    if (
      draw_boundary &&
        (n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))
    ) {
      coords_line_top <- coords |>
        dplyr::group_by(seqnames1, seqnames2) |>
        dplyr::slice_max(order_by = ymax, n = 1) |>
        dplyr::ungroup() |>
        dplyr::filter(seqnames1 == seqnames2) |>
        dplyr::select(xmax, ymax)
      coords_line_br <- coords |>
        dplyr::group_by(seqnames1, seqnames2) |>
        dplyr::slice_max(order_by = xend, n = 1) |>
        dplyr::ungroup() |>
        dplyr::filter(seqnames1 == seqnames2) |>
        dplyr::select(xend, yend)
      coords_line_bl <- coords |>
        dplyr::group_by(seqnames1, seqnames2) |>
        dplyr::slice_min(order_by = xmin, n = 1) |>
        dplyr::ungroup() |>
        dplyr::filter(seqnames1 == seqnames2) |>
        dplyr::select(xmin, ymin)
      coords_line_left <- dplyr::bind_cols(coords_line_top, coords_line_br) |>
        dplyr::filter(xend < max(xmax))
      grob_boundary_left <- grid::polylineGrob(
        x = c(coords_line_left$xmax, coords_line_left$xend),
        y = c(coords_line_left$ymax, coords_line_left$yend),
        id = rep(seq_len(nrow(coords_line_left)), 2),
        gp = grid::gpar(col = boundary_colour, lty = linetype),
        default.units = "native"
      )
      coords_line_right <- dplyr::bind_cols(coords_line_top, coords_line_bl) |>
        dplyr::filter(xmin > min(xmin))
      grob_boundary_right <- grid::polylineGrob(
        x = c(coords_line_right$xmax, coords_line_right$xmin),
        y = c(coords_line_right$ymax, coords_line_right$ymin),
        id = rep(seq_len(nrow(coords_line_right)), 2),
        gp = grid::gpar(col = boundary_colour, lty = linetype),
        default.units = "native"
      )
    }
    grob_hic <- grid::polygonGrob(
      x = c(coords$x, coords$xmin, coords$xmax, coords$xend),
      y = c(coords$y, coords$ymin, coords$ymax, coords$yend),
      id = rep(seq_len(nrow(coords)), 4),
      default.units = "native",
      gp = grid::gpar(fill = coords$fill, col = NA)
    )
    grid::gList(grob_hic, grob_boundary_left, grob_boundary_right)
  }
)

#' geom_hic
#'
#' @description A ggplot2 geom for Hi-C data.
#' @inheritParams ggplot2::geom_polygon
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param draw_boundary Whether to draw the boundary line or not when plotting
#'   multiple chromosomes. Default is `TRUE`.
#' @param boundary_colour The color of the boundary line. Default is `"black"`.
#' @param linetype  The line type of the boundary line. Default is `"dashed"`.
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqnames1
#' * start1
#' * end1
#' * seqnames2
#' * start2
#' * end2
#' * fill
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' library(gghic)
#' library(ggplot2)
#' library(dplyr)
#' library(HiCExperiment)
#' library(InteractionSet)
#' library(scales)
#' library(glue)
#' library(rappdirs)
#'
#' download_example_files()
#' dir_cache_gghic <- user_cache_dir(appname = "gghic")
#'
#' hic <- glue("{dir_cache_gghic}/chr4_11-5kb.cool") |>
#'   CoolFile() |>
#'   import(cf)
#'
#' gis <- interactions(hic)
#' gis$score <- log10(gis$balanced)
#' x <- as_tibble(gis)
#' scores <- x$score[pairdist(gis) != 0 & !is.na(pairdist(gis) != 0)]
#' scores <- scores[!is.na(scores) & !is.infinite(scores)]
#' x$score <- oob_squish(x$score, c(min(scores), max(scores)))
#'
#' x |>
#'   filter(
#'     seqnames1 == "chr11", seqnames2 == "chr11",
#'     center1 > 67000000, center1 < 67100000,
#'     center2 > 67000000, center2 < 67100000
#'   ) |>
#'   ggplot(
#'     aes(
#'       seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'       seqnames2 = seqnames2, start2 = start2, end2 = end2,
#'       fill = score
#'     )
#'   ) +
#'   geom_hic() +
#'   theme_hic()
#' }
#' @export geom_hic
#' @aliases geom_hic
geom_hic <- function(
  mapping = NULL, data = NULL, stat = StatHic, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  draw_boundary = TRUE, boundary_colour = "black", linetype = "dashed"
) {
  ggplot2::layer(
    geom = GeomHic, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      draw_boundary = draw_boundary, boundary_colour = boundary_colour,
      linetype = linetype
    )
  )
}
