#' StatHic
#' @keywords internal
#' @noRd
StatHic <- ggplot2::ggproto(
  "StatHic",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  setup_params = function(data, params) params,
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
      .calculateHicCoordinates()

    n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    env$n_annotation <- 0
    env$n_track <- 0
    env$n_concatemer <- 0
    env$n_hic_under <- 0
    env$n_hic <- 1
    env$max_y <- max(dat$ymax, na.rm = TRUE)
    env$min_y <- min(dat$y, na.rm = TRUE)
    env$MIN_Y <- min(dat$y, na.rm = TRUE)
    env$max_x <- max(dat$xend, na.rm = TRUE)
    env$min_x <- min(dat$xmin, na.rm = TRUE)
    env$res <- data$end1[1] - data$start1[1] + 1
    env$maxs_x <- dat |>
      dplyr::group_by(seqnames2) |>
      dplyr::summarize(maxs_x = max(xend)) |>
      dplyr::pull(maxs_x) |>
      stats::setNames(unique(dat$seqnames2))
    env$n_sn <- n_sn
    env$grs_range <- data |>
      dplyr::group_by(seqnames1) |>
      dplyr::summarize(
        start = min(start1),
        end = max(end1)
      ) |>
      dplyr::reframe(range = paste0(seqnames1, ":", start, "-", end)) |>
      dplyr::pull(range) |>
      GenomicRanges::GRanges()
    env$has_chr <- any(stringr::str_detect(data$seqnames1, "^chr")) ||
      any(stringr::str_detect(data$seqnames2, "^chr"))

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        .calculateAddLengths()
      chroms_sub <- data |>
        .calculateSubtractLengths()
      env$chroms_add <- chroms_add
      env$chroms_sub <- chroms_sub
    } else {
      env$chroms_add <- NULL
      env$chroms_sub <- NULL
    }

    dat
  }
)

#' GeomHic
#' @keywords internal
#' @noRd
GeomHic <- ggplot2::ggproto(
  "GeomHic",
  ggplot2::Geom,
  required_aes = c(
    "x", "xmin", "xmax", "xend", "y", "ymin", "ymax", "yend", "fill", "type"
  ),
  default_aes = ggplot2::aes(fill = NA, colour = NA, linetype = "dashed"),
  draw_key = ggplot2::draw_key_polygon,
  draw_panel = function(
    data, panel_params, coord, rasterize = TRUE, dpi = 300, dev = "cairo",
    scale = 1, draw_boundary = TRUE, boundary_colour = "black",
    linetype = "dashed"
  ) {
    coords <- coord$transform(data, panel_params)
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_sn <- env$n_sn
    grob_boundary_left <- grob_boundary_right <- grid::nullGrob()
    if (
      draw_boundary &&
        (n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))
    ) {
      coords_line_top <- coords |>
        dplyr::filter(type == "boundary") |>
        dplyr::group_by(seqnames1, seqnames2) |>
        dplyr::slice_max(order_by = ymax, n = 1, with_ties = TRUE) |>
        dplyr::slice_max(order_by = xmax, n = 1, with_ties = FALSE) |>
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

    coords_hic <- coords |>
      dplyr::filter(type != "boundary")
    grob_hic <- grid::polygonGrob(
      x = c(coords_hic$x, coords_hic$xmin, coords_hic$xmax, coords_hic$xend),
      y = c(coords_hic$y, coords_hic$ymin, coords_hic$ymax, coords_hic$yend),
      id = rep(seq_len(nrow(coords_hic)), 4),
      default.units = "native",
      gp = grid::gpar(fill = coords_hic$fill, col = NA)
    )

    if (rasterize) {
      class(grob_hic) <- c("rasteriser", class(grob_hic))
      grob_hic$dpi <- dpi
      grob_hic$dev <- dev
      grob_hic$scale <- scale
    }

    grid::gList(grob_hic, grob_boundary_left, grob_boundary_right)
  }
)

#' geom_hic
#'
#' @description
#' Creates Hi-C contact heatmap layer with diamond-shaped tiles.
#'
#' @inheritParams ggplot2::geom_polygon
#' @param mapping Aesthetic mappings from [ggplot2::aes()].
#' @param rasterize Logical. Rasterize for performance (default: TRUE).
#' @param dpi Numeric. Rasterization resolution (default: 300).
#' @param dev Character. Graphics device (default: `"cairo"`).
#' @param scale Numeric. Rasterization scaling (default: 1).
#' @param draw_boundary Logical. Draw chromosome boundaries (default: TRUE).
#' @param boundary_colour Character. Boundary color (default: `"black"`).
#' @param linetype Boundary line type (default: `"dashed"`).
#' @param ... Additional parameters (unused).
#'
#' @details
#' Transforms interaction data into diamond tiles at 45°, creating the
#' characteristic triangular Hi-C visualization.
#'
#' **Required aesthetics:** seqnames1, start1, end1, seqnames2, start2, end2,
#' fill
#'
#' @return A ggplot2 layer.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("file.cool") |> import()
#' library(ggplot2)
#' ggplot() +
#'   geom_hic(
#'     data = scaleData(cc["chr4"], "balanced", log10),
#'     aes(seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'         seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score)
#'   ) + theme_hic()
#' }
#'
#' @seealso [gghic()], [theme_hic()], [scaleData()]
#' @export
#' @aliases geom_hic
geom_hic <- function(
  mapping = NULL, data = NULL, stat = StatHic, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, rasterize = TRUE,
  dpi = 300, dev = "cairo", scale = 1, draw_boundary = TRUE,
  boundary_colour = "black", linetype = "dashed", ...
) {
  ggplot2::layer(
    geom = GeomHic, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, rasterize = rasterize, dpi = dpi, dev = dev, scale = scale,
      draw_boundary = draw_boundary, boundary_colour = boundary_colour,
      linetype = linetype, ...
    )
  )
}
