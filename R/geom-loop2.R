#' StatLoop2
#' @keywords internal
#' @noRd
StatLoop2 <- ggplot2::ggproto(
  "StatLoop2",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(data, scales) {
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic > 0) {
      n_sn <- env$n_sn
      chroms_add <- env$chroms_add
      chroms_sub <- env$chroms_sub
      MIN_Y <- env$MIN_Y
    } else {
      stop("geom_loop2() requires a HiC plot to be drawn first.")
    }

    if (!is.null(chroms_add) && !is.null(chroms_sub)) {
      data <- data |>
        dplyr::rename(seqname = seqnames1) |>
        .adjustCoordinates2(
          chroms_add, chroms_sub, c(start1 = "start1", end1 = "end1")
        ) |>
        dplyr::rename(seqname1 = seqname, seqname = seqnames2) |>
        .adjustCoordinates2(
          chroms_add, chroms_sub, c(start2 = "start2", end2 = "end2")
        ) |>
        dplyr::rename(seqname2 = seqname)
    }

    dat <- data |>
      dplyr::mutate(
        x = (start1 + end2) / 2,
        y = (x - start1),
        xmin = start1,
        xmax = end2,
        ymin = MIN_Y
      )

    dat
  }
)

#' GeomLoop2
#' @keywords internal
#' @noRd
GeomLoop2 <- ggplot2::ggproto(
  "GeomLoop2",
  ggplot2::Geom,
  required_aes = c("x", "y"),
  # draw_key = ggplot2::draw_key_point,
  default_aes = ggplot2::aes(
    colour = "black", shape = 21, fill = NA, size = 1, fontsize = 3,
    alpha = 1, stroke = 0.5
  ),
  draw_panel = function(
    data, panel_params, coord, style = "circle", n_arc_points = 50
  ) {
    coords <- coord$transform(data, panel_params)

    if (style == "arc") {
      grobs <- lapply(seq_len(nrow(coords)), function(i) {
        xmin <- coords$xmin[i]
        xmax <- coords$xmax[i]
        width <- xmax - xmin
        height <- coords$y[i] - coords$ymin[i]
        t <- seq(0, pi, length.out = n_arc_points)
        x_points <- xmin + width / 2 + (width / 2) * cos(t)
        y_points <- coords$ymin[i] + height * sin(t)

        grid::polylineGrob(
          x = x_points, y = y_points,
          default.units = "native",
          gp = grid::gpar(
            col = coords$colour[i],
            lwd = coords$stroke[i],
            alpha = coords$alpha[i]
          )
        )
      })
      do.call(grid::gList, grobs)
    } else {
      if (is.null(data$size)) {
        size <- grid::unit(1 / 80, "native")
      } else if (is.numeric(data$size)) {
        size <- grid::unit(data$size, "mm")
      } else if (grid::is.unit(data$size)) {
        size <- data$size
      }

      grid::pointsGrob(
        coords$x, coords$y,
        pch = coords$shape, size = size, default.units = "native",
        gp = grid::gpar(
          col = coords$colour, fill = coords$fill, lwd = coords$stroke
        )
      )
    }
  }
)

#' Visualize chromatin loops (a second version)
#'
#' @description
#' Another version of [geom_loop()] that requires direct data frame input
#' with loop coordinates instead of file paths. Supports both circle and arc
#' visualization styles.
#'
#' @inheritParams ggplot2::geom_point
#' @inheritParams geom_hic
#' @param style Character. Visualization style:
#'   * `"circle"`: Display loops as circular points (default)
#'   * `"arc"`: Display loops as curved arcs connecting anchors
#' @param n_arc_points Integer. Number of points for drawing smooth arcs when
#'   `style = "arc"`. Higher values = smoother curves (default: 50).
#' @param ... Additional parameters (unused).
#'
#' @details
#' ## Required aesthetics
#' * `seqnames1`, `start1`, `end1`: First loop anchor
#' * `seqnames2`, `start2`, `end2`: Second loop anchor
#'
#' ## Difference from geom_loop()
#' This version requires a Hi-C plot to be drawn first with `geom_hic()` and
#' takes data frames directly rather than file paths. Use [geom_loop()] for
#' more convenience with `gghic()` workflow.
#'
#' @return A ggplot2 layer.
#'
#' @seealso [geom_loop()], [gghic()]
#' @export
#' @examples
#' \dontrun{
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' # Prepare loop data frame
#' loop_file <- "path/to/loops.bedpe"
#' loops <- read.table(
#'   loop_file,
#'   col.names = c(
#'     "chr1", "start1", "end1", "chr2", "start2", "end2"
#'   )
#' )
#'
#' # Draw loops as circles (default)
#' gghic(cc) +
#'   geom_loop2(
#'     data = loops, aes(
#'       seqnames1 = chr1, start1 = start1, end1 = end1,
#'       seqnames2 = chr2, start2 = start2, end2 = end2
#'     )
#'   )
#'
#' # Draw loops as arcs
#' gghic(cc) +
#'   geom_loop2(
#'     data = loops, aes(
#'       seqnames1 = chr1, start1 = start1, end1 = end1,
#'       seqnames2 = chr2, start2 = start2, end2 = end2
#'     ), style = "arc"
#'   )
#' }
#' @export
#' @aliases geom_loop2
geom_loop2 <- function(
  mapping = NULL, data = NULL, stat = StatLoop2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
  style = "circle", n_arc_points = 50, ...
) {
  key_glyph <- if (style == "arc") {
    ggplot2::draw_key_path
  } else {
    ggplot2::draw_key_point
  }

  ggplot2::layer(
    geom = GeomLoop2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = TRUE,
    params = list(
      na.rm = na.rm, style = style, n_arc_points = n_arc_points, ...
    ),
    key_glyph = key_glyph
  )
}
