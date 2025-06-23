StatTrack2 <- ggplot2::ggproto(
  "StatTrack2",
  ggplot2::Stat,
  required_aes = c(
    "seqnames", "start", "end", "score", "name"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "width_ratio", "spacing_ratio", "data_range"
  ),
  dropped_aes = c("seqnames"),
  compute_panel = function(
    data, scales,
    width_ratio = 1 / 20, spacing_ratio = 0.5, data_range = c("auto", "maximum")
  ) {
    n_data <- length(unique(data$name))
    if (
      is.numeric(data_range) &&
        (length(data_range) != 1 && length(data_range) != n_data)
    ) {
      stop("data_range must be of length 1 or equal to the number of data.")
    }
    tracks <- split(data, data$name) |>
      purrr::map(
        function(x) {
          x |>
            dplyr::arrange(seqnames, start, end)
        }
      )

    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    chroms_add <- chroms_sub <- NULL
    if (env$n_hic == 1) {
      max_y <- env$max_y
      min_x <- env$min_x
      res <- env$res
      n_sn <- env$n_sn
      chroms_add <- env$chroms_add
      chroms_sub <- env$chroms_sub
    } else {
      max_y <- 10000
      res <- 0
      n_sn <- length(unique(c(data$seqnames)))

      if (n_sn > 1) {
        chroms_add <- data |>
          calculate_add_lengths_1d()
        chroms_sub <- data |>
          calculate_subtract_lengths_1d()
      }

      seqname_1st <- tracks[[1]]$seqnames[1]
      min_x <- min(data$start[data$seqnames == seqname_1st], na.rm = TRUE)

      if (n_sn > 1) {
        min_x <- min_x - chroms_sub$orignal_start[
          chroms_sub$seqname == seqname_1st
        ]
      }
    }
    if (n_sn > 1) {
      maxs_x <- data |>
        dplyr::group_by(seqnames) |>
        dplyr::summarize(maxs_x = max(end)) |>
        dplyr::pull(maxs_x) |>
        stats::setNames(unique(data$seqnames)) |>
        as.data.frame() |>
        setNames("x_max") |>
        tibble::rownames_to_column("seqname") |>
        adjust_coordinates2(chroms_add, chroms_sub, c(x_max = "x_max"))
      maxs_x <- maxs_x$x_max |>
        setNames(maxs_x$seqname)
    }
    min_y <- ifelse(
      n_annotation > n_track || n_concatemer > n_track, env$min_y, 0
    )

    .height <- max_y * width_ratio

    names_track <- unique(data$name)

    ys <- names_track |>
      tibble::as_tibble() |>
      dplyr::rename(name = value) |>
      dplyr::mutate(
        y = min_y -
          res / 2 -
          (.height * (dplyr::row_number())) -
          ((dplyr::row_number()) * (.height * spacing_ratio))
      )

    if (data_range[1] == "auto") {
      .max <- data |>
        dplyr::group_by(name) |>
        dplyr::summarize(.max = round(max(score, na.rm = TRUE), digits = 2))
    }
    if (data_range[1] == "maximum") {
      .max <- tibble::tibble(
        name = unique(data$name),
        .max = round(max(data$score, na.rm = TRUE), digits = 2)
      )
    }
    if (is.numeric(data_range)) {
      .max <- tibble::tibble(name = unique(data$name), .max = data_range)
    }

    dat_track <- data |>
      dplyr::rename(seqname = seqnames) |>
      dplyr::left_join(ys, by = "name") |>
      dplyr::left_join(.max, by = "name") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymin = y + ((.height / .max) * score),
        type = "track"
      ) |>
      dplyr::mutate(
        ymin = ifelse(ymin - y > .height, y + .height, ymin)
      )

    if (!is.null(chroms_add) && !is.null(chroms_sub)) {
      dat_track <- dat_track |>
        adjust_coordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))
    }

    dat_axis <- dat_track |>
      dplyr::distinct(name, .keep_all = TRUE) |>
      dplyr::mutate(
        ymin = y + .height,
        x = min_x,
        text_ymin = as.character(.max),
        type = "axis"
      )

    dat_text <- dat_axis |>
      dplyr::mutate(
        y = (y + ymin) / 2,
        type = "text"
      )

    dat <- dplyr::bind_rows(dat_track, dat_axis, dat_text)

    env$min_y <- min(dat$y)
    env$n_track <- env$n_track + 1

    if (n_sn > 1) {
      dat_vline <- dat |>
        dplyr::slice(seq_len(length(maxs_x))) |>
        dplyr::mutate(
          xmax = maxs_x,
          y = env$min_y,
          ymin = min_y - (.height * spacing_ratio),
          type = "vline"
        ) |>
        dplyr::slice(seq_len((dplyr::n() - 1)))

      dat <- dplyr::bind_rows(dat, dat_vline)
    }

    dat
  }
)

GeomTrack2 <- ggproto(
  "GeomTrack2",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "type", "name", "seqname"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params,
    "draw_boundary", "boundary_colour", "linetype",
    "rasterize", "dpi", "dev", "scale"
  ),
  draw_key = ggplot2::draw_key_polygon,
  default_aes = ggplot2::aes(
    fill = "black", colour = NA, alpha = 1, fontsize = 5
  ),
  draw_panel = function(
    data, panel_params, coord,
    rasterize = TRUE, dpi = 300, dev = "cairo", scale = 1,
    draw_boundary = TRUE, boundary_colour = "black", linetype = "dashed"
  ) {
    coords <- coord$transform(data, panel_params)

    coords_track <- coords |>
      dplyr::filter(type == "track")
    grob_track <- grid::polygonGrob(
      x = c(
        coords_track$x, coords_track$x, coords_track$xmax, coords_track$xmax
      ),
      y = c(
        coords_track$y, coords_track$ymin, coords_track$ymin, coords_track$y
      ),
      id = rep(seq_len(nrow(coords_track)), 4),
      gp = grid::gpar(
        fill = coords_track$fill, col = coords_track$colour
      ),
      default.units = "native"
    )

    grob_vline <- grid::nullGrob()
    if (length(unique(coords_track$seqname)) > 1 && draw_boundary) {
      coords_vline <- coords |>
        dplyr::filter(type == "vline")
      grob_vline <- grid::polylineGrob(
        x = c(coords_vline$xmax, coords_vline$xmax),
        y = c(coords_vline$ymin, coords_vline$y),
        id = rep(seq_len(nrow(coords_vline)), 2),
        gp = grid::gpar(col = boundary_colour, lty = linetype),
        default.units = "native"
      )
    }

    coords_axis <- coords |>
      dplyr::filter(type == "axis")
    grob_axis <- grid::segmentsGrob(
      x0 = coords_axis$x - (1 / 200),
      x1 = coords_axis$x - (1 / 200),
      y0 = coords_axis$ymin,
      y1 = coords_axis$y,
      gp = grid::gpar(col = "black", fill = NA),
      default.units = "native"
    )
    grob_tick <- grid::segmentsGrob(
      x0 = rep(coords_axis$x - (1 / 200), nrow(coords_axis)),
      x1 = rep(coords_axis$x - (1 / 90) - (1 / 200), nrow(coords_axis)),
      y0 = c(coords_axis$ymin, coords_axis$y),
      y1 = c(coords_axis$ymin, coords_axis$y),
      gp = grid::gpar(col = "black", fill = NA),
      default.units = "native"
    )
    grob_tick_text <- grid::textGrob(
      x = rep(coords_axis$x - (1 / 90) - (1 / 200) - (1 / 900), nrow(coords_axis) * 2),
      y = c(coords_axis$ymin, coords_axis$y),
      label = c(coords_axis$text_ymin, rep("0", nrow(coords_axis))),
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = coords_axis$fontsize),
      default.units = "native"
    )

    coords_text <- coords |>
      dplyr::filter(type == "text")
    grob_text <- grid::textGrob(
      x = coords_text$x - (1 / 200) - (1 / 180),
      y = coords_text$y,
      label = coords_text$name,
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = coords_text$fontsize),
      default.units = "native"
    )

    if (rasterize) {
      class(grob_track) <- c("rasteriser", class(grob_track))
      grob_track$dpi <- dpi
      grob_track$dev <- dev
      grob_track$scale <- scale
    }

    grid::gList(
      grob_track, grob_axis, grob_tick, grob_tick_text, grob_text, grob_vline
    )
  }
)

#' geom_track2
#'
#' @description A second version of [geom_track()].
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#'   Default is `NULL`.
#' @param width_ratio The ratio of the width of each track relative to the
#'   height of the Hi-C plot. Default is `1/20`.
#' @param spacing_ratio The ratio of the spacing between two tracks.
#'   Default is `0.5`.
#' @param data_range The range of the x axis. It can be `"auto"`, `"maximum"`,
#'   or a number (vector). Default is `"auto"`.
#' @param rasterize Whether to rasterize the plot or not. Default is `TRUE`.
#' @param dpi The resolution of the rasterised plot. Default is `300`.
#' @param dev The device to rasterise the plot. Default is `"cairo"`.
#' @param scale The scale of the rasterised plot. Default is `1`.
#' @param draw_boundary Whether to draw the boundary line or not when plotting
#'   multiple chromosomes. Default is `TRUE`.
#' @param boundary_colour The color of the boundary line. Default is `"black"`.
#' @param linetype The line type of the boundary line. Default is `"dashed"`.
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqnames
#' * start
#' * end
#' * score
#' * name
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' }
#' @export geom_track2
#' @aliases geom_track2
geom_track2 <- function(
  mapping = NULL, data = NULL, stat = StatTrack2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  width_ratio = 1 / 20, spacing_ratio = 0.5,
  data_range = c("auto", "maximum"),
  rasterize = TRUE, dpi = 300, dev = "cairo", scale = 1,
  draw_boundary = TRUE, boundary_colour = "black", linetype = "dashed"
) {
  ggplot2::layer(
    geom = GeomTrack2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio,
      data_range = data_range,
      rasterize = rasterize, dpi = dpi, dev = dev, scale = scale,
      draw_boundary = draw_boundary, boundary_colour = boundary_colour,
      linetype = linetype
    )
  )
}
