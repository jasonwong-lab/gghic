#' @rdname geom_ideogram
#' @format NULL
#' @usage NULL
#' @export
StatIdeogram <- ggplot2::ggproto(
  "StatIdeogram",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, genome = "hg19", chrom_prefix = TRUE, highlight = TRUE,
    width_ratio = 1 / 30, length_ratio = 0.8, show_coord = FALSE
  ) {
    # ======================================================================== #
    #   ^                                                                      #
    #   | +--------------------------------------------------------+           #
    #   | | (x, ymax)                                 (xmax, ymax) |           #
    #   | |                                                        |           #
    #   | | (x, y)                                       (xmax, y) |           #
    #   | +--------------------------------------------------------+           #
    #   |                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, 0)                                                               #
    # ======================================================================== #
    bands <- .retriveCytoband(data, genome, chrom_prefix)

    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic == 1) {
      max_y <- env$max_y
      max_x <- env$max_x
      min_x <- env$min_x
      res <- env$res
    } else {
      dat_hic <- data |>
        .calculateHicCoordinates()
      max_y <- max(dat_hic$ymax, na.rm = TRUE)
      max_x <- max(dat_hic$xend, na.rm = TRUE)
      min_x <- min(dat_hic$xmin, na.rm = TRUE)
      res <- data$end1[1] - data$start1[1] + 1
    }

    .height <- max_y * width_ratio
    .scale <- ((max_x - min_x) * length_ratio) / max(bands$end)

    chroms_ordered <- c(data$seqnames1, data$seqnames2) |>
      unique() |>
      sort()
    ys <- bands |>
      dplyr::distinct(seqname) |>
      dplyr::mutate(seqname = factor(seqname, levels = chroms_ordered)) |>
      dplyr::arrange(dplyr::desc(seqname)) |>
      dplyr::mutate(
        y = (.height * (dplyr::row_number() - 1)) +
          (dplyr::row_number() * (.height / 2)) +
          max_y +
          res / 2
      )

    dat_band <- bands |>
      dplyr::left_join(ys, by = "seqname") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymax = y + .height,
        type = "band",
        x_scale = x * .scale,
        xmax_scale = xmax * .scale,
        x = x_scale + min_x - dplyr::first(x_scale),
        xmax = xmax_scale + min_x - dplyr::first(x_scale)
      ) |>
      dplyr::select(x, y, xmax, ymax, type, seqname, band_fill)

    tmp <- data |>
      dplyr::group_by(seqnames1) |>
      dplyr::summarise(start1 = min(start1), end1 = max(end1)) |>
      dplyr::rename(seqname = seqnames1, start = start1, end = end1)

    dat_text <- dat_band |>
      dplyr::group_by(seqname) |>
      dplyr::filter(xmax == max(xmax)) |>
      dplyr::left_join(tmp, by = "seqname") |>
      dplyr::mutate(
        y = (ymax + y) / 2,
        x = xmax + (max(dat_band$xmax) - min(dat_band$x)) / 50,
        info = if (show_coord) {
          paste0(seqname, ":", start, "-", end)
        } else {
          seqname
        },
        type = "text",
        band_fill = "black"
      )

    dat_hl <- NULL
    if (highlight) {
      dat_hl <- dat_band |>
        dplyr::distinct(seqname, .keep_all = TRUE) |>
        dplyr::left_join(tmp, by = "seqname") |>
        dplyr::mutate(
          x = start * .scale + min_x,
          xmax = end * .scale + min_x,
          ymax = ymax + .height / 10, y = y - .height / 10,
          type = "highlight",
          band_fill = "red"
        )
    }

    dat <- dplyr::bind_rows(dat_hl, dat_band, dat_text)

    dat
  }
)

#' @rdname geom_ideogram
#' @format NULL
#' @usage NULL
#' @export
GeomIdeogram <- ggplot2::ggproto(
  "GeomIdeogram",
  ggplot2::Geom,
  required_aes = c("x", "y", "xmax", "ymax", "type", "seqname", "band_fill"),
  default_aes = ggplot2::aes(
    colour = "red", fill = "#FFE3E680", fontsize = 10
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord, fontsize = 10, colour = "red",
    fill = "#FFE3E680"
  ) {
    coords <- coord$transform(data, panel_params)

    coords_hl <- coords |>
      dplyr::filter(type == "highlight")
    grob_hl <- grid::nullGrob()
    if (nrow(coords_hl) > 0) {
      grob_hl <- grid::polygonGrob(
        x = c(coords_hl$x, coords_hl$x, coords_hl$xmax, coords_hl$xmax),
        y = c(coords_hl$ymax, coords_hl$y, coords_hl$y, coords_hl$ymax),
        id = rep(seq_len(nrow(coords_hl)), 4),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )
    }

    coords_band <- coords |>
      dplyr::filter(type == "band")
    grob_band <- grid::polygonGrob(
      x = c(coords_band$x, coords_band$x, coords_band$xmax, coords_band$xmax),
      y = c(coords_band$ymax, coords_band$y, coords_band$y, coords_band$ymax),
      id = rep(seq_len(nrow(coords_band)), 4),
      gp = grid::gpar(col = "black", fill = coords_band$band_fill),
      default.units = "native"
    )

    coords_text <- coords |>
      dplyr::filter(type == "text")
    grob_text <- grid::textGrob(
      label = coords_text$info,
      x = coords_text$x, y = coords_text$y,
      just = c("left", "center"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )

    grid::gList(grob_band, grob_hl, grob_text)
  }
)

#' geom_ideogram
#'
#' @description A ggplot2 geom for chromosome ideogram.
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param genome The genome name. Default is `"hg19"`.
#' @param chrom_prefix Whether the input data has chromosome names
#'   with prefix 'chr' or not. Default is `TRUE`.
#' @param show_coord Whether to show coordinates on the right side of the
#'   ideogram. Default is `FALSE`.
#' @param highlight Whether to highlight the boundary of the chromosome.
#'   Default is `TRUE`.
#' @param width_ratio The ratio of the width of each chromosome ideogram
#'   relative to the height of the Hi-C plot. Default is `1/30`.
#' @param length_ratio The ratio of the length of each chromosome ideogram
#'   relative to the width of the Hi-C plot. Default is `0.8`.
#' @param fontsize The font size of the chromosome name. Default is `10`.
#' @param colour The color of the chromosome boundary. Default is `"red"`.
#' @param fill The fill color of the highlighted region on the ideogram.
#'   Default is `"#FFE3E680"`.
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
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' # Add ideogram with default hg19 genome
#' gghic(cc) + geom_ideogram(genome = "hg19")
#'
#' # Highlight region with custom colors
#' gghic(cc) +
#'   geom_ideogram(
#'     genome = "hg19", highlight = TRUE, colour = "blue", fill = "#ADD8E680"
#'   )
#'
#' # Show coordinates on ideogram
#' gghic(cc) +
#'   geom_ideogram(genome = "hg19", show_coord = TRUE, fontsize = 8)
#' }
#' @export
#' @aliases geom_ideogram
geom_ideogram <- function(
  mapping = NULL, data = NULL, stat = StatIdeogram, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, genome = "hg19",
  chrom_prefix = TRUE, highlight = TRUE, show_coord = FALSE,
  width_ratio = 1 / 30, length_ratio = 0.8, fontsize = 10, colour = "red",
  fill = "#FFE3E680", ...
) {
  ggplot2::layer(
    geom = GeomIdeogram, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(
      na.rm = na.rm, genome = genome, chrom_prefix = chrom_prefix,
      highlight = highlight, show_coord = show_coord, width_ratio = width_ratio,
      length_ratio = length_ratio, fontsize = fontsize, colour = colour,
      fill = fill, ...
    )
  )
}
