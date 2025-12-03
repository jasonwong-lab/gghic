StatLoop <- ggplot2::ggproto(
  "StatLoop",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, loop_path = NULL, loop_gis = NULL, is_0_based = FALSE
  ) {
    # ======================================================================== #
    #   ^        /\                                                            #
    #   |       /  \                                                           #
    #   |      /    \                                                          #
    #   |     /      \                                                         #
    #   |    /   ()   \                                                        #
    #   |   /          \                                                       #
    #   |  /   ()       \                                                      #
    #   | /          ()  \                                                     #
    # --+--------------------------------------------------------------------> #
    #   |                                                                      #
    # ======================================================================== #
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic == 1) {
      res <- env$res
      n_sn <- env$n_sn
      MIN_Y <- env$MIN_Y
    } else {
      dat_hic <- data |>
        .calculateHicCoordinates()
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
      env$res <- res
      env$n_sn <- n_sn
      env$MIN_Y <- min(dat_hic$y, na.rm = TRUE)
      MIN_Y <- env$MIN_Y
    }

    if (!is.null(loop_path)) {
      tmp <- read.delim(
        loop_path,
        header = FALSE, col.names = c(
          "chrom1", "start1", "end1", "chrom2", "start2", "end2"
        )
      )

      if (is_0_based) {
        tmp$start1 <- tmp$start1 + 1
        tmp$start2 <- tmp$start2 + 1
      }

      anchor1 <- GenomicRanges::GRanges(
        seqnames = tmp$chrom1,
        ranges = IRanges::IRanges(tmp$start1, tmp$end1)
      )
      anchor2 <- GenomicRanges::GRanges(
        seqnames = tmp$chrom2,
        ranges = IRanges::IRanges(tmp$start2, tmp$end2)
      )
      gis_loop <- InteractionSet::GInteractions(anchor1, anchor2)
    }

    if (!is.null(loop_gis)) gis_loop <- loop_gis

    if (is.null(env$gis)) {
      gis_data <- .tbl2Gis(data)
    } else {
      gis_data <- env$gis
    }

    .first <- InteractionSet::anchors(gis_loop)$first
    .second <- InteractionSet::anchors(gis_loop)$second
    Seqinfo::seqlevels(.first) <- Seqinfo::seqlevels(gis_data)
    Seqinfo::seqlevels(.second) <- Seqinfo::seqlevels(gis_data)
    Seqinfo::seqinfo(.first) <- Seqinfo::seqinfo(gis_data)
    Seqinfo::seqinfo(.second) <- Seqinfo::seqinfo(gis_data)

    to_keep1 <- IRanges::overlapsAny(.first, gis_data)
    to_keep2 <- IRanges::overlapsAny(.second, gis_data)
    gis_loop <- gis_loop[to_keep1 & to_keep2]

    dat_loop <- gis_loop |>
      tibble::as_tibble()

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        .calculateAddLengths()
      chroms_sub <- data |>
        .calculateSubtractLengths()

      dat_loop <- dat_loop |>
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

    dat <- dat_loop |>
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

GeomLoop <- ggplot2::ggproto(
  "GeomLoop",
  ggplot2::GeomPoint,
  required_aes = c("x", "y"),
  draw_key = ggplot2::draw_key_point,
  default_aes = ggplot2::aes(
    colour = "black", shape = 21, fill = NA,
    size = grid::unit(1 / 80, "native"), stroke = 1
  ),
  draw_panel = function(
    data, panel_params, coord, colour = "black", shape = 21, fill = NA,
    size = NULL, stroke = 1, style = "circle", n_arc_points = 50
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
            col = coords$colour[i] %||% colour,
            lwd = coords$stroke[i] %||% stroke
          )
        )
      })
      do.call(grid::gList, grobs)
    } else {
      if (is.null(size)) {
        size <- grid::unit(1 / 80, "native")
      } else if (is.numeric(size)) {
        size <- grid::unit(size, "mm")
      } else if (grid::is.unit(size)) {
        size <- size
      }

      grid::pointsGrob(
        coords$x, coords$y,
        pch = shape, size = size,
        default.units = "native",
        gp = grid::gpar(col = colour, fill = fill, lwd = stroke)
      )
    }
  }
)

#' geom_loop
#'
#' @description A ggplot2 geom for drawing chromatin loops on the heatmap.
#' @inheritParams ggplot2::geom_point
#' @inheritParams geom_hic
#' @param loop_path A path to the loop file. Default is `NULL`.
#' @param loop_gis An InteractionSet object of loops. Default is `NULL`.
#' @param is_0_based Whether the loop file is 0-based or not.
#'   Default is `FALSE`.
#' @param style The style for drawing loops: `"circle"` for points/circles
#'   or `"arc"` for arcs under the Hi-C heatmap. Default is `"circle"`.
#' @param n_arc_points Number of points used to draw each arc (only used when
#'   `style = "arc"`). Default is `50`.
#' @param colour The color of the loops. Default is `"black"`.
#' @param shape The shape of the loops (only used when `style = "circle"`).
#'   Default is `21`.
#' @param fill The fill color of the loops
#'   (only used when `style = "circle"`).
#'   Default is `NA`.
#' @param size The size of the loops (only used when `style = "circle"`).
#'   Default is `grid::unit(1 / 80, "native")`.
#' @param stroke The line width of the loops. Default is `1`.
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
#' # Add loops from file
#' loop_file <- "path/to/loops.bedpe"
#' gghic(cc) + geom_loop(loop_path = loop_file)
#'
#' # Draw loops as arcs
#' gghic(cc) + geom_loop(loop_path = loop_file, style = "arc")
#' }
#' @export
#' @aliases geom_loop
geom_loop <- function(
  mapping = NULL, data = NULL, stat = StatLoop, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, loop_path = NULL,
  loop_gis = NULL, is_0_based = FALSE, style = "circle",
  n_arc_points = 50, colour = "black", shape = 21, fill = NA,
  size = grid::unit(1 / 80, "native"), stroke = 1, ...
) {
  ggplot2::layer(
    geom = GeomLoop, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, loop_path = loop_path, loop_gis = loop_gis,
      is_0_based = is_0_based, style = style,
      n_arc_points = n_arc_points, colour = colour, shape = shape, fill = fill,
      size = size, stroke = stroke, ...
    )
  )
}
