#' StatLoop
#' @keywords internal
#' @noRd
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

#' GeomLoop
#' @keywords internal
#' @noRd
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

#' Visualize chromatin loops on Hi-C heatmap
#'
#' @description 
#' Adds chromatin loop annotations to Hi-C contact maps. Loops can be
#' displayed as circles (points) or arcs. Automatically filters loops to
#' display only those within the plotted region.
#'
#' @inheritParams ggplot2::geom_point
#' @inheritParams geom_hic
#' @param loop_path Character. Path to loop file in BEDPE-like format with
#'   columns: chrom1, start1, end1, chrom2, start2, end2. Either `loop_path`
#'   or `loop_gis` must be provided (default: NULL).
#' @param loop_gis GInteractions object containing loop coordinates. Either
#'   `loop_path` or `loop_gis` must be provided (default: NULL).
#' @param is_0_based Logical. Whether input coordinates are 0-based (e.g., BED
#'   format). Set TRUE for BEDPE files (default: FALSE).
#' @param style Character. Visualization style:
#'   * `"circle"`: display loops as circular points (default)
#'   * `"arc"`: display loops as curved arcs connecting anchors
#' @param n_arc_points Integer. Number of points used to draw each arc when
#'   `style = "arc"`. Higher values produce smoother curves (default: 50).
#' @param colour Character. Color for loop markers or arcs (default: `"black"`).
#' @param shape Integer. Point shape when `style = "circle"` (default: 21 =
#'   filled circle).
#' @param fill Character. Fill color for points when `style = "circle"`
#'   (default: NA for transparent).
#' @param size Unit or numeric. Size of loop markers. Can be a grid unit or
#'   numeric value in mm (default: `unit(1/80, "native")`).
#' @param stroke Numeric. Line width for point borders or arc lines (default: 1).
#' @param ... Additional parameters passed to layer (unused).
#'
#' @details
#' ## Required aesthetics
#' Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
#' `start2`, `end2`
#'
#' ## Loop file format
#' Tab-delimited file with columns:
#' ```
#' chrom1  start1  end1  chrom2  start2  end2
#' chr1    1000000 1005000 chr1  2000000 2005000
#' ```
#'
#' ## Performance tips
#' For large loop files, pre-filter to region of interest before plotting.
#' Arc style is more computationally intensive than circle style.
#'
#' @return A ggplot2 layer that can be added to a gghic plot.
#'
#' @seealso [geom_tad()] for TAD boundaries, [gghic()] for creating Hi-C plots
#'
#' @examples
#' \dontrun{
#' # Basic usage with loop file
#' cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()
#' gghic(cc) + geom_loop(loop_path = "loops.bedpe")
#'
#' # Arc style visualization
#' gghic(cc) + 
#'   geom_loop(loop_path = "loops.bedpe", style = "arc", colour = "red")
#'
#' # Using GInteractions object
#' loops <- rtracklayer::import("loops.bedpe")
#' gghic(cc) + geom_loop(loop_gis = loops, colour = "blue", size = 2)
#'
#' # 0-based coordinates (BEDPE format)
#' gghic(cc) + geom_loop(loop_path = "loops.bedpe", is_0_based = TRUE)
#'
#' # Customized appearance
#' gghic(cc) +
#'   geom_loop(
#'     loop_path = "loops.bedpe",
#'     colour = "darkred",
#'     fill = "red",
#'     size = 3,
#'     stroke = 1.5,
#'     shape = 21
#'   )
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
