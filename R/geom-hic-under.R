StatHicUnder <- ggplot2::ggproto(
  "StatHicUnder",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(data, scales) {
    # ======================================================================== #
    #   ^                                                                      #
    #   |                                                                      #
    # --+--------------------------------------------------------------------> #
    #   |          /\ (x, y)                                                   #
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
    #   |          \/ (xmax, ymax)                                             #
    # ======================================================================== #
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    if (env$n_hic <= 0) {
      stop("geom_hic_under() requires a HiC plot to be drawn first.")
    }
    grs_range <- env$grs_range

    gis_data <- .tbl2Gis(data)
    to_keep1 <- IRanges::overlapsAny(
      InteractionSet::anchors(gis_data)$first, grs_range
    )
    to_keep2 <- IRanges::overlapsAny(
      InteractionSet::anchors(gis_data)$second, grs_range
    )
    data <- data[to_keep1 & to_keep2, ]
    dat <- data |>
      .calculateHicCoordinates(lower = TRUE)

    min_y <- ifelse(
      n_annotation > 0 || n_track > 0 || n_concatemer > 0, env$min_y, 0
    )

    res <- data$end1[1] - data$start1[1] + 1

    dat <- dat |>
      dplyr::mutate(
        y = y + min_y - res * 0.5,
        ymin = ymin + min_y - res * 0.5,
        yend = yend + min_y - res * 0.5,
        ymax = ymax + min_y - res * 0.5
      )

    env$n_hic_under <- 1
    env$min_y <- min(dat$ymax)

    dat
  }
)

#' geom_hic_under
#'
#' @description A ggplot2 layer to plot flipped Hi-C interaction data.
#' @inheritParams ggplot2::geom_polygon
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param rasterize Whether to rasterize the plot or not. Default is `TRUE`.
#' @param dpi The resolution of the rasterised plot. Default is `300`.
#' @param dev The device to rasterise the plot. Default is `"cairo"`.
#' @param scale The scale of the rasterised plot. Default is `1`.
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
#' # Load two Hi-C datasets for comparison
#' cc1 <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' # Simulate a second dataset (in practice, load a different sample)
#' cc2 <- cc1
#'
#' # Compare two Hi-C maps with different color scales
#' library(ggplot2)
#' ggplot() +
#'   geom_hic(
#'     data = scaleData(cc1, "balanced", log10),
#'     aes(
#'       seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'       seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
#'     )
#'   ) +
#'   geom_hic_under(
#'     data = scaleData(cc2, "balanced", log10),
#'     aes(
#'       seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'       seqnames2 = seqnames2, start2 = start2, end2 = end2, fill2 = score
#'     )
#'   ) |>
#'   renameGeomAes(new_aes = c("fiil" = "fill2")) +
#'   scale_fill_gradientn(
#'     aesthetics = "fill2", colors = c("white", "blue"), name = "Sample 2"
#'   ) +
#'   theme_hic()
#' }
#' @export
#' @aliases geom_hic_under
geom_hic_under <- function(
  mapping = NULL, data = NULL, stat = StatHicUnder, position = "identity",
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
