StatLoop <- ggplot2::ggproto(
  "StatLoop",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales,
    loop_path, loop_gis, is_0based
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
    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic == 1) {
      res <- env$res
      n_sn <- env$n_sn
    } else {
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
      env$res <- res
      env$n_sn <- n_sn
    }

    if (!is.null(loop_path)) {
      tmp <- loop_path |>
        vroom::vroom(
          delim = "\t", col_types = vroom::cols(),
          col_names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
        )

      if (is_0based) {
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
      gis_data <- tbl2gis(data)
    } else {
      gis_data <- env$gis
    }

    to_keep1 <- IRanges::overlapsAny(anchors(gis_loop)$first, gis_data)
    to_keep2 <- IRanges::overlapsAny(anchors(gis_loop)$second, gis_data)
    gis_loop <- gis_loop[to_keep1 & to_keep2]

    dat_loop <- gis_loop |>
      tibble::as_tibble()

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        calculate_add_lengths()
      chroms_sub <- data |>
        calculate_subtract_lengths()

      dat_loop <- dat_loop |>
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

    dat <- dat_loop |>
      dplyr::mutate(
        x = (start1 + end2) / 2,
        y = (x - start1)
      )

    dat
  }
)

GeomLoop <- ggproto(
  "GeomLoop",
  ggplot2::GeomPoint,
  required_aes = c("x", "y"),
  extra_params = c(
    ggplot2::Geom$extra_params,
    "colour", "shape", "fill", "size"
  ),
  draw_key = ggplot2::draw_key_point,
  default_aes = NULL,
  draw_panel = function(
    data, panel_params, coord,
    colour, shape, fill, size
  ) {
    coords <- coord$transform(data, panel_params)

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
      gp = grid::gpar(col = colour, fill = fill)
    )
  }
)

#' geom_loop
#'
#' @description A ggplot2 geom for drawing chromatin loops on the heatmap.
#' @inheritParams ggplot2::geom_point
#' @inheritParams geom_hic
#' @param loop_path A path to the loop file. Default is `NULL`.
#' @param loop_gis An InteractionSet object of loops. Default is `NULL`.
#' @param is_0based Whether the loop file is 0-based or not. Default is `FALSE`.
#' @param colour The color of the loops. Default is `grey`.
#' @param shape The shape of the loops. Default is `21`.
#' @param fill The fill color of the loops. Default is `NA`.
#' @param size The size of the loops. Default is `NULL`.
#'   If `NULL`, the size is `grid::unit(1 / 80, "native")`.
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
#' hic <- glue("{dir_cache_gghic}/chr4_11-100kb.cool") |>
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
#' p <- x |>
#'   filter(
#'     seqnames1 == "chr4", seqnames2 == "chr4",
#'     center1 > 11000000 & center1 < 21000000,
#'     center2 > 11000000 & center2 < 21000000
#'   ) |>
#'   gghic(expand_xaxis = TRUE)
#'
#' path_loop <- glue("{dir_cache_gghic}/loops-chr4_11.txt")
#'
#' p + geom_loop(loop_path = path_loop, is_0based = TRUE)
#' }
#' @export geom_loop
#' @aliases geom_loop
geom_loop <- function(
  mapping = NULL, data = NULL, stat = StatLoop, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  loop_path = NULL, loop_gis = NULL, is_0based = FALSE,
  colour = "black", shape = 21, fill = NA, size = NULL
) {
  ggplot2::layer(
    geom = GeomLoop, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      loop_path = loop_path, loop_gis = loop_gis, is_0based = is_0based,
      colour = colour, shape = shape, fill = fill, size = size
    )
  )
}
