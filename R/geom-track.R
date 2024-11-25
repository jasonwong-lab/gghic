StatTrack <- ggplot2::ggproto(
  "StatTrack",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "data_paths", "data_granges", "width_ratio", "spacing_ratio"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales,
    data_paths, data_granges, width_ratio, spacing_ratio, data_range
  ) {
    # ======================================================================== #
    #       ^                                                                  #
    #       | (min_x, max_y)                                                   #
    #       |                              [HiC plot]                          #
    # ------+----------------------------------------------------------------> #
    #       | (0, min_y)                                                       #
    #       |                                                                  #
    # x --+ |      |                       (x, ymin) +----+ (xmax, ymin)       #
    #     | |     | |        /\                      |    |                    #
    #     | |    |   |      /  \                     |    |                    #
    #     | |   |     \    /    \                    |    |                    #
    # 0 --+ | --       ---       --           (x, y) +----+ (xmax, y)          #
    # ======================================================================== #
    data_range <- data_range[1]
    if (is.null(data_paths) && is.null(data_granges)) {
      stop("Either data_paths or data_granges must be provided.")
    }
    if (!is.null(data_paths) && !is.null(data_granges)) {
      stop("Only one of data_paths or data_granges can be provided.")
    }
    if (!is.null(data_paths) && is.null(names(data_paths))) {
      names(data_paths) <- glue::glue("track_{seq_along(data_paths)}")
      n_data <- length(data_paths)
    }
    if (!is.null(data_granges) && is.null(names(data_granges))) {
      names(data_granges) <- glue::glue("track_{seq_along(data_granges)}")
      n_data <- length(data_granges)
    }
    if (
      is.numeric(data_range) &&
        length(data_range) != 1 &&
        length(data_range) != n_data
    ) {
      stop("data_range must be of length 1 or equal to the number of data.")
    }

    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    if (env$n_hic == 1) {
      max_y <- env$max_y
      max_x <- env$max_x
      min_x <- env$min_x
      res <- env$res
      n_sn <- env$n_sn
    } else {
      dat_hic <- data |>
        calculate_hic_coordinates()
      max_y <- max(dat_hic$ymax, na.rm = TRUE)
      max_x <- max(dat_hic$xend, na.rm = TRUE)
      min_x <- min(dat_hic$xmin, na.rm = TRUE)
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
    }
    if (n_sn > 1) {
      if (env$n_hic == 1 || env$n_annotation == 1) {
        maxs_x <- env$maxs_x
      } else {
        maxs_x <- dat_hic |>
          dplyr::group_by(seqnames1) |>
          dplyr::summarize(maxs_x = max(xend)) |>
          dplyr::pull(maxs_x) |>
          stats::setNames(unique(dat_hic$seqnames1))
      }
    }
    min_y <- ifelse(n_annotation > n_track, env$min_y, 0)

    .height <- max_y * width_ratio

    ys <- names(data_paths) |>
      tibble::as_tibble() |>
      dplyr::rename(name = value) |>
      dplyr::mutate(
        y = min_y -
          res / 2 -
          (.height * (dplyr::row_number())) -
          ((dplyr::row_number()) * (.height * spacing_ratio))
      )

    grs_range <- data |>
      dplyr::group_by(seqnames1) |>
      dplyr::summarize(
        start = min(start1),
        end = max(end1)
      ) |>
      dplyr::reframe(range = glue::glue("{seqnames1}:{start}-{end}")) |>
      dplyr::pull(range) |>
      GenomicRanges::GRanges()

    if (!is.null(data_paths)) {
      tracks <- purrr::map(
        data_paths, function(.x) {
          tmp <- rtracklayer::import(.x, which = grs_range)
          tibble::as_tibble(tmp)
        }
      ) |>
        dplyr::bind_rows(.id = "name") |>
        dplyr::rename(seqname = seqnames)
    }
    if (!is.null(data_granges)) {
      tracks <- purrr::map(
        data_granges, function(.x) {
          tmp <- GenomicRanges::pintersect(.x, grs_range, ignore.strand = TRUE)
          tmp <- tmp[width(tmp) > 0]
          tibble::as_tibble(tmp)
        }
      ) |>
        dplyr::bind_rows(.id = "name") |>
        dplyr::rename(seqname = seqnames)
    }

    if (data_range == "auto") {
      .max <- tracks |>
        dplyr::group_by(name) |>
        dplyr::summarize(.max = round(max(score, na.rm = TRUE), digits = 2))
    }
    if (data_range == "maximum") {
      .max <- tibble::tibble(
        name = unique(tracks$name),
        .max = round(max(tracks$score, na.rm = TRUE), digits = 2)
      )
    }
    if (is.numeric(data_range)) {
      .max <- tibble::tibble(name = unique(tracks$name), .max = data_range)
    }

    dat_track <- tracks |>
      dplyr::left_join(ys, by = "name") |>
      dplyr::left_join(.max, by = "name") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymin = y + .height / .max * (score),
        type = "track"
      )

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        calculate_add_lengths()
      chroms_sub <- data |>
        calculate_subtract_lengths()

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

GeomTrack <- ggproto(
  "GeomTrack",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "type", "name", "seqname"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params,
    "fill", "fontsize", "draw_boundary", "boundary_colour", "linetype"
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fill, fontsize, draw_boundary, boundary_colour, linetype
  ) {
    coords <- coord$transform(data, panel_params)

    coords_track <- coords |>
      dplyr::filter(type == "track")

    names_track <- coords_track$name |> unique()
    if (!is.null(names(fill))) {
      col_track <- tibble::tibble(
        name = names_track,
        col_track = fill[match(names_track, names(fill))]
      )
    }
    if (is.null(names(fill))) {
      col_track <- tibble::tibble(
        name = names_track,
        col_track = fill[seq_len(length(names_track))]
      )
    }

    col_track <- col_track |>
      tidyr::replace_na(list(col_track = "black"))
    coords_track <- coords_track |>
      dplyr::left_join(col_track, by = "name")

    grob_track <- grid::polygonGrob(
      x = c(
        coords_track$x, coords_track$x, coords_track$xmax, coords_track$xmax
      ),
      y = c(
        coords_track$y, coords_track$ymin, coords_track$ymin, coords_track$y
      ),
      id = rep(seq_len(nrow(coords_track)), 4),
      gp = grid::gpar(
        fill = coords_track$col_track, col = coords_track$col_track
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
      gp = grid::gpar(col = "black", fill = "black"),
      default.units = "native"
    )
    grob_tick <- grid::segmentsGrob(
      x0 = rep(coords_axis$x - (1 / 200), nrow(coords_axis)),
      x1 = rep(coords_axis$x - (1 / 90) - (1 / 200), nrow(coords_axis)),
      y0 = c(coords_axis$ymin, coords_axis$y),
      y1 = c(coords_axis$ymin, coords_axis$y),
      gp = grid::gpar(col = "black", fill = "black"),
      default.units = "native"
    )
    grob_tick_text <- grid::textGrob(
      x = rep(coords_axis$x - (1 / 90) - (1 / 200) - (1 / 900), nrow(coords_axis) * 2),
      y = c(coords_axis$ymin, coords_axis$y),
      label = c(coords_axis$text_ymin, rep("0", nrow(coords_axis))),
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )

    coords_text <- coords |>
      dplyr::filter(type == "text")
    grob_text <- grid::textGrob(
      x = coords_text$x - (1 / 200) - (1 / 180),
      y = coords_text$y,
      label = coords_text$name,
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )

    grid::gList(
      grob_track, grob_axis, grob_tick, grob_tick_text, grob_text, grob_vline
    )
  }
)

#' geom_track
#'
#' @description A ggplot2 geom for sequencing data tracks.
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param data_paths The paths to the sequencing data files. Default is `NULL`.
#' @param data_granges The GRanges object of the sequencing data.
#'   Default is `NULL`.
#' @param width_ratio The ratio of the width of each track relative to the
#'   height of the Hi-C plot. Default is `1/20`.
#' @param spacing_ratio The ratio of the spacing between two tracks.
#'   Default is `0.5`.
#' @param data_range The range of the x axis. It can be `"auto"`, `"maximum"`,
#'   or a number (vector). Default is `"auto"`.
#' @param fill The fill color of the track. Default is `"black"`.
#' @param fontsize The font size of the track names. Default is `5`.
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
#' p <- x |>
#'   filter(
#'     center1 > 10000000 & center1 < 11000000 &
#'       center2 > 10000000 & center2 < 11000000
#'   ) |>
#'   gghic(
#'     draw_boundary = TRUE,
#'
#'     ideogram = TRUE, genome = "hg19", highlight = TRUE,
#'     ideogram_fontsize = 7, ideogram_width_ratio = 0.08,
#'
#'     annotation = TRUE, include_ncrna = FALSE, gtf_path = path_gtf,
#'     style = "arrow", maxgap = 100000, annotation_fontsize = 5,
#'     annotation_width_ratio = 0.05,
#'
#'     expand_xaxis = TRUE
#'   )
#'
#' paths_track <- glue("{dir_cache_gghic}/track{1:2}.bigWig")
#'
#' p + geom_track(
#'   data_paths = paths_track, width_ratio = 0.5,
#'   fill = c("#DC0000B2", "#00A087B2"), data_range = "auto"
#' )
#' }
#' @export geom_track
#' @aliases geom_track
geom_track <- function(
  mapping = NULL, data = NULL, stat = StatTrack, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  data_paths = NULL, data_granges = NULL, width_ratio = 1 / 20,
  spacing_ratio = 0.5, data_range = c("auto", "maximum"),
  fill = "black", fontsize = 5,
  draw_boundary = TRUE, boundary_colour = "black", linetype = "dashed"
) {
  ggplot2::layer(
    geom = GeomTrack, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      data_paths = data_paths, data_granges = data_granges,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio,
      data_range = data_range, fill = fill, fontsize = fontsize,
      draw_boundary = draw_boundary, boundary_colour = boundary_colour,
      linetype = linetype
    )
  )
}
