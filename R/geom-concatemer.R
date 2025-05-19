StatConcatemer <- ggplot2::ggproto(
  "StatConcatemer",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "width_ratio", "spacing_ratio"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales,
    width_ratio, spacing_ratio, concatemer_granges, concatemer_paths,
    group_identifier
  ) {
    # ======================================================================== #
    #   ^                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, min_y)                                                           #
    #   | +-----------+        +------------------------------------+          #
    #   | |           |        | (x, y)                   (xmax, y) |          #
    #   | |           |--------|                                    |          #
    #   | |           |        | (x, ymin)             (xmax, ymin) |          #
    #   | +-----------+        +------------------------------------+          #
    # ======================================================================== #

    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
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
      if (env$n_hic == 1 || env$n_track == 1 || env$n_annotation == 1) {
        maxs_x <- env$maxs_x
      } else {
        maxs_x <- dat_hic |>
          dplyr::group_by(seqnames2) |>
          dplyr::summarize(maxs_x = max(xend)) |>
          dplyr::pull(maxs_x) |>
          stats::setNames(unique(dat_hic$seqnames2))
      }
    }
    min_y <- ifelse(
      n_annotation > n_concatemer || n_track > n_concatemer, env$min_y, 0
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

    .height <- max_y * width_ratio

    if (!is.null(concatemer_granges)) {
      concatemers <- concatemer_granges[
        IRanges::overlapsAny(concatemer_granges, grs_range, ignore.strand = TRUE)
      ] |>
        GenomicRanges::sort()
    }
    if (!is.null(concatemer_paths)) {
      concatemers <- concatemer_paths |>
        rtracklayer::import(which = grs_range) |>
        GenomicRanges::sort()
    }

    ids <- mcols(concatemers)[[group_identifier]] |>
      as.character() |>
      unique()
    # list_concatemer <- purrr::map(
    #   ids, ~ concatemers[mcols(concatemers)[[group_identifier]] == .x] |>
    #     tibble::as_tibble() |>
    #     dplyr::mutate({{ group_identifier }} := as.character(.data[[group_identifier]])) |>
    #     dplyr::rename(seqname = seqnames)
    # )
    df_concatemer <- concatemers |>
      tibble::as_tibble() |>
      dplyr::mutate({{ group_identifier }} := as.character(.data[[group_identifier]])) |>
      dplyr::rename(seqname = seqnames)

    ys <- ids |>
      tibble::as_tibble() |>
      dplyr::rename(!!group_identifier := value) |>
      dplyr::mutate(
        y = min_y -
          (.height * (dplyr::row_number() - 1)) -
          (dplyr::row_number() * (.height * spacing_ratio)) -
          res / 2
      )

    dat <- df_concatemer |>
      dplyr::left_join(ys, by = group_identifier) |>
      dplyr::group_by(dplyr::across(all_of(group_identifier))) |>
      dplyr::reframe(
        seqname = seqname,
        y = first(y),
        x = start,
        xmax = end,
        ymin = y - .height,
        type = "concatemer"
      )

    if (!any(stringr::str_detect(data$seqnames1, "^chr"))) {
      dat <- dat |>
        dplyr::mutate(seqname = stringr::str_remove(seqname, "^chr"))
    }

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        calculate_add_lengths()
      chroms_sub <- data |>
        calculate_subtract_lengths()

      dat <- dat |>
        adjust_coordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))
    }
    dat_gap <- dat |>
      split(dat[[group_identifier]]) |>
      purrr::map_dfr(
        function(.x) {
          if (nrow(.x) < 2) return(NULL)
          tibble::tibble(
            y = (.x$y[1] + .x$ymin[1]) * 0.5,
            x = .x$xmax[1],
            xmax = .x$x[length(.x$x)],
            type = "gap",
            !!group_identifier := .x[[group_identifier]][1],
            seqname = "chrN",
            ymin = .x$ymin[1],
          )
        }
      )
    dat <- dplyr::bind_rows(dat, dat_gap)

    env$min_y <- min(dat$ymin)
    env$n_concatemer <- env$n_concatemer + 1

    dat
  }
)

GeomConcatemer <- ggplot2::ggproto(
  "GeomConcatemer",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "type", "seqname"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params,
    "fill"
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fill
  ) {
    coords <- coord$transform(data, panel_params)

    coords_concatemer <- coords |>
      dplyr::filter(type == "concatemer")
    grob_concatemer <- grid::polygonGrob(
      x = c(
        coords_concatemer$x, coords_concatemer$x,
        coords_concatemer$xmax, coords_concatemer$xmax
      ),
      y = c(
        coords_concatemer$y, coords_concatemer$ymin,
        coords_concatemer$ymin, coords_concatemer$y
      ),
      id = rep(seq_len(nrow(coords_concatemer)), 4),
      gp = grid::gpar(fill = fill),
      default.units = "native"
    )

    coords_gap <- coords |>
      dplyr::filter(type == "gap")
    grob_gap <- grid::polylineGrob(
      x = c(coords_gap$x, coords_gap$xmax),
      y = c(coords_gap$y, coords_gap$y),
      id = rep(seq_len(nrow(coords_gap)), 2),
      gp = grid::gpar(fill = fill, lwd = 0.3),
      default.units = "native"
    )

    grids <- grid::gList(grob_concatemer, grob_gap)

    grids
  }
)

#' geom_concatemer
#'
#' @description A ggplot2 geom for concatemers containing multi-way contacts.
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param width_ratio The ratio of the width of each concatemer track
#'   relative to the height of the Hi-C plot. Default is `1/100`.
#' @param spacing_ratio The ratio of the spacing between two tracks.
#'   Default is `1/5`.
#' @param concatemer_granges The GRanges object of the concatemer tracks.
#'   Default is `NULL`.
#' @param concatemer_paths The paths to the concatemer files. Default is `NULL`.
#' @param group_identifier A character indicating the column name in the
#'   concatemer GRanges object that identifies the group of
#'   concatemers. Default is `NULL`.
#' @param fill The fill color of the gene model track. Default is `"black"`.
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
#' }
#' @export geom_concatemer
#' @aliases geom_concatemer
geom_concatemer <- function(
  mapping = NULL, data = NULL, stat = StatConcatemer, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  width_ratio = 1 / 100, spacing_ratio = 1 / 5,
  concatemer_granges = NULL, concatemer_paths = NULL, group_identifier = NULL,
  fill = "black"
) {
  ggplot2::layer(
    geom = GeomConcatemer, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(
      na.rm = na.rm, ...,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio,
      concatemer_granges = concatemer_granges,
      concatemer_paths = concatemer_paths, group_identifier = group_identifier,
      fill = fill
    )
  )
}
