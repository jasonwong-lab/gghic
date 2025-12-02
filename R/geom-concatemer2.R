StatConcatemer2 <- ggplot2::ggproto(
  "StatConcatemer2",
  ggplot2::Stat,
  required_aes = c("seqnames", "start", "end", "read_group"),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, check_concatemers = TRUE, width_ratio = 1 / 100,
    spacing_ratio = 1 / 5
  ) {
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    if (env$n_hic <= 0) {
      stop("geom_concatemer2() requires a HiC plot to be drawn first.")
    }
    max_y <- env$max_y
    max_x <- env$max_x
    min_x <- env$min_x
    res <- env$res
    n_sn <- env$n_sn
    chroms_add <- env$chroms_add
    chroms_sub <- env$chroms_sub
    if (n_sn > 1) {
      if (env$n_hic == 1 || env$n_track == 1 || env$n_annotation == 1) {
        maxs_x <- env$maxs_x
      } else {
        stop("geom_concatemer2() requires a HiC plot to be drawn first.")
      }
    }
    min_y <- ifelse(
      n_annotation > n_concatemer || n_track > n_concatemer, env$min_y, 0
    )

    .height <- max_y * width_ratio

    data <- data |>
      dplyr::arrange(seqnames, start, end)

    if (check_concatemers) {
      grs_range <- env$grs_range

      grs_concatemers <- data |>
        dplyr::reframe(range = paste0(seqnames, ":", start, "-", end)) |>
        dplyr::pull(range) |>
        GenomicRanges::GRanges()
      data <- data[
        suppressWarnings(
          IRanges::overlapsAny(grs_concatemers, grs_range, ignore.strand = TRUE)
        ),
      ]
      if (nrow(data) == 0) {
        return(NULL)
      }
    }

    ids <- data$read_group |>
      as.character() |>
      unique()

    data <- data |>
      dplyr::mutate(read_group = as.character(read_group)) |>
      dplyr::rename(seqname = seqnames)

    ys <- ids |>
      tibble::as_tibble() |>
      dplyr::rename(read_group = value) |>
      dplyr::mutate(
        y = min_y -
          (.height * (dplyr::row_number() - 1)) -
          (dplyr::row_number() * (.height * spacing_ratio)) -
          res / 2
      )

    data <- data |>
      dplyr::left_join(ys, by = "read_group") |>
      dplyr::group_by(read_group) |>
      dplyr::mutate(
        seqname = seqname,
        y = first(y),
        x = start,
        xmax = end,
        ymin = y - .height,
        type = "concatemer"
      )

    if (!env$has_chr) {
      data <- data |>
        dplyr::mutate(seqname = stringr::str_remove(seqname, "^chr"))
    }

    if (!is.null(chroms_add) && !is.null(chroms_sub)) {
      data <- data |>
        .adjustCoordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))
    }

    data_gap <- data |>
      split(data$read_group) |>
      purrr::map_dfr(
        function(.x) {
          if (nrow(.x) < 2) {
            return(NULL)
          }
          tmp <- tibble::tibble(
            y = (.x$y[1] + .x$ymin[1]) * 0.5,
            x = min(.x$xmax),
            xmax = max(.x$x),
            type = "gap",
            read_group = .x$read_group[1],
            seqname = "chrN",
            ymin = .x$ymin[1],
            PANEL = .x$PANEL[1],
            group = .x$group[1]
          )
          aes_other <- c(
            colour = "black", fill = NA, linetype = 1, alpha = 1, stroke = 0.5
          )
          for (name_aes in names(aes_other)) {
            if (name_aes %in% colnames(.x)) {
              tmp <- tmp |>
                dplyr::mutate(!!name_aes := .x[[name_aes]][1])
            }
          }
          tmp
        }
      )

    dat <- dplyr::bind_rows(data, data_gap)

    env$min_y <- min(dat$ymin)
    env$n_concatemer <- env$n_concatemer + 1

    dat
  }
)

GeomConcatemer2 <- ggplot2::ggproto(
  "GeomConcatemer2",
  ggplot2::Geom,
  required_aes = c("x", "y", "xmax", "ymin", "type", "seqname"),
  draw_key = ggplot2::draw_key_polygon,
  default_aes = ggplot2::aes(
    colour = "black", fill = NA, linetype = 1, alpha = 1, stroke = 0.5
  ),
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)

    if (is.null(data)) {
      return(grid::nullGrob())
    }

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
      gp = grid::gpar(
        col = coords_concatemer$colour, fill = coords_concatemer$fill,
        lwd = coords_concatemer$stroke, alpha = coords_concatemer$alpha
      ),
      default.units = "native"
    )

    coords_gap <- coords |>
      dplyr::filter(type == "gap")
    grob_gap <- grid::polylineGrob(
      x = c(coords_gap$x, coords_gap$xmax),
      y = c(coords_gap$y, coords_gap$y),
      id = rep(seq_len(nrow(coords_gap)), 2),
      gp = grid::gpar(
        col = coords_gap$colour, fill = coords_gap$fill,
        lwd = coords_gap$stroke, alpha = coords_gap$alpha
      ),
      default.units = "native"
    )

    grids <- grid::gList(grob_concatemer, grob_gap)

    grids
  }
)

#' geom_concatemer2
#'
#' @description A second version of [geom_concatemer()].
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param check_concatemers Whether to subset the concatemers according to the
#'   main heatmap's range. Default is `TRUE`.
#' @param width_ratio The ratio of the width of each concatemer track
#'   relative to the height of the Hi-C plot. Default is `1/100`.
#' @param spacing_ratio The ratio of the spacing between two tracks.
#'   Default is `1/5`.
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqnames
#' * start
#' * end
#' * read_group
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' # Load concatemer data
#' concat_df <- rtracklayer::import("path/to/concatemers.bed") |>
#'   as.data.frame()
#'
#' concat_df <- as.data.frame(concatemers)
#'
#' # Custom styling
#' gghic(cc) +
#'   geom_concatemer2(
#'     data = concat_df,
#'     aes(
#'       seqnames = seqnames,
#'       start = start,
#'       end = end,
#'       read_group = read_name,
#'       fill = read_name
#'     ),
#'     width_ratio = 1 / 80,
#'     spacing_ratio = 1 / 3
#'   )
#' }
#' @export
#' @aliases geom_concatemer2
geom_concatemer2 <- function(
  mapping = NULL, data = NULL, stat = StatConcatemer2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = FALSE,
  check_concatemers = TRUE, width_ratio = 1 / 100, spacing_ratio = 1 / 5, ...
) {
  ggplot2::layer(
    geom = GeomConcatemer2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = TRUE, params = list(
      na.rm = na.rm, check_concatemers = check_concatemers,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio, ...
    )
  )
}
