StatTad <- ggplot2::ggproto(
  "StatTad",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "tad_path", "tad_gis", "is_0based"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales,
    tad_path, tad_gis, is_0based
  ) {
    # ======================================================================== #
    #   ^        /\ (xmax, ymax)                                               #
    #   |       /  \                                                           #
    #   |      /    \                                                          #
    #   |     /      \                                                         #
    #   |    /        \                                                        #
    #   |   /          \                                                       #
    #   |  /            \                                                      #
    #   | / (xmin, ymin) \ (xend, yend)                                        #
    # --+--------------------------------------------------------------------> #
    #   |                                                                      #
    # ======================================================================== #
    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic == 1) {
      res <- env$res
      n_sn <- env$n_sn
      MIN_Y <- env$MIN_Y
    } else {
      dat_hic <- data |>
        calculate_hic_coordinates()
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
      MIN_Y <- min(dat_hic$y, na.rm = TRUE)
    }

    if (!is.null(tad_path)) {
      tmp <- tad_path |>
        vroom::vroom(
          delim = "\t",
          col_types = vroom::cols(),
          col_names = c("chrom", "start", "end")
        )

      if (is_0based) tmp$start <- tmp$start + 1

      anchor1 <- GenomicRanges::GRanges(
        seqnames = tmp$chrom,
        ranges = IRanges::IRanges(tmp$start, width = res)
      )
      anchor2 <- GenomicRanges::GRanges(
        seqnames = tmp$chrom,
        ranges = IRanges::IRanges(tmp$end - res + 1, width = res)
      )
      gis_tad <- InteractionSet::GInteractions(anchor1, anchor2)
    }
    if (!is.null(tad_gis)) gis_tad <- tad_gis

    if (is.null(env$gis)) {
      gis_data <- tbl2gis(data)
    } else {
      gis_data <- env$gis
    }

    to_keep1 <- IRanges::overlapsAny(anchors(gis_tad)$first, gis_data)
    to_keep2 <- IRanges::overlapsAny(anchors(gis_tad)$second, gis_data)
    gis_tad <- gis_tad[to_keep1 & to_keep2]

    dat_tad <- gis_tad |>
      tibble::as_tibble()

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        calculate_add_lengths()
      chroms_sub <- data |>
        calculate_subtract_lengths()

      dat_tad <- dat_tad |>
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

    dat <- dat_tad |>
      dplyr::mutate(
        xmin = start1,
        xmax = (start1 + end2) / 2,
        xend = end2,
        ymin = MIN_Y + (res * 0.5),
        ymax = (xmax - start1),
        yend = ymin
      )

    dat
  }
)

GeomTad <- ggproto(
  "GeomTad",
  ggplot2::Geom,
  required_aes = c(
    "xmin", "xmax", "xend", "ymin", "ymax", "yend"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params,
    "colour"
  ),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(
    data, panel_params, coord,
    colour
  ) {
    coords <- coord$transform(data, panel_params)
    grob1 <- grid::pathGrob(
      x = c(coords$xmin, coords$xmax),
      y = c(coords$ymin, coords$ymax),
      id = rep(seq_len(nrow(coords)), 2),
      default.units = "native",
      gp = grid::gpar(fill = NA, col = colour)
    )
    grob2 <- grid::pathGrob(
      x = c(coords$xmax, coords$xend),
      y = c(coords$ymax, coords$yend),
      id = rep(seq_len(nrow(coords)), 2),
      default.units = "native",
      gp = grid::gpar(fill = NA, col = colour)
    )

    grid::gList(grob1, grob2)
  }
)

#' geom_tad
#'
#' @description A ggplot2 geom for drawing
#'   topologically associating domains (TADs) and compartments on the heatmap.
#' @inheritParams ggplot2::geom_line
#' @inheritParams geom_hic
#' @param tad_path A path to the TAD file. Default is `NULL`.
#' @param tad_gis An InteractionSet object of TADs. Default is `NULL`.
#' @param is_0based Whether the TAD file is 0-based or not. Default is `FALSE`.
#' @param colour The color of the TADs. Default is `grey`.
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
#' dir_cache_gghic <- user_cache_dir(appname = "gghic")
#' url_file <- paste0(
#'   "https://raw.githubusercontent.com/mhjiang97/gghic-data/refs/heads/",
#'   "master/cooler/chr4_11-100kb.cool"
#' )
#' path_file <- file.path(dir_cache_gghic, "chr4_11-100kb.cool")
#' download.file(url_file, path_file)
#'
#' hic <- path_file |>
#'   CoolFile() |>
#'   import()
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
#' url_file <- paste0(
#'   "https://raw.githubusercontent.com/mhjiang97/gghic-data/refs/heads/",
#'   "master/tad/TADs_500kb-chr4_11.tsv"
#' )
#' path_tad <- glue("{dir_cache_gghic}/TADs_500kb-chr4_11.tsv")
#' download.file(url_file, path_tad)
#'
#' p + geom_tad(tad_path = path_tad, is_0based = TRUE)
#' }
#' @export geom_tad
#' @aliases geom_tad
geom_tad <- function(
  mapping = NULL, data = NULL, stat = StatTad, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  tad_path = NULL, tad_gis = NULL, is_0based = FALSE,
  colour = "grey"
) {
  ggplot2::layer(
    geom = GeomTad, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      tad_path = tad_path, tad_gis = tad_gis, is_0based = is_0based,
      colour = colour
    )
  )
}
