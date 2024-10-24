retrive_cytoband <- function(data, genome) {
  bands_all <- ensure_cytoband(genome = genome)

  chroms <- unique(c(data$seqnames1, data$seqnames2)) |>
    as.character()
  bands <- bands_all |>
    dplyr::filter(chrom %in% chroms)

  bnames <- bands$name
  indices <- is.na(bnames)
  if (any(indices)) {
    bnames[indices] <- glue::glue("band_na_{seq_len(sum(indices))}")
  }
  if (any(bnames == "")) {
    bnames[bnames == ""] <- glue::glue("band_null_{which(bnames == '')}")
  }

  cols_all <- biovizBase::getBioColor("CYTOBAND")
  cols <- c(
    cols_all[c("gneg", "stalk", "acen")],
    gpos = unname(cols_all["gpos100"]),
    gvar = unname(cols_all["gpos100"])
  )
  gpcols <- unique(grep("gpos", bands$gieStain, value = TRUE))
  crmp <- grDevices::colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
  posCols <- stats::setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
  cols <- c(cols, posCols) |>
    as.data.frame() |>
    tibble::rownames_to_column("type") |>
    dplyr::rename(band_fill = `c(cols, posCols)`)

  dat <- bands |>
    dplyr::mutate(chromStart = chromStart + 1, name = bnames) |>
    dplyr::rename(
      type = gieStain, seqname = chrom, start = chromStart, end = chromEnd
    ) |>
    dplyr::left_join(cols, by = "type")

  dat
}

StatIdeogram <- ggplot2::ggproto(
  "StatIdeogram",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "genome", "highlight", "width_ratio", "length_ratio"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales, genome, highlight, width_ratio, length_ratio
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
    bands <- retrive_cytoband(data, genome = genome)

    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic == 1) {
      max_y <- env$max_y
      max_x <- env$max_x
      min_x <- env$min_x
      res <- env$res
    } else {
      dat_hic <- data |>
        calculate_hic_coordinates()
      max_y <- max(dat_hic$ymax, na.rm = TRUE)
      max_x <- max(dat_hic$xend, na.rm = TRUE)
      min_x <- min(dat_hic$xmin, na.rm = TRUE)
      res <- data$end1[1] - data$start1[1] + 1
    }

    .height <- max_y * width_ratio
    .scale <- ((max_x - min_x) * length_ratio) / max(bands$end)

    ys <- bands |>
      dplyr::distinct(seqname) |>
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
        info = glue::glue("{seqname}:{start}-{end}"),
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

GeomIdeogram <- ggplot2::ggproto(
  "GeomIdeogram",
  ggplot2::Geom,
  required_aes = c("x", "y", "xmax", "ymax", "type", "seqname", "band_fill"),
  extra_params = c(ggplot2::Geom$extra_params, "fontsize"),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fontsize, colour, fill
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
#' library(scales)
#'
#' cf <- HiCExperiment::CoolFile(
#'   system.file("extdata", "cooler", "chr4_11-100kb.cool", package = "gghic")
#' )
#' hic <- HiCExperiment::import(cf)
#'
#' gis <- InteractionSet::interactions(hic)
#' gis$score <- log10(gis$balanced)
#' x <- tibble::as_tibble(gis)
#' scores <- x$score[
#'   InteractionSet::pairdist(gis) != 0 &
#'     !is.na(InteractionSet::pairdist(gis) != 0)
#' ]
#' scores <- scores[!is.na(scores) & !is.infinite(scores)]
#' x$score <- scales::oob_squish(x$score, c(min(scores), max(scores)))
#'
#' p <- x |>
#'   dplyr::filter(seqnames1 == "chr11", seqnames2 == "chr11") |>
#'   ggplot2::ggplot(
#'     ggplot2::aes(
#'       seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'       seqnames2 = seqnames2, start2 = start2, end2 = end2,
#'       fill = score
#'     )
#'   ) +
#'   geom_hic() +
#'   theme_hic()
#'
#' p + geom_ideogram(
#'   genome = "hg19", highlight = FALSE, length_ratio = 0.7, fontsize = 8
#' )
#' }
#' @export geom_ideogram
#' @aliases geom_ideogram
geom_ideogram <- function(
  mapping = NULL, data = NULL, stat = StatIdeogram, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  genome = "hg19", highlight = FALSE, width_ratio = 1 / 30, length_ratio = 0.8,
  fontsize = 10, colour = "red", fill = "#FFE3E680"
) {
  ggplot2::layer(
    geom = GeomIdeogram, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, ...,
      genome = genome, highlight = highlight, width_ratio = width_ratio,
      length_ratio = length_ratio, fontsize = fontsize, colour = colour,
      fill = fill
    )
  )
}
