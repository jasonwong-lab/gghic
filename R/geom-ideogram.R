#' StatIdeogram
#' @keywords internal
#' @noRd
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

#' GeomIdeogram
#' @keywords internal
#' @noRd
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

#' Add chromosome ideogram annotation to Hi-C plot
#'
#' @description
#' Displays chromosome ideograms (cytogenetic band representations) above the
#' Hi-C contact map, showing the genomic context with Giemsa staining patterns
#' and highlighting the displayed region. Ideograms provide visual reference for
#' chromosome structure, centromeres, and cytogenetic landmarks.
#'
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param genome Character. Genome assembly version for retrieving cytogenetic
#'   band information from UCSC. Supported genomes include:
#'   * `"hg38"`: Human (GRCh38/hg38)
#'   * `"hg19"`: Human (GRCh37/hg19) (default)
#'   * `"mm10"`: Mouse (GRCm38)
#'   * `"mm39"`: Mouse (GRCm39)
#'   * Other UCSC genome assemblies with cytoBand tables
#' @param chrom_prefix Logical. Whether chromosome names in the data include
#'   "chr" prefix (e.g., "chr1" vs "1"). Set FALSE for Ensembl-style naming
#'   (default: TRUE).
#' @param show_coord Logical. Display genomic coordinates (start-end) next to
#'   chromosome names on the ideogram (default: FALSE). Useful for showing exact
#'   region boundaries.
#' @param highlight Logical. Draw a colored outline around the displayed genomic
#'   region on the ideogram to emphasize the Hi-C map extent (default: TRUE).
#' @param width_ratio Numeric. Height of ideogram relative to Hi-C plot height
#'   (default: 1/30). Larger values create taller ideograms.
#' @param length_ratio Numeric. Fraction of Hi-C plot width used for ideogram
#'   length (default: 0.8 = 80%). Controls horizontal scaling to leave space for
#'   labels.
#' @param fontsize Numeric. Font size in points for chromosome labels and
#'   coordinates (default: 10).
#' @param colour Character. Border color for the highlighted region box
#'   (default: `"red"`).
#' @param fill Character. Fill color for the highlighted region box, using RGBA
#'   hex format for transparency (default: `"#FFE3E680"` = semi-transparent
#'   light red).
#' @param ... Additional parameters (unused).
#'
#' @details
#' ## Required aesthetics
#' Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
#' `start2`, `end2`
#'
#' ## Ideogram structure
#' Chromosomes are displayed horizontally above the Hi-C map with:
#' * **Cytogenetic bands**: Giemsa staining patterns (G-bands) shown with
#'   standard colors (light/dark representing staining intensity)
#' * **Centromeres**: Typically appear as darker bands near the middle
#' * **Highlighted region**: Colored box showing the exact genomic region
#'   displayed in the Hi-C map below
#' * **Labels**: Chromosome names or coordinates on the right side
#'
#' ## Multi-chromosome display
#' When visualizing multiple chromosomes, ideograms are stacked vertically in
#' the same order as they appear in the Hi-C plot.
#'
#' ## Genome assembly selection
#' Ensure the `genome` parameter matches your data's assembly. Mismatched
#' assemblies will show incorrect cytogenetic band patterns or fail to retrieve
#' band data.
#'
#' ## Performance considerations
#' Ideogram data is fetched from UCSC Genome Browser on first use per session.
#' Subsequent calls use cached data for improved performance.
#'
#' @return A ggplot2 layer that can be added to a Hi-C plot.
#'
#' @seealso [gghic()], [geom_annotation()], [geom_hic()]
#'
#' @examples
#' \dontrun{
#' # Basic usage with human genome (hg19)
#' cc <- ChromatinContacts("sample.cool", focus = "chr4") |> import()
#' gghic(cc) + geom_ideogram(genome = "hg19")
#'
#' # Use human GRCh38/hg38 assembly
#' gghic(cc) + geom_ideogram(genome = "hg38")
#'
#' # Mouse genome
#' cc_mouse <- ChromatinContacts("mouse.cool", focus = "chr1") |> import()
#' gghic(cc_mouse) + geom_ideogram(genome = "mm10")
#'
#' # Show coordinates instead of just chromosome names
#' gghic(cc) +
#'   geom_ideogram(genome = "hg19", show_coord = TRUE)
#'
#' # Customize highlight colors
#' gghic(cc) +
#'   geom_ideogram(
#'     genome = "hg19",
#'     highlight = TRUE,
#'     colour = "blue",
#'     fill = "#ADD8E680"  # Semi-transparent light blue
#'   )
#'
#' # Adjust ideogram size
#' gghic(cc) +
#'   geom_ideogram(
#'     genome = "hg19",
#'     width_ratio = 1/20,    # Taller ideogram
#'     length_ratio = 0.9,    # Wider ideogram
#'     fontsize = 12          # Larger labels
#'   )
#'
#' # Multiple chromosomes with ideograms
#' cc_multi <- ChromatinContacts("sample.cool", focus = "chr1|chr2") |>
#'   import()
#' gghic(cc_multi) +
#'   geom_ideogram(genome = "hg19", show_coord = TRUE)
#'
#' # Ensembl-style chromosome names (without "chr" prefix)
#' cc_ensembl <- ChromatinContacts("ensembl.cool", focus = "1") |> import()
#' gghic(cc_ensembl) +
#'   geom_ideogram(genome = "hg19", chrom_prefix = FALSE)
#'
#' # Disable highlighting for cleaner look
#' gghic(cc) +
#'   geom_ideogram(genome = "hg19", highlight = FALSE)
#'
#' # Complete publication-ready plot
#' gghic(cc, ideogram = FALSE) +  # Disable auto-ideogram
#'   geom_ideogram(
#'     genome = "hg19",
#'     highlight = TRUE,
#'     colour = "darkred",
#'     fill = "#FFE3E650",
#'     fontsize = 11
#'   ) +
#'   theme_hic()
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
