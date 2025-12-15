#' StatTad
#' @keywords internal
#' @noRd
StatTad <- ggplot2::ggproto(
  "StatTad",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, tad_path = NULL, tad_gis = NULL, is_0_based = FALSE
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
      MIN_Y <- min(dat_hic$y, na.rm = TRUE)
    }

    if (!is.null(tad_path)) {
      tmp <- read.delim(
        tad_path,
        header = FALSE, col.names = c("chrom", "start", "end")
      )

      if (is_0_based) tmp$start <- tmp$start + 1

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
      gis_data <- .tbl2Gis(data)
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
        .calculateAddLengths()
      chroms_sub <- data |>
        .calculateSubtractLengths()

      dat_tad <- dat_tad |>
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

#' GeomTad
#' @keywords internal
#' @noRd
GeomTad <- ggplot2::ggproto(
  "GeomTad",
  ggplot2::Geom,
  required_aes = c("xmin", "xmax", "xend", "ymin", "ymax", "yend"),
  default_aes = ggplot2::aes(colour = "grey"),
  draw_key = ggplot2::draw_key_path,
  draw_panel = function(data, panel_params, coord, colour = "grey") {
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

#' Visualize TAD boundaries on Hi-C heatmap
#'
#' @description
#' Adds topologically associating domain (TAD) boundary annotations to Hi-C
#' contact maps. TAD boundaries are displayed as triangular outlines marking
#' the diagonal extent of each domain.
#'
#' @inheritParams ggplot2::geom_line
#' @inheritParams geom_hic
#' @param tad_path Character. Path to BED-like file with TAD coordinates
#'   (columns: chrom, start, end). Either `tad_path` or `tad_gis` must be
#'   provided (default: NULL).
#' @param tad_gis GInteractions object containing TAD coordinates as
#'   interactions from start to end positions. Either `tad_path` or `tad_gis`
#'   must be provided (default: NULL).
#' @param is_0_based Logical. Whether input coordinates are 0-based (BED
#'   format). Set TRUE for BED files (default: FALSE).
#' @param colour Character. Color for TAD boundary lines (default: `"grey"`).
#' @param ... Additional parameters passed to layer (unused).
#'
#' @details
#' ## Required aesthetics
#' Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
#' `start2`, `end2`
#'
#' ## TAD file format
#' Tab-delimited BED-like file with columns:
#' ```
#' chrom  start     end
#' chr1   0         1000000
#' chr1   1000000   2500000
#' ```
#'
#' ## Interpretation
#' TAD boundaries are drawn as V-shaped lines on the Hi-C diagonal. Each TAD
#' is represented by two lines forming a triangle that encloses the domain.
#' Darker regions inside triangles represent increased intra-TAD interactions.
#'
#' @return A ggplot2 layer that can be added to a gghic plot.
#'
#' @seealso [geom_loop()] for chromatin loops, [gghic()] for creating Hi-C plots
#'
#' @examples
#' \dontrun{
#' # Basic usage with TAD file
#' cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()
#' gghic(cc) + geom_tad(tad_path = "tads.bed")
#'
#' # Custom color
#' gghic(cc) + geom_tad(tad_path = "tads.bed", colour = "blue")
#'
#' # Using GInteractions object
#' tads <- rtracklayer::import("tads.bed") |> as("GInteractions")
#' gghic(cc) + geom_tad(tad_gis = tads)
#'
#' # 0-based coordinates (BED format)
#' gghic(cc) + geom_tad(tad_path = "tads.bed", is_0_based = TRUE)
#'
#' # Combined with loops
#' gghic(cc) +
#'   geom_tad(tad_path = "tads.bed", colour = "grey") +
#'   geom_loop(loop_path = "loops.bedpe", colour = "red")
#' }
#' @export
#' @aliases geom_tad
geom_tad <- function(
  mapping = NULL, data = NULL, stat = StatTad, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, tad_path = NULL,
  tad_gis = NULL, is_0_based = FALSE, colour = "grey", ...
) {
  ggplot2::layer(
    geom = GeomTad, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, tad_path = tad_path, tad_gis = tad_gis,
      is_0_based = is_0_based, colour = colour, ...
    )
  )
}
