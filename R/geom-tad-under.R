#' StatTadUnder
#' @keywords internal
#' @noRd
StatTadUnder <- ggplot2::ggproto(
  "StatTadUnder",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, tad_path = NULL, tad_gis = NULL, is_0_based = FALSE
  ) {
    # ======================================================================== #
    #   ^                                                                      #
    #   |                                                                      #
    # --+--------------------------------------------------------------------> #
    #   | \ (xmin, ymin) / (xend, yend)                                        #
    #   |  \            /                                                      #
    #   |   \          /                                                       #
    #   |    \        /                                                        #
    #   |     \      /                                                         #
    #   |      \    /                                                          #
    #   |       \  /                                                           #
    #   |        \/ (x, y) ==> (xmax, ymax)                                    #
    # ======================================================================== #
    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    if (env$n_hic_under <= 0) {
      stop("geom_tad_under() requires a HiC under plot to be drawn first.")
    }

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

    to_keep1 <- IRanges::overlapsAny(
      InteractionSet::anchors(gis_tad)$first, gis_data
    )
    to_keep2 <- IRanges::overlapsAny(
      InteractionSet::anchors(gis_tad)$second, gis_data
    )
    gis_tad <- gis_tad[to_keep1 & to_keep2]

    dat_tad <- gis_tad |>
      tibble::as_tibble()

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- env$chroms_add
      if (is.null(chroms_add)) {
        chroms_add <- data |>
          .calculateAddLengths()
        chroms_sub <- data |>
          .calculateSubtractLengths()
      } else {
        chroms_sub <- env$chroms_sub
      }

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

    min_y <- ifelse(
      n_annotation > 0 || n_track > 0 || n_concatemer > 0, env$min_y, 0
    )

    dat <- dat_tad |>
      dplyr::mutate(
        xmin = start1,
        xmax = (start1 + end2) / 2,
        xend = end2,
        ymin = min_y - (res * 0.5),
        ymax = (xmax - start1) * -1 + min_y - (res * 0.5),
        yend = ymin
      )

    dat
  }
)

#' Visualize TAD boundaries below Hi-C heatmap
#'
#' @description
#' Adds topologically associating domain (TAD) boundary annotations below the
#' Hi-C contact map (or below the inverted Hi-C map).
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
#' ## Interpretation
#' TAD boundaries are drawn as inverted V-shaped lines below the Hi-C diagonal
#' or below other tracks. This is useful when using `geom_hic_under()` or simply
#' to display TADs in an inverted orientation.
#'
#' @return A ggplot2 layer that can be added to a gghic plot.
#'
#' @seealso [geom_tad()], [geom_hic_under()]
#'
#' @export
#' @aliases geom_tad_under
geom_tad_under <- function(
  mapping = NULL, data = NULL, stat = StatTadUnder, position = "identity",
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
