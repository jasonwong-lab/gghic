StatTad2 <- ggplot2::ggproto(
  "StatTad",
  ggplot2::Stat,
  required_aes = c("seqnames", "start", "end"),
  extra_params = c(ggplot2::Stat$extra_params),
  dropped_aes = c("seqnames", "start", "end"),
  compute_panel = function(data, scales) {
    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    if (env$n_hic <= 0) {
      stop("geom_tad2() requires a HiC plot to be drawn first.")
    }
    res <- env$res
    MIN_Y <- env$MIN_Y
    chrom_add <- env$chrom_add
    chrom_sub <- env$chrom_sub

    anchor1 <- GenomicRanges::GRanges(
      seqnames = data$seqnames,
      ranges = IRanges::IRanges(data$start, width = res)
    )
    anchor2 <- GenomicRanges::GRanges(
      seqnames = data$seqnames,
      ranges = IRanges::IRanges(data$end - res + 1, width = res)
    )
    dat_tad <- InteractionSet::GInteractions(anchor1, anchor2) |>
      tibble::as_tibble()

    aes_other <- c(
      colour = "black", fill = NA, linetype = 1, alpha = 1, stroke = 0.5
    )
    for (name_aes in names(aes_other)) {
      if (name_aes %in% colnames(data)) {
        dat_tad <- dat_tad |>
          dplyr::mutate(!!name_aes := data[[name_aes]])
      }
    }

    if (!is.null(chrom_add) && !is.null(chrom_sub)) {
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

GeomTad2 <- ggproto(
  "GeomTad",
  ggplot2::Geom,
  required_aes = c("xmin", "xmax", "xend", "ymin", "ymax", "yend"),
  extra_params = c(ggplot2::Geom$extra_params),
  draw_key = ggplot2::draw_key_path,
  default_aes = ggplot2::aes(
    colour = "black", fill = NA, linetype = 1, alpha = 1, stroke = 0.5
  ),
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)

    grob <- grid::polylineGrob(
      x = c(coords$xmin, coords$xmax, coords$xend),
      y = c(coords$ymin, coords$ymax, coords$yend),
      id = rep(seq_len(nrow(coords)), 3),
      default.units = "native",
      gp = grid::gpar(
        col = coords$colour, fill = coords$fill, lwd = coords$stroke,
        alpha = coords$alpha
      )
    )

    grid::gList(grob)
  }
)

#' geom_tad2
#'
#' @description A second version of [geom_tad()].
#' @inheritParams ggplot2::geom_line
#' @inheritParams geom_hic
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqname
#' * start
#' * end
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' # example usage
#' }
#' @export geom_tad2
#' @aliases geom_tad2
geom_tad2 <- function(
  mapping = NULL, data = NULL, stat = StatTad2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, ...
) {
  ggplot2::layer(
    geom = GeomTad2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(na.rm = na.rm, ...)
  )
}
