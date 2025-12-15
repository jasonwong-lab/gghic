#' StatTrack2
#' @keywords internal
#' @noRd
StatTrack2 <- ggplot2::ggproto(
  "StatTrack2",
  ggplot2::Stat,
  required_aes = c("seqnames", "start", "end", "score", "name"),
  dropped_aes = c("seqnames"),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales, width_ratio = 1 / 20, spacing_ratio = 0.5,
    data_range = c("auto", "maximum")
  ) {
    n_data <- length(unique(data$name))
    if (
      is.numeric(data_range) &&
        (length(data_range) != 1 && length(data_range) != n_data)
    ) {
      stop("data_range must be of length 1 or equal to the number of data.")
    }
    tracks <- split(data, data$name) |>
      purrr::map(
        function(x) {
          x |>
            dplyr::arrange(seqnames, start, end)
        }
      )

    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    chroms_add <- chroms_sub <- NULL
    if (env$n_hic == 1) {
      max_y <- env$max_y
      min_x <- env$min_x
      res <- env$res
      n_sn <- env$n_sn
      chroms_add <- env$chroms_add
      chroms_sub <- env$chroms_sub
    } else {
      max_y <- 10000
      res <- 0
      n_sn <- length(unique(c(data$seqnames)))

      if (n_sn > 1) {
        chroms_add <- data |>
          .calculateAddLengths1d()
        chroms_sub <- data |>
          .calculateSubtractLengths1d()
      }

      seqname_1st <- tracks[[1]]$seqnames[1]
      min_x <- min(data$start[data$seqnames == seqname_1st], na.rm = TRUE)

      if (n_sn > 1) {
        min_x <- min_x - chroms_sub$orignal_start[
          chroms_sub$seqname == seqname_1st
        ]
      }
    }
    if (n_sn > 1) {
      maxs_x <- data |>
        dplyr::group_by(seqnames) |>
        dplyr::summarize(maxs_x = max(end)) |>
        dplyr::pull(maxs_x) |>
        stats::setNames(unique(data$seqnames)) |>
        as.data.frame() |>
        setNames("x_max") |>
        tibble::rownames_to_column("seqname") |>
        .adjustCoordinates2(chroms_add, chroms_sub, c(x_max = "x_max"))
      maxs_x <- maxs_x$x_max |>
        setNames(maxs_x$seqname)
    }
    min_y <- ifelse(
      n_annotation > n_track || n_concatemer > n_track, env$min_y, 0
    )

    .height <- max_y * width_ratio

    names_track <- unique(data$name)

    ys <- names_track |>
      tibble::as_tibble() |>
      dplyr::rename(name = value) |>
      dplyr::mutate(
        y = min_y -
          res / 2 -
          (.height * (dplyr::row_number())) -
          ((dplyr::row_number()) * (.height * spacing_ratio))
      )

    if (data_range[1] == "auto") {
      .max <- data |>
        dplyr::group_by(name) |>
        dplyr::summarize(.max = round(max(score, na.rm = TRUE), digits = 2))
    }
    if (data_range[1] == "maximum") {
      .max <- tibble::tibble(
        name = unique(data$name),
        .max = round(max(data$score, na.rm = TRUE), digits = 2)
      )
    }
    if (is.numeric(data_range)) {
      .max <- tibble::tibble(name = unique(data$name), .max = data_range)
    }

    dat_track <- data |>
      dplyr::rename(seqname = seqnames) |>
      dplyr::left_join(ys, by = "name") |>
      dplyr::left_join(.max, by = "name") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymin = y + ((.height / .max) * score),
        type = "track"
      ) |>
      dplyr::mutate(
        ymin = ifelse(ymin - y > .height, y + .height, ymin)
      )

    if (!is.null(chroms_add) && !is.null(chroms_sub)) {
      dat_track <- dat_track |>
        .adjustCoordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))
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
        dplyr::slice(seq_along(maxs_x)) |>
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

#' GeomTrack2
#' @keywords internal
#' @noRd
GeomTrack2 <- ggplot2::ggproto(
  "GeomTrack2",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "type", "name", "seqname"
  ),
  draw_key = ggplot2::draw_key_polygon,
  default_aes = ggplot2::aes(
    fill = "black", colour = NA, alpha = 1, fontsize = 5
  ),
  draw_panel = function(
    data, panel_params, coord, rasterize = TRUE, dpi = 300, dev = "cairo",
    scale = 1, draw_boundary = TRUE, boundary_colour = "black",
    linetype = "dashed"
  ) {
    coords <- coord$transform(data, panel_params)

    coords_track <- coords |>
      dplyr::filter(type == "track")
    grob_track <- grid::polygonGrob(
      x = c(
        coords_track$x, coords_track$x, coords_track$xmax, coords_track$xmax
      ),
      y = c(
        coords_track$y, coords_track$ymin, coords_track$ymin, coords_track$y
      ),
      id = rep(seq_len(nrow(coords_track)), 4),
      gp = grid::gpar(
        fill = coords_track$fill, col = coords_track$colour
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
      gp = grid::gpar(col = "black", fill = NA),
      default.units = "native"
    )
    grob_tick <- grid::segmentsGrob(
      x0 = rep(coords_axis$x - (1 / 200), nrow(coords_axis)),
      x1 = rep(coords_axis$x - (1 / 90) - (1 / 200), nrow(coords_axis)),
      y0 = c(coords_axis$ymin, coords_axis$y),
      y1 = c(coords_axis$ymin, coords_axis$y),
      gp = grid::gpar(col = "black", fill = NA),
      default.units = "native"
    )
    grob_tick_text <- grid::textGrob(
      x = rep(
        coords_axis$x - (1 / 90) - (1 / 200) - (1 / 900),
        nrow(coords_axis) * 2
      ),
      y = c(coords_axis$ymin, coords_axis$y),
      label = c(coords_axis$text_ymin, rep("0", nrow(coords_axis))),
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = coords_axis$fontsize),
      default.units = "native"
    )

    coords_text <- coords |>
      dplyr::filter(type == "text")
    grob_text <- grid::textGrob(
      x = coords_text$x - (1 / 200) - (1 / 180),
      y = coords_text$y,
      label = coords_text$name,
      just = c("right", "centre"),
      gp = grid::gpar(col = "black", fontsize = coords_text$fontsize),
      default.units = "native"
    )

    if (rasterize) {
      class(grob_track) <- c("rasteriser", class(grob_track))
      grob_track$dpi <- dpi
      grob_track$dev <- dev
      grob_track$scale <- scale
    }

    grid::gList(
      grob_track, grob_axis, grob_tick, grob_tick_text, grob_text, grob_vline
    )
  }
)

#' Add genomic signal tracks with direct data input (a second version)
#'
#' @description
#' Another version of [geom_track()] that accepts pre-loaded genomic signal
#' data as a data frame or tibble with explicit aesthetic mappings. Displays
#' continuous genomic signals (ChIP-seq, ATAC-seq, RNA-seq, etc.) as filled area
#' tracks below Hi-C heatmaps. Unlike [geom_track()], this function uses
#' standard ggplot2 aesthetic mappings, making it more flexible for custom
#' styling and integration with ggplot2 scales and themes.
#'
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param width_ratio Numeric value controlling the height of each track relative
#'   to the Hi-C plot height. Smaller values create shorter tracks. Default is
#'   `1/20` (5% of Hi-C plot height).
#' @param spacing_ratio Numeric value controlling the vertical spacing between
#'   tracks as a proportion of track height. Default is `0.5` (50% of track height).
#' @param data_range Character string or numeric value(s) controlling y-axis scaling:
#'   * `"auto"`: Each track scaled independently to its own maximum value
#'   * `"maximum"`: All tracks share the same y-axis scale (global maximum)
#'   * Numeric value: Fixed maximum for all tracks (e.g., `100`)
#'   * Numeric vector: Individual maximum for each track (length must match number of tracks)
#'   Default is `"auto"`.
#' @param rasterize Logical indicating whether to rasterize the track polygons for
#'   faster rendering and smaller file sizes. Recommended for high-resolution data.
#'   Default is `TRUE`.
#' @param dpi Numeric value specifying the resolution (dots per inch) for rasterized
#'   tracks. Higher values increase quality but also file size. Default is `300`.
#' @param dev Character string specifying the graphics device for rasterization.
#'   Options include `"cairo"`, `"ragg"`, or other devices. Default is `"cairo"`.
#' @param scale Numeric scaling factor for rasterized output. Values > 1 increase
#'   resolution. Default is `1`.
#' @param draw_boundary Logical indicating whether to draw vertical boundary lines
#'   between chromosomes in multi-chromosome displays. Default is `TRUE`.
#' @param boundary_colour Character string specifying the color of chromosome
#'   boundary lines. Default is `"black"`.
#' @param linetype Character string or integer specifying the line type for
#'   chromosome boundaries. Default is `"dashed"`.
#' @param ... Additional parameters (currently ignored).
#'
#' @details
#' ## Required Aesthetics
#' This geom requires the following aesthetics to be mapped:
#' * `seqnames`: Chromosome name for each genomic bin
#' * `start`: Start position of each bin
#' * `end`: End position of each bin
#' * `score`: Signal intensity value (numeric)
#' * `name`: Track identifier/label (groups bins into tracks)
#'
#' ## Comparison with geom_track()
#' | Feature | geom_track() | geom_track2() |
#' |---------|-------------|---------------|
#' | Input | File paths or GRanges list | Data frame with aesthetics |
#' | Flexibility | Simple, automatic | Full ggplot2 integration |
#' | Color control | `fill` parameter | `aes(fill = ...)` + scales |
#' | Data loading | Automatic from files | Manual pre-loading |
#' | Best for | Quick visualization | Custom styling & scales |
#'
#' ## Input Data Format
#' Data should be a data frame or tibble where each row represents one genomic bin:
#' ```
#' seqnames  start    end  score      name
#' chr1      10000  10100   5.2    "ChIP-seq"
#' chr1      10100  10200   6.1    "ChIP-seq"
#' chr1      10200  10300   4.8    "ChIP-seq"
#' ```
#'
#' ## Track Visualization
#' Each track displays:
#' * **Signal area**: Filled polygon showing signal intensity
#' * **Y-axis**: Left side with minimum (0) and maximum value labels
#' * **Track label**: Name (from `name` aesthetic) displayed on the left
#' * **Baseline**: Zero line at bottom of each track
#'
#' ## Color and Fill Aesthetics
#' This function fully integrates with ggplot2's aesthetic system:
#' * Map `fill` aesthetic to track names for automatic coloring
#' * Use `scale_fill_manual()`, `scale_fill_brewer()`, etc. for custom colors
#' * Supports all standard ggplot2 color specifications
#' * Can combine with other aesthetic mappings (alpha, etc.)
#'
#' ## Y-Axis Scaling
#' Control track scaling via `data_range`:
#' * **Independent** (`"auto"`): Compare patterns within each track
#' * **Shared** (`"maximum"`): Compare absolute signal between tracks
#' * **Fixed** (numeric): Consistent scaling across multiple plots
#'
#' ## Performance
#' * Rasterization recommended for dense genomic data
#' * Adjust `dpi` to balance quality and file size
#' * Pre-filter data to displayed region for faster rendering
#'
#' ## Multi-Chromosome Support
#' When data includes multiple chromosomes:
#' * Coordinates automatically adjusted for proper alignment
#' * Vertical boundaries separate chromosomes (if `draw_boundary = TRUE`)
#' * Works seamlessly with Hi-C multi-chromosome displays
#'
#' @return A ggplot2 layer object that can be added to a gghic plot.
#'
#' @examples
#' \dontrun{
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' # Load and prepare track data
#' library(rtracklayer)
#' track1 <- import("path/to/H3K27ac.bw")
#' track1$name <- "H3K27ac"
#' track_df <- as.data.frame(track1)
#'
#' # Basic usage with aesthetic mappings
#' library(ggplot2)
#' gghic(cc) +
#'   geom_track2(
#'     data = track_df,
#'     aes(
#'       seqnames = seqnames,
#'       start = start,
#'       end = end,
#'       score = score,
#'       name = name
#'     )
#'   )
#'
#' # Multiple tracks with color mapping
#' track2 <- import("path/to/ATAC.bw")
#' track2$name <- "ATAC-seq"
#' tracks_df <- rbind(
#'   as.data.frame(track1),
#'   as.data.frame(track2)
#' )
#'
#' gghic(cc) +
#'   geom_track2(
#'     data = tracks_df,
#'     aes(
#'       seqnames = seqnames,
#'       start = start,
#'       end = end,
#'       score = score,
#'       name = name,
#'       fill = name  # Color by track name
#'     )
#'   ) +
#'   scale_fill_manual(values = c("H3K27ac" = "#E63946", "ATAC-seq" = "#457B9D"))
#'
#' # Using scale_fill_brewer for automatic colors
#' gghic(cc) +
#'   geom_track2(
#'     data = tracks_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name, fill = name)
#'   ) +
#'   scale_fill_brewer(palette = "Set2")
#'
#' # Shared y-axis scaling
#' gghic(cc) +
#'   geom_track2(
#'     data = tracks_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name, fill = name),
#'     data_range = "maximum"
#'   )
#'
#' # Custom track dimensions
#' gghic(cc) +
#'   geom_track2(
#'     data = track_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name),
#'     width_ratio = 1 / 15,    # Taller tracks
#'     spacing_ratio = 1         # More spacing
#'   )
#'
#' # Fixed y-axis maximum
#' gghic(cc) +
#'   geom_track2(
#'     data = track_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name),
#'     data_range = 50  # Fix max to 50
#'   )
#'
#' # Multiple tracks with different y-axis limits
#' gghic(cc) +
#'   geom_track2(
#'     data = tracks_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name),
#'     data_range = c(100, 50)  # Different max for each track
#'   )
#'
#' # High-resolution rasterization
#' gghic(cc) +
#'   geom_track2(
#'     data = track_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name),
#'     rasterize = TRUE,
#'     dpi = 600,
#'     dev = "ragg"
#'   )
#'
#' # Combine with other geoms
#' gghic(cc) +
#'   geom_track2(
#'     data = track_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name),
#'     fill = "darkblue"
#'   ) +
#'   geom_loop("path/to/loops.bedpe")
#'
#' # Multi-chromosome display
#' cc_multi <- ChromatinContacts("path/to/cooler.cool",
#'                               focus = c("chr4", "chr8")) |>
#'   import()
#' gghic(cc_multi) +
#'   geom_track2(
#'     data = tracks_df,
#'     aes(seqnames = seqnames, start = start, end = end,
#'         score = score, name = name, fill = name),
#'     draw_boundary = TRUE
#'   )
#' }
#'
#' @seealso
#' * [geom_track()] for file-based input with automatic data loading
#' * [gghic()] for creating the base Hi-C plot
#' * [geom_annotation()] for gene annotation tracks
#' * [geom_hic()] for the Hi-C heatmap layer
#'
#' @export
#' @aliases geom_track2
geom_track2 <- function(
  mapping = NULL, data = NULL, stat = StatTrack2, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, width_ratio = 1 / 20,
  spacing_ratio = 0.5, data_range = c("auto", "maximum"), rasterize = TRUE,
  dpi = 300, dev = "cairo", scale = 1, draw_boundary = TRUE,
  boundary_colour = "black", linetype = "dashed", ...
) {
  ggplot2::layer(
    geom = GeomTrack2, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE,
    params = list(
      na.rm = na.rm, width_ratio = width_ratio, spacing_ratio = spacing_ratio,
      data_range = data_range, rasterize = rasterize, dpi = dpi, dev = dev,
      scale = scale, draw_boundary = draw_boundary,
      boundary_colour = boundary_colour, linetype = linetype, ...
    )
  )
}
