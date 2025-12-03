#' Scale Hi-C interaction data
#' @description
#' Transforms and scales chromatin interaction data for visualization. Applies
#' a scaling function to a specified column and handles missing values.
#' @param data A `ChromatinContacts`, `GInteractions`, `data.frame`, or `tibble`
#'   object containing chromatin interaction data.
#' @param scale_column Character string. Name of the column to scale (e.g.,
#'   `"balanced"`, `"count"`).
#' @param scale_method Function to apply for scaling. Common choices include
#'   `log10`, `log2`, or identity function `function(x) x`.
#' @param remove_na Logical. If `TRUE`, removes rows with `NA` or infinite
#'   values in the score column. Default is `FALSE`.
#' @return A tibble with columns: `seqnames1`, `start1`, `end1`, `seqnames2`,
#'   `start2`, `end2`, and `score` (the scaled values).
#' @details
#' The function:
#' 1. Converts input data to a tibble format
#' 2. Applies the scaling method to the specified column
#' 3. Creates a `score` column with the transformed values
#' 4. Optionally removes missing values
#' 5. Applies out-of-bounds squishing to ensure values are within range
#' @examples
#' \dontrun{
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool") |>
#'   import()
#'
#' # Scale using log10 transformation (most common)
#' scaled_data <- scaleData(cc, "balanced", log10)
#' head(scaled_data)
#'
#' # Use raw counts without transformation
#' scaled_raw <- scaleData(cc, "count", function(x) x)
#'
#' # Scale with log2 and remove missing values
#' scaled_clean <- scaleData(cc, "balanced", log2, remove_na = TRUE)
#'
#' # From GInteractions object
#' gis <- interactions(cc)
#' scaled_gis <- scaleData(gis, "balanced", log10)
#'
#' # From data frame (must have required columns)
#' df <- as.data.frame(gis)
#' scaled_df <- scaleData(df, "count", log10)
#'
#' # Percentile rank transformation
#' percentile_rank <- function(x) {
#'   rank(x, na.last = "keep") / sum(!is.na(x))
#' }
#' scaled_pct <- scaleData(cc, "balanced", percentile_rank)
#'
#' # Use with ggplot2 directly
#' library(ggplot2)
#' ggplot(scaled_data, aes(x = (start1 + end1) / 2, y = score)) +
#'   geom_point(alpha = 0.1) +
#'   labs(title = "Distance decay", x = "Position", y = "Log10(balanced)")
#' }
#' @export
scaleData <- function(data, scale_column, scale_method, remove_na = FALSE) {
  if (methods::is(data, "ChromatinContacts")) {
    x <- data |>
      as_tibble()
  }
  if (methods::is(data, "GInteractions")) {
    x <- data |>
      InteractionSet::swapAnchors(mode = "order") |>
      tibble::as_tibble()
  }
  if (tibble::is_tibble(data) || methods::is(data, "data.frame")) {
    x <- data
  }

  x <- x |>
    dplyr::mutate(score = scale_method(.data[[scale_column]])) |>
    dplyr::select(seqnames1, start1, end1, seqnames2, start2, end2, score)

  if (remove_na) {
    x <- x |>
      # dplyr::filter(
      #   InteractionSet::pairdist(gis) != 0,
      #   !is.na(InteractionSet::pairdist(gis) != 0)
      # ) |>
      dplyr::filter(!is.na(score), !is.infinite(score))
  }

  x <- x |>
    dplyr::mutate(
      score = scales::oob_squish(
        score, c(min(score, na.rm = TRUE), max(score, na.rm = TRUE))
      )
    )

  x
}

.checkDataType <- function(data, ...) {
  name_pkg <- .getPkgName()
  env <- get(".env", envir = asNamespace(name_pkg))

  if (methods::is(data, "GInteractions")) {
    x <- scaleData(data, ...)
    env$gis <- data
  } else if (methods::is(data, "ChromatinContacts")) {
    gis <- interactions(data)
    x <- scaleData(gis, ...)
    env$gis <- interactions(data)
  } else if (methods::is(data, "HiCExperiment")) {
    gis <- InteractionSet::interactions(data)
    x <- scaleData(gis, ...)
    env$gis <- InteractionSet::interactions(data)
  } else if (methods::is(data, "data.frame")) {
    cols_required <- c(
      "seqnames1", "seqnames2", "start1", "end1", "start2", "end2"
    )
    cols_missing <- setdiff(cols_required, colnames(data))
    if (length(cols_missing) > 0) {
      stop(
        "data must have the following columns: ",
        paste(cols_missing, collapse = ", ")
      )
    }
    x <- scaleData(data, ...)
  } else {
    stop("data must be a HiCExperiment object or a tibble/data.frame")
  }

  x
}

.calculateXRange <- function(data) {
  n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))

  if (n_sn == 1 || (n_sn == 2 && all(data$seqnames1 != data$seqnames2))) {
    data <- data
  } else {
    data <- data |>
      .adjustCoordinates(list(c(start1 = "start1", end1 = "end1")), FALSE)
  }

  c(min(data$start1), max(data$end1))
}

#' Create Hi-C visualization plot
#' @name gghic
#' @aliases gghic,ChromatinContacts-method
#' @description
#' High-level wrapper function to create publication-ready Hi-C contact map
#' visualizations with optional genomic features. Automatically handles data
#' transformation, feature integration, and theme application.
#' @param x A `ChromatinContacts` object with imported interaction data, or a
#'   `GInteractions` object, or a `data.frame`/`tibble` with interaction data.
#' @param scale_method Function to apply for data transformation. Common
#'   choices: `log10` (default), `log2`, `function(x) x` (no transformation).
#' @param scale_column Character string. Name of the column to use for scaling
#'   when input is `GInteractions` or `data.frame`. Default is `"balanced"`.
#' @param ideogram Logical. Add chromosome ideogram track. Default is `FALSE`.
#' @param ideogram_width_ratio Numeric. Height of ideogram relative to heatmap
#'   height. Default is `1/30`.
#' @param ideogram_fontsize Numeric. Font size for ideogram labels. Default is
#'   `10`.
#' @param ideogram_colour Character. Color for highlighted region on ideogram.
#'   Default is `"red"`.
#' @param ideogram_fill Character. Fill color for highlighted region. Default
#'   is `"#FFE3E680"` (transparent red).
#' @param annotation Logical. Add gene annotation track. Requires `gtf_path` in
#'   `...`. Default is `FALSE`.
#' @param annotation_style Character. Style for gene annotation: `"basic"` or
#'   `"arrow"`. Default is `"basic"`.
#' @param annotation_width_ratio Numeric. Height of annotation track relative
#'   to heatmap. Default is `1/50`.
#' @param annotation_spacing_ratio Numeric. Spacing between genes. Default is
#'   `1/3`.
#' @param annotation_fontsize Numeric. Font size for gene labels. Default is
#'   `10`.
#' @param annotation_colour Character. Color for gene features. Default is
#'   `"#48CFCB"`.
#' @param annotation_fill Character. Fill color for gene features. Default is
#'   `"#48CFCB"`.
#' @param track Logical. Add genomic signal tracks (e.g., ChIP-seq). Requires
#'   `tracks` in the `ChromatinContacts` object or `data_paths` in `...`.
#'   Default is `FALSE`.
#' @param track_width_ratio Numeric. Height of track area relative to heatmap.
#'   Default is `1/20`.
#' @param track_spacing_ratio Numeric. Spacing between multiple tracks. Default
#'   is `1/2`.
#' @param track_fill Character or vector. Colors for track signals. Default is
#'   `"black"`.
#' @param track_fontsize Numeric. Font size for track labels. Default is `5`.
#' @param tad Logical. Add TAD (topologically associating domain) boundaries.
#'   Requires `TADs` in the `ChromatinContacts` object or `tad_path` in `...`.
#'   Default is `FALSE`.
#' @param tad_colour Character. Color for TAD boundaries. Default is `"grey"`.
#' @param loop Logical. Add chromatin loop arcs. Requires `loops` in the
#'   `ChromatinContacts` object or `loop_path` in `...`. Default is `FALSE`.
#' @param loop_style Character. Style for loops: `"circle"` or `"arc"`.
#'   Default is `"circle"`.
#' @param loop_colour Character. Color for loop arcs. Default is `"black"`.
#' @param loop_fill Fill color for loop arcs. Default is `NA`.
#' @param concatemer Logical. Add concatemer visualization for multi-way
#'   contacts. Requires `multi_contacts` in object. Default is `FALSE`.
#' @param concatemer_width_ratio Numeric. Height of concatemer track. Default
#'   is `1/100`.
#' @param concatemer_spacing_ratio Numeric. Spacing between concatemers.
#'   Default is `1/5`.
#' @param expand_xaxis Logical. Expand x-axis to show full chromosome context.
#'   Default is `FALSE`.
#' @param expand_left Numeric. Left expansion in base pairs.
#'   If `NULL` (default),
#'   uses 10× resolution.
#' @param expand_right Numeric. Right expansion in base pairs. If `NULL`
#'   (default), uses 10× resolution.
#' @param ... Additional arguments passed to individual geom functions:
#'   * For `geom_ideogram()`: `genome`, `highlight`, `length_ratio`
#'   * For `geom_annotation()`: `gtf_path`, `style`, `maxgap`, `gene_symbols`,
#'     `include_ncrna`
#'   * For `geom_track()`: `data_paths`, `data_range`, `rasterize`
#'   * For `geom_tad()`: `tad_path` (path to BED file), `tad_is_0_based`
#'     (logical, default TRUE), `stroke`
#'   * For `geom_loop()`: `loop_path` (path to BEDPE file), `loop_is_0_based`
#'     (logical, default TRUE), `stroke`
#'   * For `geom_hic()`: `draw_boundary`, `rasterize`
#'
#' @return A `ggplot2` object that can be further customized with additional
#'   ggplot2 layers and functions.
#'
#' @details
#' `gghic()` provides a high-level interface for creating Hi-C visualizations.
#' It automatically:
#'
#' * Scales interaction data using the specified method
#' * Applies appropriate coordinate transformations
#' * Integrates genomic features from the `ChromatinContacts` object
#' * Sets up axis labels and breaks
#' * Applies a clean default theme
#'
#' **For more control**, use individual `geom_*()` functions with `ggplot2`:
#'
#' * `geom_hic()` - Base heatmap layer
#' * `geom_ideogram()` - Chromosome ideogram
#' * `geom_annotation()` - Gene annotations
#' * `geom_track()` - Genomic signal tracks
#' * `geom_tad()` - TAD boundaries
#' * `geom_loop()` - Chromatin loops
#' * `geom_concatemer()` - Multi-way contacts
#'
#' @examples
#' \dontrun{
#' # === Basic Usage ===
#'
#' # Load and visualize Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool") |>
#'   import()
#' gghic(cc)
#'
#' # Focus on specific region
#' cc["chr4:0-50000000"] |>
#'   gghic()
#'
#' # === With Ideogram ===
#'
#' # Add chromosome ideogram for context
#' gghic(cc, ideogram = TRUE, genome = "hg19")
#'
#' # Customize ideogram appearance
#' gghic(cc,
#'   ideogram = TRUE,
#'   ideogram_width_ratio = 1 / 25,
#'   ideogram_colour = "blue",
#'   ideogram_fill = "#ADD8E680"
#' )
#'
#' # === Adding Genomic Features ===
#'
#' # Add TAD boundaries
#' tad_file <- "path/to/tads.bed"
#' gghic(cc, tad = TRUE, tad_path = tad_file, tad_colour = "darkgreen")
#'
#' # Add chromatin loops
#' loop_file <- "path/to/loops.bedpe"
#' gghic(cc, loop = TRUE, loop_path = loop_file, loop_colour = "red")
#'
#' # Add signal tracks (e.g., ChIP-seq)
#' track1 <- "path/to/track1.bw"
#' track2 <- "path/to/track2.bw"
#' tracks <- GRangesList(
#'   H3K27ac_rep1 = rtracklayer::import(track1),
#'   H3K27ac_rep2 = rtracklayer::import(track2)
#' )
#' features(cc, "tracks") <- tracks
#'
#' gghic(cc, track = TRUE, track_fill = c("blue", "red"))
#'
#' # === Axis Expansion ===
#'
#' # Expand x-axis to show genomic context
#' gghic(
#'   cc,
#'   expand_xaxis = TRUE,
#'   expand_left = 1e6, # Extend 1 Mb to the left
#'   expand_right = 1e6 # Extend 1 Mb to the right
#' )
#' }
#' @seealso
#' * [gghic::ChromatinContacts] - Main data class
#' * [gghic::geom_hic()] - Base heatmap layer
#' * [gghic::theme_hic()] - Default theme
#' * [gghic::scaleData()] - Data transformation
#' @export
methods::setMethod(
  "gghic", "ChromatinContacts",
  function(
    x,
    scale_method = log10,
    ideogram = FALSE,
    ideogram_width_ratio = 1 / 30,
    ideogram_fontsize = 10,
    ideogram_colour = "red",
    ideogram_fill = "#FFE3E680",
    annotation = FALSE,
    annotation_style = "basic",
    annotation_width_ratio = 1 / 50,
    annotation_spacing_ratio = 1 / 3,
    annotation_fontsize = 10,
    annotation_colour = "#48CFCB",
    annotation_fill = "#48CFCB",
    track = FALSE,
    track_width_ratio = 1 / 20,
    track_spacing_ratio = 1 / 2,
    track_fill = "black",
    track_fontsize = 5,
    tad = FALSE,
    tad_colour = "grey",
    loop = FALSE,
    loop_style = "circle",
    loop_colour = "black",
    loop_fill = NA,
    concatemer = FALSE,
    concatemer_width_ratio = 1 / 100,
    concatemer_spacing_ratio = 1 / 5,
    expand_xaxis = FALSE,
    expand_left = NULL,
    expand_right = NULL,
    ...
  ) {
    dat <- x |>
      interactions() |>
      .checkDataType(scale_column = "balanced", scale_method = scale_method)

    p <- dat |>
      ggplot2::ggplot(
        ggplot2::aes(
          seqnames1 = seqnames1, start1 = start1, end1 = end1,
          seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
        )
      )

    p <- p + geom_hic(...)

    if (ideogram) {
      p <- p + geom_ideogram(
        width_ratio = ideogram_width_ratio, fontsize = ideogram_fontsize,
        colour = ideogram_colour, fill = ideogram_fill, ...
      )
    }

    if (annotation) {
      p <- p + geom_annotation(
        width_ratio = annotation_width_ratio,
        spacing_ratio = annotation_spacing_ratio, style = annotation_style,
        fontsize = annotation_fontsize, colour = annotation_colour,
        fill = annotation_fill, ...
      )
    }

    if (track && length(x@tracks) > 0) {
      p <- p + geom_track(
        track_grs = x@tracks, width_ratio = track_width_ratio,
        spacing_ratio = track_spacing_ratio, fill = track_fill,
        fontsize = track_fontsize, ...
      )
    }

    if (tad && length(x@TADs) > 0) {
      p <- p + geom_tad(
        tad_gis = .grs2Gis(x@TADs, x@resolution), colour = tad_colour, ...
      )
    }

    if (loop && length(x@loops) > 0) {
      p <- p + geom_loop(
        loop_gis = x@loops, colour = loop_colour, fill = loop_fill,
        style = loop_style, ...
      )
    }

    if (concatemer && length(x@multi_contacts) > 0) {
      p <- p + geom_concatemer(
        concatemer_granges = x@multi_contacts,
        width_ratio = concatemer_width_ratio,
        spacing_ratio = concatemer_spacing_ratio, ...
      )
    }

    breaks_labels <- .getBreaksLabels(dat)
    breaks <- breaks_labels$breaks
    labels <- breaks_labels$labels

    range_x <- NULL
    if (expand_xaxis) {
      range_x <- .calculateXRange(dat)

      if (is.null(expand_left)) expand_left <- resolution(x) * 10
      if (is.null(expand_right)) expand_right <- resolution(x) * 10
    }

    p <- p +
      theme_hic(
        breaks = breaks, labels = labels, xmin = range_x[1], xmax = range_x[2],
        expand_x = c(expand_left, expand_right)
      )

    p
  }
)

.function <- function(
  x,
  scale_column = "balanced",
  scale_method = log10,
  ideogram = FALSE,
  ideogram_width_ratio = 1 / 30,
  ideogram_fontsize = 10,
  ideogram_colour = "red",
  ideogram_fill = "#FFE3E680",
  annotation = FALSE,
  annotation_style = "basic",
  annotation_width_ratio = 1 / 50,
  annotation_spacing_ratio = 1 / 3,
  annotation_fontsize = 10,
  annotation_colour = "#48CFCB",
  annotation_fill = "#48CFCB",
  track = FALSE,
  track_width_ratio = 1 / 20,
  track_spacing_ratio = 1 / 2,
  track_fill = "black",
  track_fontsize = 5,
  tad = FALSE,
  tad_is_0_based = FALSE,
  tad_colour = "grey",
  loop = FALSE,
  loop_style = "circle",
  loop_is_0_based = FALSE,
  loop_colour = "black",
  loop_fill = NA,
  concatemer = FALSE,
  concatemer_width_ratio = 1 / 100,
  concatemer_spacing_ratio = 1 / 5,
  expand_xaxis = FALSE,
  expand_left = NULL,
  expand_right = NULL,
  ...
) {
  dat <- x |>
    .checkDataType(scale_column = scale_column, scale_method = scale_method)

  p <- dat |>
    ggplot2::ggplot(
      ggplot2::aes(
        seqnames1 = seqnames1, start1 = start1, end1 = end1,
        seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
      )
    )

  p <- p + geom_hic(...)

  if (ideogram) {
    p <- p + geom_ideogram(
      width_ratio = ideogram_width_ratio, fontsize = ideogram_fontsize,
      colour = ideogram_colour, fill = ideogram_fill, ...
    )
  }

  if (annotation) {
    p <- p + geom_annotation(
      width_ratio = annotation_width_ratio,
      spacing_ratio = annotation_spacing_ratio, style = annotation_style,
      fontsize = annotation_fontsize, colour = annotation_colour,
      fill = annotation_fill, ...
    )
  }

  if (track) {
    p <- p + geom_track(
      width_ratio = track_width_ratio, spacing_ratio = track_spacing_ratio,
      fill = track_fill, fontsize = track_fontsize, ...
    )
  }

  if (tad) {
    p <- p + geom_tad(is_0_based = tad_is_0_based, colour = tad_colour, ...)
  }

  if (loop) {
    p <- p + geom_loop(
      is_0_based = loop_is_0_based, colour = loop_colour, fill = loop_fill,
      style = loop_style, ...
    )
  }

  if (concatemer) {
    p <- p + geom_concatemer(
      width_ratio = concatemer_width_ratio,
      spacing_ratio = concatemer_spacing_ratio, ...
    )
  }

  breaks_labels <- .getBreaksLabels(dat)
  breaks <- breaks_labels$breaks
  labels <- breaks_labels$labels

  range_x <- NULL
  if (expand_xaxis) {
    range_x <- .calculateXRange(dat)

    .resolution <- dat |>
      dplyr::slice(1) |>
      dplyr::mutate(res = end1 - start1) |>
      dplyr::pull(res)

    if (is.null(expand_left)) expand_left <- .resolution * 10
    if (is.null(expand_right)) expand_right <- .resolution * 10
  }

  p <- p +
    theme_hic(
      breaks = breaks, labels = labels, xmin = range_x[1], xmax = range_x[2],
      expand_x = c(expand_left, expand_right)
    )

  p
}

#' @rdname gghic
#' @param scale_column Character string. Name of the column to use for scaling
#'   (e.g., `"balanced"`, `"count"`). Used when input is `GInteractions` or
#'   `data.frame`.
#' @param tad_is_0_based Logical. Whether TAD coordinates are 0-based. Default
#'   is `TRUE`.
#' @param loop_is_0_based Logical. Whether loop coordinates are 0-based. Default
#'   is `TRUE`.
#' @export
methods::setMethod("gghic", "data.frame", .function)

#' @rdname gghic
#' @export
methods::setMethod("gghic", "GInteractions", .function)
