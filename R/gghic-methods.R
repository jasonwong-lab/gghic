#' Scale and transform Hi-C interaction data for visualization
#'
#' @description
#' Transforms and scales chromatin interaction data to prepare it for
#' visualization. Applies user-defined scaling functions (e.g., log
#' transformation) to interaction scores and handles missing values.
#'
#' @param data Input data in one of these formats:
#'   * ChromatinContacts object with imported interactions
#'   * GInteractions object with score metadata
#'   * data.frame or tibble with columns: seqnames1, start1, end1, seqnames2,
#'     start2, end2, plus score column
#' @param scale_column Character. Name of column containing values to scale.
#'   Common options:
#'   * `"balanced"`: ICE-normalized counts (recommended)
#'   * `"count"`: raw contact counts
#'   * Any other numeric metadata column
#' @param scale_method Function to apply for transformation. Common options:
#'   * `log10`: log10 transformation (default for most Hi-C data)
#'   * `log2`: log2 transformation
#'   * `log1p`: log(x + 1) transformation (handles zeros)
#'   * `identity` or `function(x) x`: no transformation
#'   * Custom function: any function that takes numeric vector and returns
#'     numeric vector
#' @param remove_na Logical. Whether to remove rows with NA or infinite values
#'   after scaling (default: FALSE). Set TRUE to remove missing data that could
#'   cause visualization issues.
#'
#' @return Tibble (data frame) with standardized columns:
#'   * `seqnames1`, `start1`, `end1`: First anchor coordinates
#'   * `seqnames2`, `start2`, `end2`: Second anchor coordinates
#'   * `score`: Transformed and scaled values
#'
#' @details
#' ## Processing steps
#' 1. Convert input to tibble format
#' 2. Apply `scale_method` function to `scale_column`
#' 3. Create new `score` column with transformed values
#' 4. Optionally remove NA/infinite values
#' 5. Squish extreme outliers to prevent visualization artifacts
#'
#' ## Recommended scaling
#' For typical Hi-C data visualization:
#' * Use `"balanced"` column with `log10` transformation
#' * Set `remove_na = TRUE` to handle bins with no coverage
#'
#' ## Custom transformations
#' You can provide any function for scaling:
#' ```r
#' # Square root transformation
#' scaleData(cc, "count", sqrt)
#'
#' # Custom transformation
#' scaleData(cc, "balanced", function(x) log2(x + 0.1))
#' ```
#'
#' @examples
#' \dontrun{
#' # Load Hi-C data
#' cc <- ChromatinContacts("file.cool") |> import()
#'
#' # Standard log10 scaling of balanced data
#' scaled_data <- scaleData(cc, "balanced", log10)
#'
#' # Raw counts without transformation
#' scaled_raw <- scaleData(cc, "count", function(x) x)
#'
#' # Log2 scaling with NA removal
#' scaled_clean <- scaleData(cc, "balanced", log2, remove_na = TRUE)
#'
#' # Use with plotting
#' library(ggplot2)
#' ggplot() +
#'   geom_hic(data = scaleData(cc, "balanced", log10),
#'            aes(seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'                seqnames2 = seqnames2, start2 = start2, end2 = end2,
#'                fill = score))
#' }
#'
#' @seealso [gghic()], [geom_hic()], [ChromatinContacts()]
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
#'
#' @name gghic
#' @aliases gghic,ChromatinContacts-method
#'
#' @description
#' High-level wrapper for publication-ready Hi-C visualizations with genomic
#' features.
#'
#' @param x ChromatinContacts, GInteractions, or data.frame with interactions.
#' @param scale_method Function for data transformation (default: log10).
#' @param scale_column Character. Column to scale for GInteractions/data.frame
#'   (default: `"balanced"`).
#' @param ideogram Logical. Add chromosome ideogram (default: FALSE).
#' @param ideogram_width_ratio Numeric. Ideogram height ratio (default: 1/30).
#' @param ideogram_fontsize Numeric. Ideogram font size (default: 10).
#' @param ideogram_colour Character. Highlight color (default: `"red"`).
#' @param ideogram_fill Character. Highlight fill (default: `"#FFE3E680"`).
#' @param annotation Logical. Add gene annotations (default: FALSE). Requires
#'   `gtf_path` in `...`.
#' @param annotation_style Character. `"basic"` or `"arrow"` (default:
#'   `"basic"`).
#' @param annotation_width_ratio Numeric. Track height (default: 1/50).
#' @param annotation_spacing_ratio Numeric. Gene spacing (default: 1/3).
#' @param annotation_fontsize Numeric. Label size (default: 10).
#' @param annotation_colour Character. Feature color (default: `"#48CFCB"`).
#' @param annotation_fill Character. Feature fill (default: `"#48CFCB"`).
#' @param track Logical. Add genomic tracks (default: FALSE). Requires `tracks`
#'   in object or `data_paths` in `...`.
#' @param track_width_ratio Numeric. Track height (default: 1/20).
#' @param track_spacing_ratio Numeric. Track spacing (default: 1/2).
#' @param track_fill Character. Track colors (default: `"black"`).
#' @param track_fontsize Numeric. Label size (default: 5).
#' @param tad Logical. Add TAD boundaries (default: FALSE). Requires `TADs` in
#'   object or `tad_path` in `...`.
#' @param tad_colour Character. TAD color (default: `"grey"`).
#' @param loop Logical. Add loops (default: FALSE). Requires `loops` in object
#'   or `loop_path` in `...`.
#' @param loop_style Character. `"circle"` or `"arc"` (default: `"circle"`).
#' @param loop_colour Character. Loop color (default: `"black"`).
#' @param loop_fill Fill color (default: NA).
#' @param concatemer Logical. Add multi-way contacts (default: FALSE).
#' @param concatemer_width_ratio Numeric. Track height (default: 1/100).
#' @param concatemer_spacing_ratio Numeric. Spacing (default: 1/5).
#' @param expand_xaxis Logical. Expand x-axis (default: FALSE).
#' @param expand_left Numeric. Left expansion in bp (default: 10×resolution).
#' @param expand_right Numeric. Right expansion in bp (default: 10×resolution).
#' @param ... Additional arguments for geom functions.
#' @return ggplot2 object.
#'
#' @details
#' Automatically handles: data scaling, coordinate transformation, feature
#' integration, axis labels, default theme.
#'
#' Use individual `geom_*()` functions for more control: `geom_hic()`,
#' `geom_ideogram()`, `geom_annotation()`, `geom_track()`, `geom_tad()`,
#' `geom_loop()`, `geom_concatemer()`.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' cc <- ChromatinContacts("file.cool") |> import()
#' gghic(cc)
#' gghic(cc["chr4:0-50000000"])
#'
#' # With features
#' gghic(cc, ideogram = TRUE, genome = "hg19")
#' gghic(cc, tad = TRUE, tad_path = "tads.bed")
#' gghic(cc, loop = TRUE, loop_path = "loops.bedpe")
#'
#' # With tracks
#' features(cc, "tracks") <- GRangesList(
#'   H3K27ac = rtracklayer::import("track.bw")
#' )
#' gghic(cc, track = TRUE, track_fill = "blue")
#' }
#'
#' @seealso
#' [ChromatinContacts], [geom_hic()], [theme_hic()], [scaleData()]
#'
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
#'
#' @param scale_column Character string. Name of the column to use for scaling
#'   (e.g., `"balanced"`, `"count"`). Used when input is `GInteractions` or
#'   `data.frame`.
#' @param tad_is_0_based Logical. Whether TAD coordinates are 0-based. Default
#'   is `TRUE`.
#' @param loop_is_0_based Logical. Whether loop coordinates are 0-based. Default
#'   is `TRUE`.
#'
#' @export
methods::setMethod("gghic", "data.frame", .function)

#' @rdname gghic
#'
#' @export
methods::setMethod("gghic", "GInteractions", .function)
