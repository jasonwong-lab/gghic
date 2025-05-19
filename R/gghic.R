scale_data <- function(gis, score_column, scale_method) {
  x <- tibble::as_tibble(gis) |>
    dplyr::mutate(score = scale_method(.data[[score_column]])) |>
    # dplyr::filter(
    #   InteractionSet::pairdist(gis) != 0,
    #   !is.na(InteractionSet::pairdist(gis) != 0)
    # ) |>
    dplyr::filter(!is.na(score), !is.infinite(score)) |>
    dplyr::mutate(score = scales::oob_squish(score, c(min(score), max(score))))

  x
}

check_data_type <- function(data, ...) {
  name_pkg <- get_pkg_name()
  env <- get(".env", envir = asNamespace(name_pkg))

  if (methods::is(data, "GInteractions")) {
    x <- scale_data(data, ...)
    env$gis <- data
  } else if (methods::is(data, "HiCExperiment")) {
    gis <- InteractionSet::interactions(data)
    x <- scale_data(gis, ...)
    env$gis <- InteractionSet::interactions(data)
  } else if (tibble::is_tibble(data) || methods::is(data, "data.frame")) {
    cols_required <- c(
      "seqnames1", "seqnames2", "start1", "end1", "start2", "end2", "score"
    )
    cols_missing <- setdiff(cols_required, colnames(data))
    if (length(cols_missing) > 0) {
      stop(
        "data must have the following columns: ",
        paste(cols_missing, collapse = ", ")
      )
    }
    x <- data
  } else {
    stop("data must be a HiCExperiment object or a tibble/data.frame")
  }
  x
}

calculate_xrange <- function(data) {
  n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))

  if (n_sn == 1 || (n_sn == 2 && all(data$seqnames1 != data$seqnames2))) {
    data <- data
  } else {
    data <- data |>
      adjust_coordinates(list(c(start1 = "start1", end1 = "end1")), FALSE)
  }

  c(min(data$start1), max(data$end1))
}

#' gghic
#'
#' @description A ggplot2 wrapper of [geom_hic()], [geom_ideogram()],
#'   [geom_annotation()], [geom_track()], [geom_tad()], and [geom_loop()]
#'   for easy visualisation of Hi-C/-like data.
#' @inheritParams geom_hic
#' @inheritParams geom_ideogram
#' @inheritParams geom_annotation
#' @inheritParams geom_track
#' @inheritDotParams geom_tad
#' @inheritDotParams geom_loop
#' @param data Either a HiCExperiment object, or a GInteractions object, or a
#'   tibble/data.frame.
#' @param score_column The column name of which the score is calculated.
#'   Default is `"balanced"`.
#' @param scale_method The function to scale the score. Default is `log10`.
#' @param ideogram Whether to add ideogram or not. Default is `FALSE`.
#' @param ideogram_width_ratio The width ratio of the ideogram.
#'   Default is `1/30`.
#' @param ideogram_fontsize The font size of the ideogram. Default is `10`.
#' @param ideogram_colour The color of the ideogram. Default is `"red"`.
#' @param ideogram_fill The fill color of the highlighted region on the
#'   ideogram. Default is `"#FFE3E680"`.
#' @param annotation Whether to add annotation or not. Default is `FALSE`.
#' @param annotation_width_ratio The width ratio of the annotation.
#'   Default is `1/50`.
#' @param annotation_spacing_ratio The spacing ratio of the annotation.
#'   Default is `1/3`.
#' @param  annotation_fontsize The font size of the annotation. Default is `10`.
#' @param annotation_colour The color of the annotation. Default is `"#48CFCB"`.
#' @param annotation_fill The fill color of the annotation.
#'   Default is `"#48CFCB"`.
#' @param track Whether to add track or not. Default is `FALSE`.
#' @param track_width_ratio The width ratio of the track. Default is `1/20`.
#' @param track_spacing_ratio The spacing ratio of the track. Default is `1/2`.
#' @param track_fill The fill color of the track. Default is `"black"`.
#' @param track_fontsize The font size of the track. Default is `5`.
#' @param tad Whether to add TAD or not. Default is `FALSE`.
#' @param tad_is_0based Whether the TAD coordinates are 0-based or not.
#'   Default is `FALSE`.
#' @param tad_colour The color of the TAD. Default is `"grey"`.
#' @param loop Whether to add loop or not. Default is `FALSE`.
#' @param loop_is_0based Whether the loop coordinates are 0-based or not.
#'   Default is `FALSE`.
#' @param loop_colour The color of the loop. Default is `"black"`.
#' @param loop_fill The fill color of the loop. Default is `NA`.
#' @param expand_xaxis Whether to expand the x-axis or not. Default is `FALSE`.
#' @param expand_left The left expansion of the x-axis. Default is `NULL`.
#' @param expand_right The right expansion of the x-axis. Default is `NULL`.
#' @inheritDotParams geom_hic -mapping
#' @inheritDotParams geom_ideogram -mapping -width_ratio -fontsize -colour -fill
#' @inheritDotParams geom_annotation -mapping -width_ratio -spacing_ratio
#'   -fontsize -colour -fill
#' @inheritDotParams geom_track -mapping -width_ratio -spacing_ratio -fill
#'   -fontsize
#' @inheritDotParams geom_tad -mapping -is_0based -colour
#' @inheritDotParams geom_loop -mapping -is_0based -colour -fill
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' library(gghic)
#' library(dplyr)
#' library(HiCExperiment)
#' library(InteractionSet)
#' library(scales)
#' library(scales)
#'
#' dir_cache_gghic <- user_cache_dir(appname = "gghic")
#' url_file <- "https://raw.githubusercontent.com/mhjiang97/gghic-data/refs/heads/master/cooler/chr4_11-5kb.cool"
#' path_file <- file.path(dir_cache_gghic, "chr4_11-5kb.cool")
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
#' url_file <- "https://raw.githubusercontent.com/mhjiang97/gghic-data/refs/heads/master/gtf/gencode-chr4_11.gtf.gz"
#' path_gtf <- glue("{dir_cache_gghic}/gencode-chr4_11.gtf.gz")
#' download.file(url_file, path_gtf)
#'
#' x |>
#'   filter(
#'     center1 > 10000000 & center1 < 11000000 &
#'       center2 > 10000000 & center2 < 11000000
#'   ) |>
#'   gghic(
#'     draw_boundary = TRUE,
#'
#'     ideogram = TRUE, genome = "hg19", highlight = TRUE,
#'     ideogram_fontsize = 7, ideogram_width_ratio = 0.08,
#'
#'     annotation = TRUE, include_ncrna = FALSE, gtf_path = path_gtf,
#'     style = "arrow", maxgap = 100000, annotation_fontsize = 5,
#'     annotation_width_ratio = 0.05,
#'
#'     expand_xaxis = TRUE
#'   )
#' }
#' @export gghic
#' @aliases gghic
gghic <- function(
  data = NULL,
  score_column = "balanced",
  scale_method = log10,

  ideogram = FALSE,
  ideogram_width_ratio = 1 / 30,
  ideogram_fontsize = 10,
  ideogram_colour = "red",
  ideogram_fill = "#FFE3E680",

  annotation = FALSE,
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
  tad_is_0based = FALSE,
  tad_colour = "grey",

  loop = FALSE,
  loop_is_0based = FALSE,
  loop_colour = "black",
  loop_fill = NA,

  expand_xaxis = FALSE,
  expand_left = NULL,
  expand_right = NULL,

  ...
) {
  dat <- data |>
    check_data_type(score_column = score_column, scale_method = scale_method) |>
    tidyr::drop_na(score)

  p <- dat |>
    ggplot2::ggplot(
      ggplot2::aes(
        seqnames1 = seqnames1, start1 = start1, end1 = end1,
        seqnames2 = seqnames2, start2 = start2, end2 = end2,
        fill = score
      )
    )

  p <- p + geom_hic(...)

  if (ideogram) {
    p <- p +
      geom_ideogram(
        width_ratio = ideogram_width_ratio,
        fontsize = ideogram_fontsize,
        colour = ideogram_colour,
        fill = ideogram_fill,
        ...
      )
  }

  if (annotation) {
    p <- p +
      geom_annotation(
        width_ratio = annotation_width_ratio,
        spacing_ratio = annotation_spacing_ratio,
        fontsize = annotation_fontsize,
        colour = annotation_colour,
        fill = annotation_fill,
        ...
      )
  }

  if (track) {
    p <- p +
      geom_track(
        width_ratio = track_width_ratio,
        spacing_ratio = track_spacing_ratio,
        fill = track_fill,
        fontsize = track_fontsize,
        ...
      )
  }

  if (tad) {
    p <- p +
      geom_tad(
        is_0based = tad_is_0based,
        colour = tad_colour,
        ...
      )
  }

  if (loop) {
    p <- p +
      geom_loop(
        is_0based = loop_is_0based,
        colour = loop_colour,
        fill = loop_fill,
        ...
      )
  }

  range_x <- NULL
  if (expand_xaxis) {
    range_x <- dat |>
      calculate_xrange()

    if (is.null(expand_left)) expand_left <- 150000
    if (is.null(expand_right)) expand_right <- 150000
  }

  breaks <- ggplot2::waiver()
  labels <- scales::unit_format(unit = "M", scale = 1e-6)
  n_sn <- length(unique(c(dat$seqnames1, dat$seqnames2)))
  if (n_sn != 1 && (n_sn != 2 || any(dat$seqnames1 == dat$seqnames2))) {
    tmp <- dat |>
      dplyr::group_by(seqnames1) |>
      dplyr::reframe(
        bin = seq(from = min(end1), to = max(end1), length.out = 5)
      ) |>
      dplyr::rename(seqname = seqnames1)

    chroms_add <- dat |> calculate_add_lengths()
    chroms_sub <- dat |> calculate_subtract_lengths()

    tmp2 <- tmp |>
      adjust_coordinates2(chroms_add, chroms_sub, c(bin = "bin"))

    breaks <- tmp2$bin
    labels <- labels(tmp$bin)

    indices <- match(unique(tmp$seqname), tmp$seqname)[-1]
    labels[indices] <- paste0("\n", labels[indices])
  }

  p <- p +
    theme_hic(
      breaks = breaks, labels = labels,
      xmin = range_x[1], xmax = range_x[2],
      expand_x = c(expand_left, expand_right)
    )

  p
}
