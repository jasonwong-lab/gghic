#' Convert concatemer reads to pairwise interactions
#'
#' @description
#' Converts multi-way contact reads (concatemers) from Pore-C or other
#' multi-contact sequencing to pairwise interactions in a `GInteractions`
#' object. Useful for comparing multi-way data with traditional Hi-C.
#'
#' @param grs A `GRanges` object containing concatemer reads. Must have a
#'   metadata column (default `read_name`) identifying which fragments belong
#'   to the same read.
#' @param region A `GRanges` object specifying the genomic region of interest.
#'   Default is `NULL`, which uses the entire range of `grs`.
#' @param bin_size Integer. The size of bins (in base pairs) to divide the
#'   genomic region into. Default is `100000` (100 kb).
#' @param read_group Character string. The name of the metadata column in `grs`
#'   that defines read groups (typically read names or IDs). Default is
#'   `"read_name"`.
#'
#' @return A `GInteractions` object containing all pairwise combinations of
#'   fragments within each concatemer read, binned at the specified resolution.
#'
#' @details
#' For each multi-way contact (concatemer), this function generates all possible
#' pairwise combinations of the fragments and bins them to the specified
#' resolution. This allows visualization of Pore-C or other multi-contact data
#' alongside traditional Hi-C data.
#'
#' The resulting interaction counts represent the number of concatemer reads
#' supporting each pairwise interaction, which differs from traditional Hi-C
#' counting.
#'
#' **Workflow:**
#' 1. Groups fragments by read identifier
#' 2. For each read with multiple fragments, generates all pairwise combinations
#' 3. Bins fragment coordinates to specified resolution
#' 4. Counts interactions per bin pair
#'
#' @examples
#' \dontrun{
#' # Load concatemer data
#' concatemers <- readRDS("concatemers.rds")
#'
#' # Convert to pairwise interactions at 100kb resolution
#' gis <- concatemers2Gis(concatemers, bin_size = 100000)
#'
#' # Visualize
#' gghic(gis)
#'
#' # Focus on specific region
#' region <- GRanges("chr1:1000000-5000000")
#' gis_sub <- concatemers2Gis(concatemers, region = region, bin_size = 10000)
#' }
#'
#' @export
#' @aliases concatemers2Gis
concatemers2Gis <- function(
  grs, region = NULL, bin_size = 100000, read_group = "read_name"
) {
  if (!is.null(region) && !inherits(region, "GRanges")) {
    stop("'region' must be a GRanges object or NULL.")
  }

  if (!is.null(region)) {
    hits <- GenomicRanges::findOverlaps(grs, region)
    grs <- grs[S4Vectors::queryHits(hits)]
    df_region <- tibble::as_tibble(region)
  }

  if (length(grs) == 0) {
    stop("No concatemer reads overlap the specified region.")
  }

  df <- tibble::as_tibble(grs)

  if (!is.null(region)) {
    df_query <- df_region
  } else {
    df_query <- df
  }

  chrs <- df_query$seqnames |>
    unique() |>
    as.character()
  names(chrs) <- chrs
  starts <- purrr::map_int(chrs, ~ min(df_query$start[df_query$seqnames == .x]))
  ends <- purrr::map_int(chrs, ~ max(df_query$end[df_query$seqnames == .x]))

  starts_bin <- purrr::map2(
    starts, ends, ~ seq(from = .x, to = .y - bin_size + 1, by = bin_size)
  )

  grs_bin <- purrr::map2_dfr(chrs, starts_bin, function(x, y) {
    tibble::tibble(
      seqnames = x,
      start = y,
      end = y + bin_size - 1
    )
  }) |>
    methods::as("GRanges")

  df_bin <- tibble::as_tibble(grs_bin) |>
    dplyr::mutate(bin_id = dplyr::row_number())

  pairs_all <- expand.grid(
    bin1_id = df_bin$bin_id,
    bin2_id = df_bin$bin_id
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(bin1_id <= bin2_id) |>
    dplyr::arrange(bin1_id, bin2_id)

  hits <- GenomicRanges::findOverlaps(grs, grs_bin)
  grs_with_bins <- grs[S4Vectors::queryHits(hits)]
  grs_with_bins$bin_id <- S4Vectors::subjectHits(hits)

  df_with_bins <- tibble::as_tibble(grs_with_bins)

  concatemer_interactions <- split(df_with_bins, df_with_bins[[read_group]]) |>
    purrr::map_dfr(function(x) {
      if (nrow(x) < 2) {
        return(NULL)
      }
      bin_ids <- unique(x$bin_id)
      if (length(bin_ids) < 2) {
        return(NULL)
      }
      .pairs <- combn(bin_ids, 2, simplify = FALSE)
      purrr::map_dfr(.pairs, function(pair) {
        tibble::tibble(bin1_id = min(pair), bin2_id = max(pair))
      })
    }) |>
    dplyr::group_by(bin1_id, bin2_id) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  if (nrow(concatemer_interactions) == 0) {
    counts_interaction <- pairs_all |>
      dplyr::mutate(count = 0)
  } else {
    counts_interaction <- pairs_all |>
      dplyr::left_join(concatemer_interactions, by = c("bin1_id", "bin2_id")) |>
      dplyr::mutate(count = tidyr::replace_na(count, 0))
  }

  anchors1 <- GenomicRanges::GRanges(
    seqnames = df_bin$seqnames[counts_interaction$bin1_id],
    ranges = IRanges::IRanges(
      df_bin$start[counts_interaction$bin1_id],
      df_bin$end[counts_interaction$bin1_id]
    )
  )

  anchors2 <- GenomicRanges::GRanges(
    seqnames = df_bin$seqnames[counts_interaction$bin2_id],
    ranges = IRanges::IRanges(
      df_bin$start[counts_interaction$bin2_id],
      df_bin$end[counts_interaction$bin2_id]
    )
  )

  gis <- InteractionSet::GInteractions(anchors1, anchors2)
  S4Vectors::mcols(gis) <- counts_interaction[
    , c("count", "bin1_id", "bin2_id")
  ]

  gis
}
