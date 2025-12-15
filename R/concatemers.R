#' Convert multi-way concatemer reads to pairwise Hi-C interactions
#'
#' @description
#' Transforms multi-way chromatin contact reads (concatemers from Pore-C,
#' Tri-C, Multi-C, etc.) into traditional pairwise Hi-C format. Generates all
#' pairwise combinations of fragments within each read, bins to specified
#' resolution, and creates a GInteractions object compatible with standard
#' Hi-C analysis tools.
#'
#' @param grs GRanges object containing concatemer fragment coordinates. Must
#'   include metadata column identifying which fragments belong to the same
#'   read (specified by `read_group` parameter).
#' @param region GRanges object defining genomic region(s) of interest. Only
#'   reads overlapping this region are included (default: NULL for genome-wide
#'   analysis).
#' @param bin_size Integer. Genomic bin size in base pairs for aggregating
#'   contacts. Typical values: 5000-1000000 (default: 100000 = 100kb).
#' @param read_group Character. Name of metadata column in `grs` that contains
#'   read identifiers grouping fragments from the same sequencing read
#'   (default: `"read_name"`).
#'
#' @return GInteractions object containing:
#'   * Binned pairwise interactions from all concatemers
#'   * `count` metadata column: number of reads supporting each interaction
#'   * Compatible with [gghic()], [ChromatinContacts()], and other Hi-C tools
#'
#' @details
#' ## Algorithm
#' 1. Group fragments by read identifier (`read_group`)
#' 2. For each read with ≥2 fragments, generate all pairwise combinations
#' 3. Bin both anchors of each pair to specified resolution
#' 4. Aggregate identical bin pairs, counting supporting reads
#' 5. Remove self-interactions (same bin pairs)
#'
#' ## Multi-way to pairwise conversion
#' A read contacting N fragments generates N×(N-1)/2 pairwise interactions.
#' For example:
#' * 2-way contact → 1 pair
#' * 3-way contact → 3 pairs
#' * 4-way contact → 6 pairs
#'
#' This inflation should be considered when comparing multi-way and traditional
#' Hi-C datasets.
#'
#' ## Memory considerations
#' For large datasets, use `region` parameter to process chromosomes
#' individually and reduce memory usage.
#'
#' @examples
#' \dontrun{
#' # Load concatemer data
#' concatemers <- readRDS("concatemers.rds")
#'
#' # Convert to pairwise at 100kb resolution
#' gis <- concatemers2Gis(concatemers, bin_size = 100000)
#'
#' # Visualize as Hi-C map
#' gghic(gis)
#'
#' # Focus on specific region at higher resolution
#' region <- GRanges("chr1:1000000-5000000")
#' gis_region <- concatemers2Gis(
#'   concatemers,
#'   region = region,
#'   bin_size = 10000
#' )
#'
#' # Process by chromosome to manage memory
#' chr_gis <- lapply(paste0("chr", 1:22), function(chr) {
#'   region <- GRanges(seqnames = chr,
#'                     ranges = IRanges(1, 250000000))
#'   concatemers2Gis(concatemers, region, bin_size = 100000)
#' })
#'
#' # Use with ChromatinContacts
#' gis_all <- do.call(c, chr_gis)
#' # Then visualize or analyze
#' }
#'
#' @seealso [ChromatinContacts()], [gghic()], [MultiWayContacts()]
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
