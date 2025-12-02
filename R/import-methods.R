#' Import chromatin contact data
#'
#' @name import-ChromatinContacts
#' @aliases import,ChromatinContacts,ANY,ANY-method
#'
#' @description
#' Imports chromatin interaction data from a cooler file into a
#' ChromatinContacts object. The method loads bin-level interactions and
#' applies balancing weights if available.
#'
#' @param con A `ChromatinContacts` object containing the path to the cooler
#'   file and optional focus regions.
#' @param format Character string specifying the balance column name to use
#'   for normalization. This parameter is used as the second positional
#'   argument for convenience. Default is `"weight"`.
#' @param text Not used. Included for compatibility with the generic.
#' @param balance_column Character string specifying the balance column name
#'   to use for normalization. If `format` is provided, it overrides this
#'   parameter. Default is `"weight"`.
#' @param ... Additional arguments (currently not used).
#'
#' @return A `ChromatinContacts` object with the `interactions` slot populated
#'   by a `GInteractions` object containing:
#'   * bin_id1: Bin ID for the first anchor
#'   * bin_id2: Bin ID for the second anchor
#'   * count: Raw interaction count
#'   * balanced: Normalized interaction count (if balance weights are available)
#'
#' @details
#' The function loads chromatin interaction data from a cooler file. If a
#' focus region is specified in the `ChromatinContacts` object, only
#' interactions within or between those regions are loaded. Otherwise, all
#' genome-wide interactions are imported.
#'
#' The balance column (e.g., "weight", "KR", "VC") is used to normalize
#' the raw counts. If the specified balance column is not found in the
#' cooler file, a warning is issued and unbalanced counts are returned.
#'
#' For convenience, the balance column name can be provided as the second
#' positional argument: `import(cc, "KR")` is equivalent to
#' `import(cc, balance_column = "KR")`.
#'
#' @examples
#' \dontrun{
#' # Basic import from cooler file
#' cc <- ChromatinContacts("sample.cool")
#' cc <- import(cc)
#'
#' # Using pipe operator
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # With specific balance column
#' cc <- import(cc, balance_column = "KR")
#'
#' # Convenient syntax for balance column
#' cc <- import(cc, "VC")
#'
#' # Import with focus region (more efficient)
#' cc <- ChromatinContacts(
#'   "sample.cool",
#'   focus = "chr1:1000000-5000000"
#' ) |> import()
#'
#' # Multi-resolution cooler
#' cc <- ChromatinContacts(
#'   "sample.mcool",
#'   resolution = 10000
#' ) |> import()
#'
#' # Check imported data
#' interactions(cc)
#' length(interactions(cc))
#' }
#'
#' @export
methods::setMethod(
  "import", "ChromatinContacts", function(
    con, format, text, balance_column = "weight", ...
  ) {
    if (!missing(format)) {
      balance_column <- format
    }

    path <- con@cooler_path
    resolution <- con@resolution
    focus <- con@focus

    bins <- .getBins(path, resolution)

    if (is.null(focus)) {
      df <- tibble::tibble(
        bin1_id = .queryCool(path, "pixels/bin1_id", resolution),
        bin2_id = .queryCool(path, "pixels/bin2_id", resolution),
        count = .queryCool(path, "pixels/count", resolution)
      )
    } else {
      max_offset <- .queryCool(path, "indexes/bin1_offset", resolution) |>
        max(na.rm = TRUE)

      df <- purrr::map_dfr(seq_along(focus), function(i) {
        gi <- focus[i]
        anchor_1 <- InteractionSet::anchors(gi, "first")
        anchor_2 <- InteractionSet::anchors(gi, "second")
        is_equal <- IRanges::overlapsAny(anchor_1, anchor_2, type = "equal")

        if (is_equal) {
          result <- .getChunks(path, resolution, bins, anchor_1, max_offset, TRUE)

          .getPixels(path, resolution, result$chunks) |>
            dplyr::filter(bin2_id %in% (result$idx - 1))
        } else {
          result_1 <- .getChunks(path, resolution, bins, anchor_1, max_offset)
          result_2 <- .getChunks(path, resolution, bins, anchor_2, max_offset)

          chunks_valid <- .queryCool(
            path, "pixels/bin1_id", resolution, result_2$chunks
          ) |> unique()

          .getPixels(path, resolution, result_1$chunks) |>
            dplyr::filter(bin2_id %in% chunks_valid)
        }
      })
    }

    gis <- InteractionSet::GInteractions(
      bins[df$bin1_id + 1], bins[df$bin2_id + 1],
      count = df$count
    )
    gis$bin_id1 <- df$bin1_id
    gis$bin_id2 <- df$bin2_id

    col_balance_1 <- paste0("anchor1.", balance_column)
    col_balance_2 <- paste0("anchor2.", balance_column)
    balance_1 <- mcols(gis)[[col_balance_1]]
    balance_2 <- mcols(gis)[[col_balance_2]]
    if (!is.null(balance_1) && !is.null(balance_2)) {
      gis$balanced <- gis$count * balance_1 * balance_2
    } else {
      warning(
        "Balance column '", balance_column,
        "' not found in the cooler file. Returning unbalanced counts."
      )
      gis$balanced <- gis$count
    }

    S4Vectors::mcols(gis) <- S4Vectors::mcols(gis)[
      , c("bin_id1", "bin_id2", "count", "balanced")
    ]

    con@interactions <- gis

    if (is.null(con@resolution)) {
      message("Setting resolution based on the bin size.")
      con@resolution <- width(bins)[1]
    }

    con
  }
)
