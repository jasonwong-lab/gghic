#' Import chromatin contact data
#'
#' @name import-ChromatinContacts
#' @aliases import,ChromatinContacts,ANY,ANY-method
#'
#' @description
#' Imports interaction data from cooler file into ChromatinContacts object.
#'
#' @param con ChromatinContacts object with cooler file path.
#' @param format Character. Balance column name for normalization (default:
#'   `"weight"`). Can be provided as second argument for convenience.
#' @param text Not used (compatibility).
#' @param balance_column Character. Balance column name (overridden by
#'   `format` if provided).
#' @param ... Additional arguments (not used).
#'
#' @return ChromatinContacts with `interactions` slot populated by GInteractions
#'   containing bin_id1, bin_id2, count, and balanced columns.
#'
#' @details
#' Loads Hi-C data from cooler file. If focus regions are specified, only those
#' interactions are loaded. Balance weights (e.g., "weight", "KR", "VC") are
#' used to normalize raw counts. Warning issued if balance column not found.
#'
#' Shorthand: `import(cc, "KR")` is equivalent to
#' `import(cc, balance_column = "KR")`.
#'
#' @examples
#' \dontrun{
#' # Basic import
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # With specific balance column
#' cc <- import(cc, "KR")
#'
#' # Focus region (more efficient)
#' cc <- ChromatinContacts("sample.cool", focus = "chr1:1-5e6") |> import()
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

#' Import multi-way contact data
#'
#' @name import-MultiWayContacts
#' @aliases import,MultiWayContacts,ANY,ANY-method
#'
#' @description
#' Imports multi-way chromatin contact data from a .pairs file into a
#' MultiWayContacts object. The method loads pairwise contacts from long-read
#' sequencing experiments (e.g., Pore-C).
#'
#' @param con A `MultiWayContacts` object containing the path to the pairs file
#'   and optional focus regions.
#' @param format Logical indicating whether to include inter-chromosomal
#'   contacts. This parameter is used as the second positional argument for
#'   convenience. Default is `FALSE`.
#' @param text Not used. Included for compatibility with the generic.
#' @param inter_chrom Logical. If `TRUE`, includes inter-chromosomal contacts.
#'   If `FALSE` (default), only intra-chromosomal contacts are loaded. If
#'   `format` is provided, it overrides this parameter.
#' @param ... Additional arguments (currently not used).
#'
#' @return A `MultiWayContacts` object with the `pairs` slot populated by a
#'   data frame containing:
#'   * read_name: Unique identifier for each read
#'   * chrom1, chrom2: Chromosomes for each contact
#'   * pos1, pos2: Genomic positions for each contact
#'
#' @details
#' The function loads pairwise contact data from a .pairs file format (typically
#' .pairs.gz). If a focus chromosome is specified in the `MultiWayContacts`
#' object, only contacts involving that chromosome are loaded.
#'
#' When `inter_chrom = FALSE`, reads that contain any inter-chromosomal
#' contacts are completely removed from the dataset to ensure clean
#' intra-chromosomal analysis.
#'
#' The pairs file is read using optimized C code for efficient parsing of
#' large files.
#'
#' For convenience, the inter_chrom parameter can be provided as the second
#' positional argument: `import(mc, TRUE)` is equivalent to
#' `import(mc, inter_chrom = TRUE)`.
#'
#' @examples
#' \dontrun{
#' # Basic import (intra-chromosomal only)
#' mc <- MultiWayContacts("sample.pairs.gz")
#' mc <- import(mc)
#'
#' # Using pipe operator
#' mc <- MultiWayContacts("sample.pairs.gz") |> import()
#'
#' # Include inter-chromosomal contacts
#' mc <- import(mc, inter_chrom = TRUE)
#'
#' # Convenient syntax
#' mc <- import(mc, TRUE)
#'
#' # Import with focus on specific chromosome
#' mc <- MultiWayContacts(
#'   "sample.pairs.gz",
#'   focus = "chr1"
#' ) |> import()
#'
#' # Check imported data
#' head(mc@pairs)
#' }
#'
#' @export
methods::setMethod(
  "import", "MultiWayContacts", function(
    con, format, text, inter_chrom = FALSE, ...
  ) {
    if (!missing(format)) {
      inter_chrom <- format
    }
    con@pairs <- .loadPairsData(NULL, con@pairs_path, con@focus, inter_chrom)

    con
  }
)
