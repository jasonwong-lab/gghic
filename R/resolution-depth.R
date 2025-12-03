#' Calculate Resolution and Depth for Pairs Data
#'
#' Functions to analyze genomic interaction data resolution and coverage depth.
#'
#' @param pairs A data frame or tibble with columns:
#'   \code{read_name}, \code{chrom1}, \code{pos1}, \code{chrom2}, \code{pos2}.
#'   Can also have additional columns like \code{strand1}, \code{strand2}, etc.
#'   For large files, use `pairs_file` parameter instead to read in chunks.
#' @param pairs_file Character path to pairs file for chunked reading.
#'   File can be gzip compressed (.gz).
#' @param cache External pointer to cached pairs data (from `readPairsCache()`).
#'   If provided, data is read from cache instead of file for much faster
#'   repeated operations.
#' @param return_cache Logical, if TRUE, `findOptResChunked()`
#'   returns a list with both the optimal bin size and the cache pointer.
#'   Default: FALSE.
#' @param bin_size Integer bin size in base pairs.
#' @param min_contacts Integer minimum number of contacts required for a bin.
#'   Default: 1000.
#' @param min_bin Integer minimum bin size for search. Default: 1000.
#' @param max_bin Integer maximum bin size for search. Default: 5000000.
#' @param target_coverage Numeric target coverage fraction (0-1). Default: 0.8.
#'
#' @return
#' - `calculateResolutionDepth()`: A tibble with columns \code{chrom},
#'   \code{bin}, and \code{count}.
#' - `calculateGenomeCoverage()`: A numeric value representing the fraction
#'   of bins with >= min_contacts.
#' - `findOptimalResolution()`: An integer representing the optimal bin size.
#' - `calcResDepthChunked()`: Same as above, reads file in chunks.
#' - `calcGenomeCovChunked()`: Coverage calculation from file.
#' - `findOptResChunked()`: Optimal resolution for large files.
#'
#' @details
#' These functions analyze Hi-C pairs data to determine optimal genomic resolution.
#' `calculateResolutionDepth()` bins genomic positions and counts unique
#' interactions per bin. `calculateGenomeCoverage()` computes the fraction of
#' bins meeting a contact threshold. `findOptimalResolution()` uses binary
#' search to find the smallest bin size achieving target coverage.
#'
#' For large files that don't fit in memory, use the `_chunked` variants with
#' `pairs_file` parameter. These read the file in chunks and aggregate results.
#'
#' For large datasets, C-accelerated computation is used when available.
#'
#' @examples
#' \dontrun{
#' # Load pairs data
#' pairs <- read.table("contact_matrix.txt",
#'   header = FALSE,
#'   col.names = c(
#'     "read_name", "strand1", "chrom1", "pos1", "frag1",
#'     "strand2", "chrom2", "pos2", "frag2", "mapq1", "mapq2"
#'   )
#' )
#'
#' # Calculate resolution and depth
#' depth_10kb <- calculateResolutionDepth(pairs, bin_size = 10000)
#'
#' # For large file, read in chunks
#' depth_10kb <- calcResDepthChunked(
#'   pairs_file = "contact_matrix.txt",
#'   bin_size = 10000
#' )
#'
#' # Calculate coverage
#' coverage <- calculateGenomeCoverage(pairs, bin_size = 10000)
#'
#' # Find optimal resolution
#' opt_bin <- findOptimalResolution(pairs, target_coverage = 0.8)
#'
#' # For large file - method 1: let function read file each time
#' opt_bin <- findOptResChunked(
#'   pairs_file = "contact_matrix.txt",
#'   target_coverage = 0.8
#' )
#'
#' # For large file - method 2: cache once, reuse many times (FASTER!)
#' cache <- readPairsCache("contact_matrix.txt")
#'
#' # Find optimal resolution using cache
#' opt_bin <- findOptResChunked(cache = cache, target_coverage = 0.8)
#'
#' # Reuse cache for multiple analyses (no file I/O!)
#' depth_10kb <- calcResDepthChunked(cache = cache, bin_size = 10000)
#' depth_50kb <- calcResDepthChunked(cache = cache, bin_size = 50000)
#' coverage <- calcGenomeCovChunked(cache = cache, bin_size = 10000)
#'
#' # Method 3: Get cache back from findOptResChunked
#' result <- findOptResChunked(
#'   pairs_file = "contact_matrix.txt",
#'   target_coverage = 0.8,
#'   return_cache = TRUE
#' )
#' opt_bin <- result$bin_size
#' cache <- result$cache
#'
#' # Now reuse the cache
#' depth <- calcResDepthChunked(cache = cache, bin_size = opt_bin)
#' }
#'
#' @name resolution-depth
#' @rdname resolution-depth
NULL

#' @export
#' @rdname resolution-depth
calculateResolutionDepth <- function(pairs, bin_size) {
  required_cols <- c("read_name", "chrom1", "pos1", "chrom2", "pos2")
  if (!all(required_cols %in% names(pairs))) {
    stop("pairs must contain columns: ", paste(required_cols, collapse = ", "))
  }

  if (.Call("has_resolution_depth_c", PACKAGE = "gghic")) {
    result <- .Call(
      "calculate_resolution_depth_c",
      pairs$read_name,
      pairs$chrom1,
      pairs$pos1,
      pairs$chrom2,
      pairs$pos2,
      as.integer(bin_size),
      PACKAGE = "gghic"
    )
    return(tibble::as_tibble(result))
  }

  df_bins <- pairs |>
    dplyr::mutate(
      bin1 = ceiling(pos1 / bin_size),
      bin2 = ceiling(pos2 / bin_size)
    ) |>
    dplyr::select("read_name", "chrom1", "bin1", "chrom2", "bin2") |>
    dplyr::distinct()

  df1 <- df_bins |>
    dplyr::select(chrom = "chrom1", bin = "bin1")

  df2 <- df_bins |>
    dplyr::select(chrom = "chrom2", bin = "bin2")

  result <- dplyr::bind_rows(df1, df2) |>
    dplyr::group_by(chrom, bin) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    tibble::as_tibble()

  result
}

#' @export
#' @rdname resolution-depth
calculateGenomeCoverage <- function(pairs, bin_size, min_contacts = 1000) {
  counts <- calculateResolutionDepth(pairs, bin_size)

  total_bins <- nrow(counts)
  if (total_bins == 0) {
    return(0.0)
  }

  bins_above_threshold <- sum(counts$count >= min_contacts)

  bins_above_threshold / total_bins
}

#' @export
#' @rdname resolution-depth
findOptimalResolution <- function(
  pairs, min_bin = 1000, max_bin = 5000000, target_coverage = 0.8,
  min_contacts = 1000
) {
  left <- min_bin
  right <- max_bin
  best_bin <- right

  iteration <- 0
  while (left <= right) {
    iteration <- iteration + 1
    mid <- (left + right) %/% 2
    coverage <- calculateGenomeCoverage(pairs, mid, min_contacts)

    message(sprintf(
      "Iteration %d: Testing bin size %s bp: %.2f%% coverage",
      iteration, format(mid, big.mark = ","), coverage * 100
    ))

    if (coverage >= target_coverage) {
      best_bin <- mid
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }

  best_bin
}

#' Plot Resolution vs Coverage
#'
#' Generate a boxplot showing coverage of bins for different bin sizes.
#'
#' @param pairs A data frame or tibble with columns:
#'   \code{read_name}, \code{chrom1}, \code{pos1}, \code{chrom2}, \code{pos2}.
#' @param bin_sizes Integer vector of bin sizes to test (in base pairs).
#'   Default: c(1e3, 2e3, 5e3, 1e4, 2e4, 2.5e4, 5e4, 1e5, 1.5e5, 2e5, 2.5e5, 5e5, 1e6)
#' @param min_contacts Integer minimum number of contacts. Default: 1000.
#' @param title Character title for the plot. Default: "Coverage by Resolution".
#'
#' @return A ggplot2 object showing coverage distribution across bin sizes.
#'
#' @details
#' Creates a boxplot visualization where each panel represents a bin size,
#' and the distribution shows the number of contacts per bin. This helps
#' visualize how resolution affects coverage across the genome.
#'
#' @examples
#' \dontrun{
#' # Plot coverage for default bin sizes
#' p <- plotResolutionCoverage(pairs)
#' print(p)
#'
#' # Custom bin sizes
#' p <- plotResolutionCoverage(pairs,
#'   bin_sizes = c(10000, 50000, 100000, 500000)
#' )
#' print(p)
#' }
#'
#' @export
plotResolutionCoverage <- function(
  pairs, bin_sizes = c(
    1e3, 2e3, 5e3, 1e4, 2e4, 2.5e4, 5e4, 1e5, 1.5e5, 2e5, 2.5e5, 5e5, 1e6
  ), min_contacts = 1000, title = "Coverage by Resolution"
) {
  # Calculate depths for all bin sizes
  results <- purrr::map_df(
    bin_sizes,
    function(bs) {
      calculateResolutionDepth(pairs, bs) |>
        dplyr::mutate(bin_size = bs)
    }
  )

  results <- results |>
    dplyr::mutate(
      bin_size_label = dplyr::case_when(
        bin_size >= 1e6 ~ sprintf("%.1fMb", bin_size / 1e6),
        bin_size >= 1e3 ~ sprintf("%.0fKb", bin_size / 1e3),
        TRUE ~ sprintf("%.0fbp", bin_size)
      ),
      bin_size_label = factor(bin_size_label,
        levels = unique(bin_size_label)
      )
    )

  # Create boxplot
  p <- ggplot2::ggplot(results, ggplot2::aes(x = bin_size_label, y = count)) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7) +
    ggplot2::geom_hline(
      yintercept = min_contacts, linetype = "dashed",
      color = "red", alpha = 0.5
    ) +
    ggplot2::labs(
      title = title,
      x = "Bin Size",
      y = "Number of Contacts per Bin"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    )

  p
}

#' Plot Coverage Fraction by Resolution
#'
#' Create a line plot showing how the fraction of bins meeting contact
#' threshold changes with resolution.
#'
#' @param pairs A data frame or tibble with pairs data.
#' @param bin_sizes Integer vector of bin sizes to test. See [plotResolutionCoverage].
#' @param min_contacts Integer minimum contacts threshold. Default: 1000.
#' @param title Character plot title.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' p <- plotCoverageCurve(pairs)
#' print(p)
#' }
#'
#' @export
plotCoverageCurve <- function(
  pairs, bin_sizes = c(
    1e3, 2e3, 5e3, 1e4, 2e4, 2.5e4, 5e4, 1e5, 1.5e5, 2e5, 2.5e5, 5e5, 1e6
  ), min_contacts = 1000, title = "Genome Coverage by Resolution"
) {
  coverage_data <- purrr::map_df(
    bin_sizes,
    function(bs) {
      tibble::tibble(
        bin_size = bs,
        coverage = calculateGenomeCoverage(pairs, bs, min_contacts)
      )
    }
  )

  coverage_data <- coverage_data |>
    dplyr::mutate(
      bin_size_kb = bin_size / 1e3,
      bin_size_label = dplyr::case_when(
        bin_size >= 1e6 ~ sprintf("%.1fMb", bin_size / 1e6),
        bin_size >= 1e3 ~ sprintf("%.0fKb", bin_size / 1e3),
        TRUE ~ sprintf("%.0fbp", bin_size)
      )
    )

  p <- ggplot2::ggplot(
    coverage_data, ggplot2::aes(x = bin_size_kb, y = coverage)
  ) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue", size = 3) +
    ggplot2::labs(
      title = title,
      x = "Bin Size (Kb)",
      y = "Fraction of Bins with >= min_contacts"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray90"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(labels = scales::percent)

  p
}

#' Read Pairs File into Cache
#'
#' Load a pairs file into C memory for fast repeated analysis.
#' The cache persists until garbage collected or explicitly cleared.
#'
#' @param pairs_file Character path to pairs file.
#'
#' @return An external pointer to the cached data.
#'
#' @examples
#' \dontrun{
#' # Create cache once
#' cache <- readPairsCache("data.pairs.gz")
#'
#' # Reuse cache for multiple operations
#' depth1 <- calcResDepthChunked(cache = cache, bin_size = 10000)
#' depth2 <- calcResDepthChunked(cache = cache, bin_size = 50000)
#' coverage <- calcGenomeCovChunked(cache = cache, bin_size = 10000)
#'
#' # Cache is automatically freed when R session ends
#' # Or explicitly remove it:
#' rm(cache)
#' gc()
#' }
#'
#' @export
#' @rdname resolution-depth
readPairsCache <- function(pairs_file) {
  pairs_file <- normalizePath(pairs_file, mustWork = TRUE)
  message("Reading positions from file into C memory...")

  cache_ptr <- tryCatch(
    {
      .Call("read_pairs_cache_c", pairs_file, PACKAGE = "gghic")
    },
    error = function(e) {
      stop("Failed to read positions from file: ", e$message)
    }
  )

  message("Cache created successfully")
  cache_ptr
}

#' @export
#' @rdname resolution-depth
calcResDepthChunked <- function(
  pairs_file = NULL, bin_size, cache = NULL
) {
  # Use cache if provided, otherwise read file
  if (!is.null(cache)) {
    if (!inherits(cache, "externalptr")) {
      stop("cache must be an external pointer created by readPairsCache()")
    }

    # Call C function that returns full depth table from cache
    tryCatch(
      {
        result <- .Call(
          "process_pairs_from_cache_c",
          cache,
          as.integer(bin_size),
          PACKAGE = "gghic"
        )
        return(tibble::as_tibble(result))
      },
      error = function(e) {
        stop("Failed to calculate depth from cache: ", e$message)
      }
    )
  }

  # Original file-based implementation
  if (is.null(pairs_file)) {
    stop("Either pairs_file or cache must be provided")
  }

  pairs_file <- normalizePath(pairs_file, mustWork = TRUE)

  tryCatch(
    {
      message("Using optimized C implementation for file reading...")
      result <- .Call(
        "process_pairs_file_c",
        pairs_file,
        as.integer(bin_size),
        PACKAGE = "gghic"
      )
      return(tibble::as_tibble(result))
    },
    error = function(e) {
      message(
        paste0(
          "C implementation unavailable, ",
          "falling back to R-based chunked reading..."
        )
      )

      NULL
    }
  )
}

#' @export
#' @rdname resolution-depth
calcGenomeCovChunked <- function(
  pairs_file = NULL, bin_size, min_contacts = 1000, cache = NULL
) {
  # Use cache if provided for fast calculation
  if (!is.null(cache)) {
    if (!inherits(cache, "externalptr")) {
      stop("cache must be an external pointer created by readPairsCache()")
    }

    coverage <- .Call(
      "calculate_coverage_from_cache_c",
      cache,
      as.integer(bin_size),
      as.integer(min_contacts),
      PACKAGE = "gghic"
    )
    return(coverage)
  }

  # Original file-based implementation
  if (is.null(pairs_file)) {
    stop("Either pairs_file or cache must be provided")
  }

  depth <- calcResDepthChunked(
    pairs_file = pairs_file, bin_size = bin_size
  )

  if (is.null(depth) || nrow(depth) == 0) {
    return(0.0)
  }

  total_bins <- nrow(depth)
  bins_above_threshold <- sum(depth$count >= min_contacts, na.rm = TRUE)

  bins_above_threshold / total_bins
}

#' @export
#' @rdname resolution-depth
findOptResChunked <- function(
  pairs_file = NULL, min_bin = 1000, max_bin = 5000000, target_coverage = 0.8,
  min_contacts = 1000, cache = NULL, return_cache = FALSE
) {
  # Use provided cache or create new one
  if (is.null(cache)) {
    if (is.null(pairs_file)) {
      stop("Either pairs_file or cache must be provided")
    }
    pairs_file <- normalizePath(pairs_file, mustWork = TRUE)
    message("Reading positions from file into C memory (one time only)...")
    cache_ptr <- tryCatch(
      {
        .Call("read_pairs_cache_c", pairs_file, PACKAGE = "gghic")
      },
      error = function(e) {
        stop("Failed to read positions from file: ", e$message)
      }
    )
  } else {
    if (!inherits(cache, "externalptr")) {
      stop("cache must be an external pointer created by readPairsCache()")
    }
    cache_ptr <- cache
    message("Using provided cache")
  }

  left <- min_bin
  right <- max_bin
  best_bin <- right

  iteration <- 0
  while (left <= right) {
    iteration <- iteration + 1
    mid <- (left + right) %/% 2

    coverage <- .Call(
      "calculate_coverage_from_cache_c",
      cache_ptr,
      as.integer(mid),
      as.integer(min_contacts),
      PACKAGE = "gghic"
    )

    message(sprintf(
      "Iteration %d: Testing bin size %s bp: %.2f%% coverage",
      iteration, format(mid, big.mark = ","), coverage * 100
    ))

    if (coverage >= target_coverage) {
      best_bin <- mid
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }

  if (return_cache) {
    return(list(bin_size = best_bin, cache = cache_ptr))
  }

  best_bin
}
