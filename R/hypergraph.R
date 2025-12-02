#' Build Hypergraph from Multi-way Contacts
#'
#' Analyze multi-way contacts from Pore-C/HiPore-C data by constructing
#' hypergraph representations where genomic bins are nodes and reads are
#' hyperedges connecting multiple bins.
#' @param pairs A data frame or tibble with columns:
#'   \code{read_name}, \code{chrom1}, \code{pos1}, \code{chrom2}, \code{pos2}.
#'   For large files, use \code{pairs_file} instead.
#' @param pairs_file Character path to pairs file for chunked reading.
#' @param bin_size Integer bin size in base pairs.
#' @param chrom Character chromosome name(s) to analyze (e.g., "chr22" or
#'   c("chr21", "chr22")). Use NULL for genome-wide analysis. Default: NULL.
#' @param min_contacts Integer minimum number of contacts per bin pair.
#'   Filters out sparse pairwise interactions. Default: NULL (no filtering).
#' @param quantile Numeric quantile threshold (0-1) for filtering pairwise
#'   contacts. Keeps only bin pairs above this quantile. Default: 0.85.
#' @param min_multiway Integer minimum number of bins a read must contact
#'   to be included in hypergraph. Default: 3.
#' @param inter_chrom Logical, whether to include inter-chromosomal contacts.
#'   When FALSE, only intra-chromosomal contacts are kept.
#'   When TRUE and multiple
#'   chromosomes are specified in \code{chrom}, includes both intra- and inter-
#'   chromosomal contacts between the specified chromosomes only (not with other
#'   chromosomes in the genome). Default: FALSE (intra-chromosomal only).
#' @return A list with class "hypergraph" containing:
#'   \item{incidence}{Sparse incidence matrix (bins x reads)}
#'   \item{bins}{Integer vector of bin IDs (sorted)}
#'   \item{reads}{Character vector of read names}
#'   \item{contacts_per_read}{Integer vector of contact counts per read}
#'   \item{bin_size}{The bin size used}
#'   \item{chrom}{The chromosome analyzed}
#' @details
#' The function performs the following steps:
#' 1. Bins genomic positions into fixed-size bins
#' 2. Removes duplicate pairwise contacts within each read
#' 3. Filters bin pairs by count quantile or minimum threshold
#' 4. Constructs incidence matrix (bins as rows, reads as columns)
#' 5. Filters reads by minimum multi-way contact degree
#'
#' The resulting hypergraph can be visualized with \code{geom_hypergraph()}.
#' @examples
#' \dontrun{
#' # Build hypergraph from pairs data
#' hg <- buildHypergraph(
#'   pairs = pairs_df,
#'   bin_size = 100000,
#'   chrom = "chr22",
#'   quantile = 0.85,
#'   min_multiway = 3
#' )
#'
#' # From large file
#' hg <- buildHypergraph(
#'   pairs_file = "contact_matrix.txt",
#'   bin_size = 100000,
#'   chrom = "chr22"
#' )
#'
#' # Multi-chromosome with inter-chromosomal contacts
#' hg <- buildHypergraph(
#'   pairs = pairs_df,
#'   bin_size = 100000,
#'   chrom = c("chr21", "chr22"),
#'   inter_chrom = TRUE
#' )
#'
#' # Genome-wide with inter-chromosomal contacts
#' hg <- buildHypergraph(
#'   pairs = pairs_df,
#'   bin_size = 100000,
#'   chrom = NULL,
#'   inter_chrom = TRUE
#' )
#'
#' # Visualize
#' plotHypergraph(hg)
#' }
#' @export
buildHypergraph <- function(
  pairs = NULL,
  pairs_file = NULL,
  bin_size,
  chrom = NULL,
  min_contacts = NULL,
  quantile = 0.85,
  min_multiway = 3,
  inter_chrom = FALSE
) {
  if (is.null(pairs) && is.null(pairs_file)) {
    stop("Either 'pairs' or 'pairs_file' must be provided")
  }

  if (!is.null(pairs) && !is.null(pairs_file)) {
    stop("Provide only one of 'pairs' or 'pairs_file', not both")
  }

  # Load or filter data
  if (!is.null(pairs_file)) {
    if (!is.null(chrom)) {
      message("Reading pairs from file using C implementation...")
      # C function now handles both single and multiple chromosomes
      # When inter_chrom = TRUE, FILTER_CHROM mode ensures both chromosomes
      # are in the specified list (for multiple chromosomes)
      pairs <- .readPairsChromC(pairs_file, chrom, inter_chrom)
    } else {
      # Genome-wide: read intra- or all contacts based on inter_chrom
      if (inter_chrom) {
        message(
          paste0(
            "Reading all contacts (intra- and inter-chromosomal, ",
            "this may take a while)..."
          )
        )
        pairs <- .read_pairs_all_c(pairs_file)
      } else {
        message(
          "Reading all intra-chromosomal contacts (this may take a while)..."
        )
        pairs <- .readPairsAllIntraC(pairs_file)
      }
    }
  } else {
    # Filter to chromosome(s)
    if (!is.null(chrom)) {
      if (length(chrom) == 1) {
        # Single chromosome
        if (inter_chrom) {
          # Keep all pairs involving this chromosome (inter- and intra-)
          pairs <- pairs |>
            dplyr::filter(chrom1 == chrom | chrom2 == chrom)
        } else {
          # Keep only intra-chromosomal pairs for this chromosome
          pairs <- pairs |>
            dplyr::filter(chrom1 == chrom & chrom2 == chrom)
        }
      } else {
        # Multiple chromosomes
        if (inter_chrom) {
          # Keep pairs where BOTH chromosomes are in the specified list
          # (intra- and inter-chromosomal contacts
          # between specified chromosomes only)
          pairs <- pairs |>
            dplyr::filter(chrom1 %in% chrom & chrom2 %in% chrom)
        } else {
          # Keep only intra-chromosomal pairs for these chromosomes
          pairs <- pairs |>
            dplyr::filter(
              chrom1 %in% chrom & chrom2 %in% chrom & chrom1 == chrom2
            )
        }
      }
    } else {
      # Filter to intra-chromosomal only if inter_chrom = FALSE
      if (!inter_chrom) {
        pairs <- pairs |>
          dplyr::filter(chrom1 == chrom2)
      }
      # Otherwise keep all contacts (no filtering)
    }
  }

  if (nrow(pairs) == 0) {
    contact_type <- if (inter_chrom) {
      "contacts"
    } else {
      "intra-chromosomal contacts"
    }
    if (!is.null(chrom)) {
      stop(
        "No ", contact_type, " found for ",
        paste(chrom, collapse = ", ")
      )
    } else {
      stop("No ", contact_type, " found")
    }
  }

  chrom_label <- if (!is.null(chrom)) {
    if (length(chrom) == 1) chrom else paste0(length(chrom), " chromosomes")
  } else {
    "genome-wide"
  }

  message(sprintf(
    "Processing %s contacts (%s)",
    format(nrow(pairs), big.mark = ","), chrom_label
  ))

  # Additional filtering: remove reads with inter-chromosomal contacts
  # This is necessary because a single read may have contacts on different
  # chromosomes across multiple rows, even after filtering pairs
  if (!inter_chrom) {
    # Identify reads that have contacts on multiple chromosomes
    read_chroms <- pairs |>
      dplyr::select(read_name, chrom1, chrom2) |>
      tidyr::pivot_longer(cols = c(chrom1, chrom2), values_to = "chrom") |>
      dplyr::distinct(read_name, chrom) |>
      dplyr::count(read_name, name = "n_chroms")

    # Keep only reads that contact a single chromosome
    reads_to_keep <- read_chroms |>
      dplyr::filter(n_chroms == 1) |>
      dplyr::pull(read_name)

    n_reads_before <- length(unique(pairs$read_name))
    pairs <- pairs |>
      dplyr::filter(read_name %in% reads_to_keep)
    n_reads_after <- length(unique(pairs$read_name))

    if (n_reads_before > n_reads_after) {
      message(sprintf(
        "Removed %d reads with inter-chromosomal contacts (%d reads remaining)",
        n_reads_before - n_reads_after, n_reads_after
      ))
    }

    if (nrow(pairs) == 0) {
      stop("No reads remain after filtering inter-chromosomal contacts")
    }
  }

  # Bin positions - handle chromosome info for multi-chrom
  dat <- pairs |>
    dplyr::mutate(
      bin1 = ceiling(pos1 / bin_size),
      bin2 = ceiling(pos2 / bin_size),
      # Create unique bin IDs across chromosomes
      bin1_id = paste0(chrom1, ":", bin1),
      bin2_id = paste0(chrom2, ":", bin2)
    ) |>
    dplyr::select(read_name, chrom1, bin1, bin1_id, chrom2, bin2, bin2_id) |>
    dplyr::distinct()

  # Count pairwise contacts and filter
  counts_pw <- dat |>
    dplyr::count(chrom1, bin1_id, chrom2, bin2_id, name = "count")

  if (!is.null(min_contacts)) {
    threshold <- min_contacts
  } else {
    threshold <- quantile(counts_pw$count, quantile)
  }

  message(sprintf(
    "Filtering bin pairs with >= %d contacts (%.0f%% quantile)",
    threshold, quantile * 100
  ))

  counts_pw_fil <- counts_pw |>
    dplyr::filter(count >= threshold)

  # Join to get filtered data
  dat_fil <- dat |>
    dplyr::inner_join(
      counts_pw_fil |> dplyr::select(-count),
      by = c("chrom1", "bin1_id", "chrom2", "bin2_id")
    )

  if (nrow(dat_fil) == 0) {
    stop(
      "No contacts remain after filtering. Try lower quantile or min_contacts"
    )
  }

  message(sprintf(
    "Retained %s contacts from %s reads",
    format(nrow(dat_fil), big.mark = ","),
    format(length(unique(dat_fil$read_name)), big.mark = ",")
  ))

  # Create indices - use bin_id for unique identification
  bin1_levels <- unique(dat_fil$bin1_id)
  bin2_levels <- unique(dat_fil$bin2_id)
  all_bin_ids <- unique(c(bin1_levels, bin2_levels))

  # Sort bin IDs by chromosome and position
  bin_info <- tibble::tibble(bin_id = all_bin_ids) |>
    dplyr::mutate(
      chrom = sub(":.*", "", bin_id),
      bin_num = as.integer(sub(".*:", "", bin_id))
    ) |>
    dplyr::arrange(chrom, bin_num)

  all_bin_ids_sorted <- bin_info$bin_id

  read_levels <- unique(dat_fil$read_name)

  dat_indexed <- dat_fil |>
    dplyr::mutate(
      bin1_idx = match(bin1_id, all_bin_ids_sorted),
      bin2_idx = match(bin2_id, all_bin_ids_sorted),
      read_idx = match(read_name, read_levels)
    )

  # Build incidence matrix using Matrix package
  n_bins <- length(all_bin_ids_sorted)
  n_reads <- length(read_levels)

  # Create sparse matrix for bin1
  mat1 <- Matrix::sparseMatrix(
    i = dat_indexed$bin1_idx,
    j = dat_indexed$read_idx,
    x = rep(1, nrow(dat_indexed)),
    dims = c(n_bins, n_reads)
  )

  # Create sparse matrix for bin2
  mat2 <- Matrix::sparseMatrix(
    i = dat_indexed$bin2_idx,
    j = dat_indexed$read_idx,
    x = rep(1, nrow(dat_indexed)),
    dims = c(n_bins, n_reads)
  )

  # Combine (logical OR, then convert to integer)
  incidence <- (mat1 + mat2 > 0) * 1

  # Count contacts per read
  contacts_per_read <- Matrix::colSums(incidence)

  # Filter by minimum multi-way contacts
  keep_reads <- which(contacts_per_read >= min_multiway)

  if (length(keep_reads) == 0) {
    stop(sprintf(
      "No reads with >= %d contacts. Try lower min_multiway", min_multiway
    ))
  }

  incidence_fil <- incidence[, keep_reads, drop = FALSE]

  # Remove empty bins after filtering
  bin_sums <- Matrix::rowSums(incidence_fil)
  keep_bins <- which(bin_sums > 0)

  incidence_final <- incidence_fil[keep_bins, , drop = FALSE]

  message(sprintf(
    "Final hypergraph: %d bins, %d reads (min %d-way contacts)",
    nrow(incidence_final), ncol(incidence_final), min_multiway
  ))

  # Get chromosome info for kept bins
  bin_info_final <- bin_info[keep_bins, ]

  structure(
    list(
      incidence = incidence_final,
      bin_ids = all_bin_ids_sorted[keep_bins],
      bin_info = bin_info_final,
      reads = read_levels[keep_reads],
      contacts_per_read = contacts_per_read[keep_reads],
      bin_size = bin_size,
      chrom = chrom
    ),
    class = "hypergraph"
  )
}

#' Convert Hypergraph to Long Format for Visualization
#'
#' @param hg A hypergraph object from \code{buildHypergraph()}.
#' @param max_reads Integer maximum number of reads to include.
#'   Default: NULL (all).
#' @return A tibble with columns:
#'   \code{read_name}, \code{read_idx}, \code{bin_id}, \code{chrom}, \code{bin},
#'   \code{bin_idx}, \code{n_contacts}.
#' @examples
#' \dontrun{
#' # Load Pore-C pairs data
#' pairs_df <- readr::read_tsv("path/to/porec_pairs.txt")
#'
#' # Build hypergraph from Pore-C data
#' hg <- buildHypergraph(
#'   pairs = pairs_df,
#'   bin_size = 100000,
#'   chrom = "chr22",
#'   quantile = 0.85,
#'   min_multiway = 3
#' )
#'
#' # Convert to tidy format
#' tidy_hg <- tidyHypergraph(hg)
#'
#' # Limit to top 50 reads by contact count
#' tidy_hg_top50 <- tidyHypergraph(hg, max_reads = 50)
#'
#' # Use with ggplot2
#' library(ggplot2)
#' ggplot(tidy_hg_top50, aes(x = read_idx, y = bin_idx, group = read_idx)) +
#'   geom_line(aes(color = n_contacts)) +
#'   geom_point() +
#'   theme_minimal()
#' }
#' @export
tidyHypergraph <- function(hg, max_reads = NULL) {
  if (!inherits(hg, "hypergraph")) {
    stop("Input must be a hypergraph object from buildHypergraph()")
  }

  # Optionally subset reads
  if (!is.null(max_reads) && max_reads < length(hg$reads)) {
    # Sort by number of contacts descending, take top
    order_reads <- order(hg$contacts_per_read, decreasing = TRUE)
    keep_idx <- order_reads[seq_len(max_reads)]

    incidence_sub <- hg$incidence[, keep_idx, drop = FALSE]
    reads_sub <- hg$reads[keep_idx]
    contacts_sub <- hg$contacts_per_read[keep_idx]

    # Remove empty bins
    bin_sums <- Matrix::rowSums(incidence_sub)
    keep_bins <- which(bin_sums > 0)

    incidence_sub <- incidence_sub[keep_bins, , drop = FALSE]
    bin_ids_sub <- hg$bin_ids[keep_bins]
    bin_info_sub <- hg$bin_info[keep_bins, ]
  } else {
    incidence_sub <- hg$incidence
    reads_sub <- hg$reads
    contacts_sub <- hg$contacts_per_read
    bin_ids_sub <- hg$bin_ids
    bin_info_sub <- hg$bin_info
  }

  # Convert to long format
  mat_summary <- Matrix::summary(incidence_sub)

  tibble::tibble(
    bin_idx = mat_summary$i,
    read_idx = mat_summary$j,
    bin_id = bin_ids_sub[bin_idx],
    chrom = bin_info_sub$chrom[bin_idx],
    bin = bin_info_sub$bin_num[bin_idx],
    read_name = reads_sub[read_idx],
    n_contacts = contacts_sub[read_idx]
  ) |>
    dplyr::arrange(chrom, bin, read_idx)
}

#' Plot Hypergraph Visualization
#'
#' Creates a hypergraph visualization showing bins as nodes (circles) on the
#' y-axis and reads as vertical lines connecting multiple bins.
#'
#' @param hg A hypergraph object from \code{buildHypergraph()}.
#' @param chrom Character vector of chromosome names to plot. If NULL (default),
#'   plots all chromosomes in the hypergraph. Use this to subset chromosomes
#'   from a multi-chromosome hypergraph.
#' @param max_reads Integer maximum number of reads to display. Default: 100.
#' @param point_size Numeric size of bin points. Default: 2.
#' @param line_width Numeric width of read connection lines. Default: 0.3.
#' @param line_alpha Numeric transparency of lines. Default: 0.6.
#' @param color_by Character, how to color reads: "n_contacts" (default) or
#'   "none".
#' @param palette Character color palette. Default: "viridis".
#' @param facet_chrom Logical, whether to facet by chromosome for multi-chrom
#'   hypergraphs. Default: TRUE.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' hg <- buildHypergraph(pairs, bin_size = 1e5, chrom = "chr22")
#' plotHypergraph(hg, max_reads = 50)
#'
#' # Multi-chromosome
#' hg <- buildHypergraph(pairs, bin_size = 1e5, chrom = c("chr21", "chr22"))
#' plotHypergraph(hg, max_reads = 100, facet_chrom = TRUE)
#'
#' # Plot only specific chromosomes from multi-chromosome hypergraph
#' plotHypergraph(hg, chrom = "chr21", max_reads = 50)
#' plotHypergraph(hg, chrom = c("chr21", "chr22"), facet_chrom = FALSE)
#' }
#'
#' @export
plotHypergraph <- function(
  hg,
  chrom = NULL,
  max_reads = 100,
  point_size = 2,
  line_width = 0.3,
  line_alpha = 0.6,
  color_by = "n_contacts",
  palette = "viridis",
  facet_chrom = TRUE
) {
  df <- tidyHypergraph(hg, max_reads = max_reads)

  # Filter by chromosome if specified
  if (!is.null(chrom)) {
    available_chroms <- unique(df$chrom)
    invalid_chroms <- setdiff(chrom, available_chroms)

    if (length(invalid_chroms) > 0) {
      stop(
        "Chromosome(s) not found in hypergraph: ",
        paste(invalid_chroms, collapse = ", "),
        "\nAvailable chromosomes: ",
        paste(available_chroms, collapse = ", ")
      )
    }

    df <- df |> dplyr::filter(chrom %in% !!chrom)

    if (nrow(df) == 0) {
      stop("No data remaining after filtering by chromosome")
    }
  }

  # Create a unique x position for each read
  df <- df |>
    dplyr::mutate(
      x = as.numeric(factor(read_idx, levels = unique(read_idx)))
    )

  # Determine if multi-chromosome
  n_chroms <- length(unique(df$chrom))

  # For multi-chromosome, create composite y-axis position
  if (n_chroms > 1 && !facet_chrom) {
    # Sort chromosomes and assign offsets
    # chrom_order <- sort(unique(df$chrom))

    chrom_info <- df |>
      dplyr::group_by(chrom) |>
      dplyr::summarise(
        max_bin = max(bin_idx),
        min_bin = min(bin_idx),
        .groups = "drop"
      ) |>
      dplyr::arrange(chrom)

    # Calculate cumulative offsets
    chrom_info$offset <- 0
    chrom_info$y_start <- 0
    chrom_info$y_end <- chrom_info$max_bin - chrom_info$min_bin + 1

    if (nrow(chrom_info) > 1) {
      for (i in 2:nrow(chrom_info)) {
        chrom_info$offset[i] <- chrom_info$y_end[i - 1] +
          max(1, 0.05 * chrom_info$y_end[i - 1])
        chrom_info$y_start[i] <- chrom_info$offset[i]
        chrom_info$y_end[i] <- chrom_info$y_start[i] +
          (chrom_info$max_bin[i] - chrom_info$min_bin[i] + 1)
      }
    }

    # Apply offsets to create composite y-axis
    df <- df |>
      dplyr::left_join(
        chrom_info |> dplyr::select(chrom, offset, min_bin),
        by = "chrom"
      ) |>
      dplyr::mutate(
        y_composite = bin_idx - min_bin + offset
      )

    y_var <- "y_composite"

    # Create chromosome labels for y-axis
    chrom_breaks <- chrom_info |>
      dplyr::mutate(
        mid_point = (y_start + y_end) / 2
      )
  } else {
    y_var <- "bin_idx"
    chrom_breaks <- NULL
  }

  # Base plot with appropriate y variable
  if (y_var == "y_composite") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_composite))
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = bin_idx))
  }

  # Add background shading for multi-chromosome plots (must be first layer)
  if (n_chroms > 1 && !facet_chrom && !is.null(chrom_breaks)) {
    chrom_rects <- chrom_breaks |>
      dplyr::mutate(
        fill_color = dplyr::if_else(
          dplyr::row_number() %% 2 == 0, "gray95", "white"
        )
      )

    p <- p +
      ggplot2::geom_rect(
        data = chrom_rects,
        ggplot2::aes(
          xmin = -Inf, xmax = Inf,
          ymin = y_start, ymax = y_end,
          fill = fill_color
        ),
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_identity()
  }

  # Add lines connecting bins for each read
  if (color_by == "n_contacts") {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(group = read_idx, color = n_contacts),
        linewidth = line_width, alpha = line_alpha
      ) +
      ggplot2::scale_color_viridis_c(
        option = palette, name = "Contacts\nper read"
      )
  } else {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(group = read_idx),
        linewidth = line_width, alpha = line_alpha, color = "gray50"
      )
  }

  p <- p +
    ggplot2::geom_point(size = point_size, color = "gray30", alpha = 0.8)

  # Labels
  chrom_label <- if (!is.null(hg$chrom)) {
    if (length(hg$chrom) == 1) {
      hg$chrom
    } else {
      paste0(length(hg$chrom), " chromosomes")
    }
  } else {
    "genome-wide"
  }

  # Y-axis label
  y_label <- if (n_chroms == 1) {
    sprintf("Genomic position (%s)", unique(df$chrom))
  } else if (n_chroms > 1 && facet_chrom) {
    "Genomic position"
  } else {
    "Genomic position (by chromosome)"
  }

  p <- p +
    ggplot2::labs(
      title = sprintf(
        "Multi-way Contacts: %s (bin size: %s)",
        chrom_label,
        format(hg$bin_size, big.mark = ",")
      ),
      x = "Read",
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # Add y-axis labels and formatting when not faceting
  if (!facet_chrom) {
    if (n_chroms > 1 && !is.null(chrom_breaks)) {
      # Multi-chromosome: show both bin positions and chromosome labels

      # Add horizontal lines between chromosomes
      if (nrow(chrom_breaks) > 1) {
        p <- p +
          ggplot2::geom_hline(
            data = chrom_breaks[2:nrow(chrom_breaks), ],
            ggplot2::aes(yintercept = y_start),
            linetype = "solid",
            color = "gray40",
            linewidth = 0.8,
            alpha = 0.7
          )
      }

      # Generate bin position breaks for each chromosome
      all_breaks <- c()
      all_labels <- c()

      for (i in seq_len(nrow(chrom_breaks))) {
        chrom_data <- df |> dplyr::filter(chrom == chrom_breaks$chrom[i])
        y_range <- range(chrom_data$y_composite)
        bin_range <- range(chrom_data$bin)

        # Create 3-5 breaks per chromosome
        n_breaks <- min(5, max(3, length(unique(chrom_data$bin)) %/% 10))
        bin_breaks <- seq(bin_range[1], bin_range[2], length.out = n_breaks)

        # Map bin numbers to y_composite positions
        for (b in bin_breaks) {
          closest_row <- chrom_data |>
            dplyr::filter(abs(bin - b) == min(abs(bin - b))) |>
            dplyr::slice(1)
          all_breaks <- c(all_breaks, closest_row$y_composite[1])
          all_labels <- c(all_labels, format(round(b), big.mark = ","))
        }
      }

      # Custom y-axis with bin position labels
      p <- p +
        ggplot2::scale_y_continuous(
          breaks = all_breaks,
          labels = all_labels,
          expand = ggplot2::expansion(mult = c(0.02, 0.02))
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0.02, 0.15))
        )

      # Add chromosome labels on the right side
      p <- p +
        ggplot2::annotate(
          "text",
          x = max(df$x) * 1.05,
          y = chrom_breaks$mid_point,
          label = chrom_breaks$chrom,
          hjust = 0,
          vjust = 0.5,
          fontface = "bold",
          size = 4,
          color = "gray30"
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(5.5, 40, 5.5, 5.5))
    } else if (n_chroms == 1) {
      # Single chromosome: show bin positions
      bin_range <- range(df$bin)
      n_breaks <- min(10, max(5, length(unique(df$bin)) %/% 10))
      bin_breaks <- seq(bin_range[1], bin_range[2], length.out = n_breaks)

      # Map bin numbers to bin_idx positions
      y_breaks <- c()
      y_labels <- c()

      for (b in bin_breaks) {
        closest_row <- df |>
          dplyr::filter(abs(bin - b) == min(abs(bin - b))) |>
          dplyr::slice(1)
        y_breaks <- c(y_breaks, closest_row$bin_idx[1])
        y_labels <- c(y_labels, format(round(b), big.mark = ","))
      }

      p <- p +
        ggplot2::scale_y_continuous(
          breaks = y_breaks,
          labels = y_labels,
          expand = ggplot2::expansion(mult = c(0.02, 0.02))
        )
    }
  }

  # Facet by chromosome if multi-chrom and requested
  if (n_chroms > 1 && facet_chrom) {
    p <- p +
      ggplot2::facet_grid(chrom ~ ., scales = "free_y", space = "free_y") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
  }

  p
}

#' Print method for hypergraph objects
#'
#' @param x A hypergraph object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.hypergraph <- function(x, ...) {
  cat("Hypergraph object\n")
  cat("================\n")

  chrom_label <- if (!is.null(x$chrom)) {
    if (length(x$chrom) == 1) {
      x$chrom
    } else {
      sprintf(
        "%s (%d chromosomes)",
        paste(head(x$chrom, 3), collapse = ", "),
        length(x$chrom)
      )
    }
  } else {
    "genome-wide"
  }

  cat(sprintf("Chromosome(s): %s\n", chrom_label))
  cat(sprintf("Bin size: %s bp\n", format(x$bin_size, big.mark = ",")))
  cat(sprintf("Bins: %d", nrow(x$incidence)))

  if (nrow(x$bin_info) > 0) {
    bin_range <- range(x$bin_info$bin_num)
    cat(sprintf(" (range: %d-%d)\n", bin_range[1], bin_range[2]))
  } else {
    cat("\n")
  }

  cat(sprintf("Reads: %d\n", length(x$reads)))
  cat(sprintf("Total edges: %d\n", Matrix::nnzero(x$incidence)))
  cat(sprintf(
    "Contacts per read: %d-%d (median: %d)\n",
    min(x$contacts_per_read),
    max(x$contacts_per_read),
    median(x$contacts_per_read)
  ))

  invisible(x)
}
