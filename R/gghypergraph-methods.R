#' Visualize multi-way chromatin contacts as hypergraph network
#'
#' @name gghypergraph
#' @aliases gghypergraph,MultiWayContacts-method
#'
#' @description
#' Creates an intuitive ggplot2 visualization of multi-way chromatin contacts
#' from Pore-C or similar long-read technologies. Displays genomic bins as
#' points and multi-way contacts as connecting lines (hyperedges), with line
#' width and color encoding contact properties.
#'
#' @param x MultiWayContacts object with selected hyperedges (after running
#'   [select()]).
#' @param point_size Numeric. Size of genomic bin points (default: 2).
#' @param line_width Numeric. Base line width for hyperedges. Actual width
#'   scaled by hyperedge weight (default: 0.3).
#' @param line_alpha Numeric between 0 and 1. Line transparency to reduce
#'   overplotting (default: 0.6).
#' @param color_by Character. Variable for line coloring:
#'   * `"n_multiways"`: Color by contact order (2-way, 3-way, etc.) (default)
#'   * Other values: Uniform gray coloring
#' @param palette Character. Color palette name:
#'   * Viridis options: `"viridis"`, `"magma"`, `"plasma"`, `"inferno"`,
#'     `"cividis"`
#'   * ColorBrewer palettes: `"YlOrRd"`, `"Blues"`, etc.
#'   * Default: `"viridis"`
#' @param facet_chrom Logical. Layout strategy:
#'   * `TRUE`: Separate facet panel per chromosome (default)
#'   * `FALSE`: Composite y-axis with all chromosomes stacked
#' @param weight_normalization Character. Deprecated parameter (weights now set
#'   in [tidy()] method).
#'
#' @return A ggplot2 object that can be further customized with additional
#'   ggplot2 layers and themes.
#'
#' @details
#' ## Visualization structure
#' * **X-axis**: Hyperedge index (each read/contact)
#' * **Y-axis**: Genomic bin position (genomic coordinate)
#' * **Points**: Individual genomic bins
#' * **Lines**: Connect all bins in a single multi-way contact
#' * **Line width**: Scaled by hyperedge weight
#' * **Line color**: Contact order (number of fragments per read)
#'
#' ## Interpretation
#' * Vertical lines: Bins contacted by same read
#' * Line thickness: Relative importance/weight
#' * Color gradient: Contact complexity (darker = higher-order)
#' * Patterns reveal genome organization and contact preferences
#'
#' ## Prerequisites
#' Must run complete workflow: `import()` → `build()` → `tidy()` → `select()`
#' before visualization.
#'
#' @examples
#' \dontrun{
#' # Complete workflow
#' mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1") |>
#'   import() |>
#'   build(bin_size = 1000000L, quantile = 0.85) |>
#'   tidy(weight_normalization = "log") |>
#'   select(n_intra = 10, n_inter = 5)
#'
#' # Basic visualization
#' gghypergraph(mc)
#'
#' # Composite layout with custom colors
#' gghypergraph(mc, facet_chrom = FALSE, palette = "magma")
#'
#' # Emphasize line patterns
#' gghypergraph(mc, line_width = 0.5, line_alpha = 0.8, point_size = 3)
#'
#' # Further customization
#' gghypergraph(mc) +
#'   theme_minimal() +
#'   labs(title = "Multi-way Contacts on Chromosome 1")
#' }
#'
#' @seealso [MultiWayContacts()], [select()], [tidy()], [build()]
#' @export
methods::setMethod(
  "gghypergraph", "MultiWayContacts",
  function(
    x, point_size = 2, line_width = 0.3, line_alpha = 0.6,
    color_by = "n_multiways", palette = "viridis", facet_chrom = TRUE,
    weight_normalization = "none"
  ) {
    if (is.null(x@select_hypergraph)) {
      warning(
        "The MultiWayContacts object has no selected hyperedges. ",
        "Please run the select() method first. ",
        "gghypergraph() will plot all tidied hyperedges."
      )
      df <- x@tidied_hypergraph
    } else {
      df <- x@select_hypergraph
    }

    df <- df |>
      dplyr::mutate(
        bin_idx = as.integer(factor(bin_id, levels = unique(bin_id))),
        x = as.integer(factor(hyperedge_idx, levels = unique(hyperedge_idx)))
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
          n_bins = dplyr::n_distinct(bin_idx),
          .groups = "drop"
        ) |>
        dplyr::arrange(chrom)

      # Calculate cumulative offsets based on number of bins
      # Add spacing between chromosomes
      gap_size <- 1
      chrom_info$offset <- 0
      chrom_info$y_start <- 1
      chrom_info$y_end <- chrom_info$n_bins

      if (nrow(chrom_info) > 1) {
        for (i in 2:nrow(chrom_info)) {
          chrom_info$offset[i] <- chrom_info$y_end[i - 1] + gap_size
          chrom_info$y_start[i] <- chrom_info$offset[i] + 1
          chrom_info$y_end[i] <- chrom_info$y_start[i] + chrom_info$n_bins[i] - 1
        }
      }

      # Apply offsets to create composite y-axis
      df <- df |>
        dplyr::left_join(
          chrom_info |> dplyr::select(chrom, offset, min_bin),
          by = "chrom"
        ) |>
        dplyr::mutate(
          y_composite = bin_idx - min_bin + offset + 1
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
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y_composite)) +
        ggplot2::scale_y_continuous(
          breaks = df$y_composite,
          labels = df$bin,
          expand = ggplot2::expansion(mult = c(0.02, 0.02))
        )
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = bin_idx)) +
        ggplot2::scale_y_continuous(
          breaks = df$bin_idx,
          labels = df$bin,
          expand = ggplot2::expansion(mult = c(0.02, 0.02))
        )
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
            xmin = -Inf, xmax = Inf, ymin = y_start, ymax = y_end,
            fill = fill_color
          ),
          alpha = 0.3, inherit.aes = FALSE
        ) +
        ggplot2::scale_fill_identity()
    }

    weight_legend <- "Hyperedge weight"

    if (color_by == "n_multiways") {
      # Check if palette is a viridis option
      viridis_options <- c(
        "viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket",
        "turbo", "A", "B", "C", "D", "E", "F", "G", "H"
      )

      if (palette %in% viridis_options) {
        # Use viridis scale (binned for integer values)
        p <- p +
          ggplot2::geom_line(
            ggplot2::aes(
              group = hyperedge_idx, color = n_multiways, linewidth = weight
            ),
            alpha = line_alpha
          ) +
          ggplot2::scale_color_viridis_b(option = palette, name = "Multiways") +
          ggplot2::geom_point(
            ggplot2::aes(color = n_multiways),
            size = point_size, alpha = 0.8
          )
      } else {
        p <- p +
          ggplot2::geom_line(
            ggplot2::aes(
              group = hyperedge_idx, color = n_multiways, linewidth = weight
            ),
            alpha = line_alpha
          ) +
          ggplot2::scale_color_fermenter(
            palette = palette, name = "Multiways", direction = 1
          ) +
          ggplot2::geom_point(
            ggplot2::aes(color = n_multiways),
            size = point_size, alpha = 0.8
          )
      }
    } else {
      p <- p +
        ggplot2::geom_line(
          ggplot2::aes(group = hyperedge_idx, linewidth = weight),
          alpha = line_alpha, color = "gray50"
        ) +
        ggplot2::geom_point(size = point_size, color = "gray30", alpha = 0.8)
    }

    p <- p +
      ggplot2::scale_linewidth_continuous(
        name = weight_legend,
        range = c(line_width * 0.5, line_width * 3)
      )

    # Labels
    chrom_label <- if (!is.null(x@focus)) {
      if (length(x@focus) == 1) {
        x@focus
      } else {
        paste0(length(x@focus), " chromosomes")
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
          format(x@bin_size, big.mark = ",")
        ),
        x = "Hyperedge",
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

        # Custom y-axis with bin position labels
        p <- p +
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
)
