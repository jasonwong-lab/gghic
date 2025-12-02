#' @rdname geom_hypergraph
#' @format NULL
#' @usage NULL
#' @export
GeomHypergraphLine <- ggplot2::ggproto(
  "GeomHypergraphLine",
  ggplot2::Geom,
  required_aes = c("x", "y", "group"),
  optional_aes = c("colour"),
  draw_key = ggplot2::draw_key_path,
  default_aes = ggplot2::aes(colour = "gray50", linewidth = 0.3, alpha = 0.6),
  draw_panel = function(
    data, panel_params, coord, line_width, line_alpha, colour_by
  ) {
    coords <- coord$transform(data, panel_params)

    grid::polylineGrob(
      x = coords$x,
      y = coords$y,
      id = coords$group,
      gp = grid::gpar(
        col = coords$colour,
        lwd = line_width * ggplot2::.pt,
        alpha = line_alpha
      )
    )
  }
)

#' @rdname geom_hypergraph
#' @format NULL
#' @usage NULL
#' @export
GeomHypergraphPoint <- ggplot2::ggproto(
  "GeomHypergraphPoint",
  ggplot2::Geom,
  required_aes = c("x", "y"),
  draw_key = ggplot2::draw_key_point,
  default_aes = ggplot2::aes(
    colour = "gray30", size = 2, alpha = 0.8, shape = 19
  ),
  draw_panel = function(data, panel_params, coord, point_size, point_alpha) {
    coords <- coord$transform(data, panel_params)

    grid::pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = coords$shape,
      gp = grid::gpar(
        col = coords$colour,
        alpha = point_alpha
      ),
      size = grid::unit(point_size, "char")
    )
  }
)

#' Visualize Multi-way Contacts as Hypergraph
#'
#' Creates a hypergraph visualization for multi-way chromatin contacts where
#' genomic bins are displayed as points (nodes) and reads are shown as vertical
#' lines (hyperedges) connecting multiple bins.
#'
#' @param mapping Set of aesthetic mappings created by `aes()`. Only
#'   `x`, `y`, and `group` are used.
#' @param data A data frame from `hypergraph_to_tidy()`.
#' @param stat The statistical transformation to use on the data.
#'   Default: "identity".
#' @param position Position adjustment. Default: "identity".
#' @param na.rm If `FALSE` (default),
#'   missing values are removed with a warning.
#' @param show.legend Logical. Should this layer be included in legends?
#' @param inherit.aes If `FALSE`, overrides the default aesthetics.
#' @param line_width Numeric line width for hyperedges. Default: 0.3.
#' @param line_alpha Numeric transparency for lines. Default: 0.6.
#' @param point_size Numeric size of bin points. Default: 2.
#' @param point_alpha Numeric transparency for points. Default: 0.8.
#' @param colour_by Character: "n_contacts" or "none". Default: "n_contacts".
#' @param palette Character: ggplot2 colour palette. Default: "viridis".
#' @param ... Additional arguments passed to layer.
#'
#' @details
#' This geom creates a specialized visualization for multi-way chromatin
#' contacts from Pore-C/HiPore-C data. The visualization shows:
#'   - **Y-axis**: Genomic bins (sorted by position)
#'   - **X-axis**: Individual reads (ordered by first appearance)
#'   - **Lines**: Connect all bins contacted by a single read
#'   - **Points**: Individual bin positions
#'   - **colour**: Number of contacts per read (optional)
#'
#' The visualization is similar to a Sankey diagram but specialized for
#' multi-way genomic interactions.
#'
#' @section Aesthetics:
#' `geom_hypergraph()` understands the following aesthetics:
#'     - `x` (required): Read index or position
#'     - `y` (required): Bin position
#'     - `group` (required): Read identifier
#'     - `colour` (optional): Number of contacts per read
#'
#' @examples
#' \dontrun{
#' # Build hypergraph
#' hg <- build_hypergraph(
#'   pairs = pairs_df,
#'   bin_size = 100000,
#'   chrom = "chr22",
#'   quantile = 0.85,
#'   min_multiway = 3
#' )
#'
#' # Convert to tidy format
#' df <- hypergraph_to_tidy(hg, max_reads = 100)
#'
#' # Create visualization
#' ggplot(df, aes(x = read_idx, y = bin, group = read_idx)) +
#'   geom_hypergraph(aes(colour = n_contacts)) +
#'   labs(title = "Multi-way Contacts on chr22") +
#'   theme_minimal()
#'
#' # Without colouring
#' ggplot(df, aes(x = read_idx, y = bin, group = read_idx)) +
#'   geom_hypergraph(colour_by = "none") +
#'   theme_minimal()
#' }
#'
#' @export
geom_hypergraph <- function(
  mapping = NULL, data = NULL, stat = "identity", position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, line_width = 0.3,
  line_alpha = 0.6, point_size = 2, point_alpha = 0.8, colour_by = "n_contacts",
  palette = "viridis", ...
) {
  list(
    ggplot2::layer(
      geom = GeomHypergraphLine, mapping = mapping, data = data, stat = stat,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm, line_width = line_width, line_alpha = line_alpha,
        colour_by = colour_by, ...
      )
    ),
    ggplot2::layer(
      geom = GeomHypergraphPoint, mapping = mapping, data = data, stat = stat,
      position = position, show.legend = FALSE, inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm, point_size = point_size, point_alpha = point_alpha, ...
      )
    )
  )
}
