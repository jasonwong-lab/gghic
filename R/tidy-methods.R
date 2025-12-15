#' Tidy hypergraph to long format
#'
#' @name tidy
#' @aliases tidy,MultiWayContacts-method
#'
#' @description
#' Converts hypergraph incidence matrix to long-format data frame for analysis
#' and visualization.
#'
#' @param x MultiWayContacts with built hypergraph.
#' @param max_hyperedges Integer or NULL. Maximum hyperedges to include
#'   (default: NULL for all). Selects top by contact order.
#' @param weight_normalization Character. Weight scheme: `"none"` (raw),
#'   `"log"`, `"by_order"`, or `"minmax"` (default: `"none"`).
#'
#' @return MultiWayContacts with `tidied_hypergraph` slot containing: bin_idx,
#'   hyperedge_idx, bin_id, chrom, bin, n_multiways, weight.
#'
#' @details
#' Converts sparse matrix to long format. When `max_hyperedges` is set, ranks
#' by contact order and retains top N. Unused bins are removed.
#'
#' @examples
#' \dontrun{
#' # Tidy all with raw weights
#' mc <- MultiWayContacts("sample.pairs.gz") |> import() |> build() |> tidy()
#'
#' # Top 100 with log weights
#' mc <- mc |> tidy(max_hyperedges = 100, weight_normalization = "log")
#' }
#'
#' @export
methods::setMethod(
  "tidy", "MultiWayContacts",
  function(x, max_hyperedges = NULL, weight_normalization = "none") {
    valid_norms <- c("none", "log", "by_order", "minmax")
    if (!weight_normalization %in% valid_norms) {
      stop(
        "weight_normalization must be one of: ",
        paste(valid_norms, collapse = ", ")
      )
    }

    weight_col <- switch(weight_normalization,
      none = "raw",
      log = "log",
      by_order = "by_order",
      minmax = "minmax"
    )

    weights_all <- x@weights[[weight_col]]

    if (!is.null(max_hyperedges) && max_hyperedges < ncol(x@hypergraph)) {
      # Sort by number of contacts descending, take top
      order_hyperedges <- order(x@multiways, decreasing = TRUE)
      keep_idx <- order_hyperedges[seq_len(max_hyperedges)]

      incidence_sub <- x@hypergraph[, keep_idx, drop = FALSE]
      weights_sub <- weights_all[keep_idx]
      contacts_sub <- x@multiways[keep_idx]

      bin_sums <- Matrix::rowSums(incidence_sub)
      keep_bins <- which(bin_sums > 0)

      incidence_sub <- incidence_sub[keep_bins, , drop = FALSE]
      bin_info_sub <- x@bin_info[keep_bins, ]
      bin_ids_sub <- bin_info_sub$bin_id
    } else {
      incidence_sub <- x@hypergraph
      weights_sub <- weights_all
      contacts_sub <- x@multiways
      bin_info_sub <- x@bin_info
      bin_ids_sub <- bin_info_sub$bin_id
    }

    mat_summary <- Matrix::summary(incidence_sub)

    x@tidied_hypergraph <- tibble::tibble(
      bin_idx = mat_summary$i,
      hyperedge_idx = mat_summary$j,
      bin_id = bin_ids_sub[bin_idx],
      chrom = bin_info_sub$chrom[bin_idx],
      bin = bin_info_sub$bin_num[bin_idx],
      n_multiways = contacts_sub[hyperedge_idx],
      weight = weights_sub[hyperedge_idx]
    ) |>
      dplyr::arrange(chrom, bin, hyperedge_idx)

    x
  }
)
