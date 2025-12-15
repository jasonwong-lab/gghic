#' Select top hyperedges for visualization
#'
#' @name select
#' @aliases select,MultiWayContacts-method
#'
#' @description
#' Selects top-weighted intra- and inter-chromosomal hyperedges per chromosome
#' for visualization.
#'
#' @param x MultiWayContacts with tidied hypergraph.
#' @param n_intra Integer. Top intra-chromosomal hyperedges per chromosome
#'   (default: 5).
#' @param n_inter Integer. Top inter-chromosomal hyperedges per chromosome
#'   (default: 5). Quota increases if fewer intra hyperedges available.
#' @param n_multiways_filter Integer vector or NULL. Filter by specific contact
#'   orders, e.g., `c(3,4)` for 3-way and 4-way only (default: NULL).
#' @param chroms Character or NULL. Filter by chromosome(s) (default: NULL).
#' @param append Logical. Append to existing selection (TRUE) or replace
#'   (FALSE, default: TRUE).
#' @param ... Additional arguments (not used).
#'
#' @return MultiWayContacts with `select_hypergraph` slot containing selected
#'   hyperedges plus `type` (intra/inter) and `source_chrom` columns.
#'
#' @details
#' Per chromosome: select top `n_intra` intra-chromosomal and `n_inter`
#' inter-chromosomal hyperedges by weight. If insufficient intra hyperedges,
#' quota transferred to inter. Hyperedges deduplicated across chromosomes.
#'
#' @examples
#' \dontrun{
#' # Default: top 5 intra + 5 inter per chromosome
#' mc <- MultiWayContacts("sample.pairs.gz") |>
#'   import() |> build() |> tidy() |> select()
#'
#' # More hyperedges
#' mc <- mc |> select(n_intra = 10, n_inter = 10, append = FALSE)
#'
#' # Filter by contact order
#' mc <- mc |> select(n_multiways_filter = c(3, 4), append = FALSE)
#'
#' # Specific chromosomes
#' mc <- mc |> select(chroms = c("chr1", "chr2"), append = FALSE)
#' }
#'
#' @export
methods::setMethod(
  "select", "MultiWayContacts",
  function(
    x, n_intra = 5, n_inter = 5, n_multiways_filter = NULL, chroms = NULL, append = TRUE, ...
  ) {
    if (is.null(x@tidied_hypergraph)) {
      stop(
        "The MultiWayContacts object has not been tidied yet. Please run the tidy() method first."
      )
    }
    hg_tidy <- x@tidied_hypergraph

    if (!is.null(n_multiways_filter)) {
      hg_tidy <- hg_tidy |>
        dplyr::filter(n_multiways %in% n_multiways_filter)

      if (nrow(hg_tidy) == 0) {
        stop(
          "No hyperedges with the specified contact orders: ",
          paste(n_multiways_filter, collapse = ", ")
        )
      }
    }

    hyperedge_types <- hg_tidy |>
      dplyr::group_by(hyperedge_idx) |>
      dplyr::summarise(
        chroms = list(unique(chrom)),
        n_chroms = dplyr::n_distinct(chrom),
        weight = dplyr::first(weight),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        type = dplyr::if_else(n_chroms == 1, "intra", "inter")
      )

    # Get all chromosomes involved in each hyperedge
    hyperedge_chroms <- hg_tidy |>
      dplyr::select(hyperedge_idx, chrom) |>
      dplyr::distinct() |>
      dplyr::group_by(hyperedge_idx) |>
      dplyr::summarise(chroms = list(unique(chrom)), .groups = "drop")

    # Join type information
    hyperedge_types <- hyperedge_types |>
      dplyr::select(hyperedge_idx, type, weight)

    # Filter by specified chromosomes if provided
    if (!is.null(chroms)) {
      # Get hyperedges that involve any of the specified chromosomes
      hyperedges_to_keep <- hyperedge_chroms |>
        dplyr::filter(purrr::map_lgl(chroms, ~ any(.x %in% !!chroms))) |>
        dplyr::pull(hyperedge_idx)

      if (length(hyperedges_to_keep) == 0) {
        stop(
          "No hyperedges found involving the specified chromosomes: ",
          paste(chroms, collapse = ", ")
        )
      }

      hg_tidy <- hg_tidy |>
        dplyr::filter(hyperedge_idx %in% hyperedges_to_keep)

      hyperedge_types <- hyperedge_types |>
        dplyr::filter(hyperedge_idx %in% hyperedges_to_keep)

      hyperedge_chroms <- hyperedge_chroms |>
        dplyr::filter(hyperedge_idx %in% hyperedges_to_keep)
    }

    # For each chromosome, select top hyperedges
    all_chroms <- if (!is.null(chroms)) chroms else unique(hg_tidy$chrom)

    selected_hyperedges <- purrr::map_df(all_chroms, function(chr) {
      # Get hyperedges involving this chromosome
      hyperedges_chr <- hyperedge_chroms |>
        dplyr::filter(purrr::map_lgl(chroms, ~ chr %in% .x)) |>
        dplyr::pull(hyperedge_idx)

      # Get type and weight info
      hyperedges_info <- hyperedge_types |>
        dplyr::filter(hyperedge_idx %in% hyperedges_chr)

      # Select top intra-chromosomal
      intra_selected <- hyperedges_info |>
        dplyr::filter(type == "intra") |>
        dplyr::arrange(dplyr::desc(weight)) |>
        dplyr::slice_head(n = n_intra)

      n_intra_found <- nrow(intra_selected)
      n_inter_needed <- n_inter + (n_intra - n_intra_found)

      # Select top inter-chromosomal
      inter_selected <- hyperedges_info |>
        dplyr::filter(type == "inter") |>
        dplyr::arrange(dplyr::desc(weight)) |>
        dplyr::slice_head(n = n_inter_needed)

      # Combine
      dplyr::bind_rows(intra_selected, inter_selected) |>
        dplyr::mutate(source_chrom = chr)
    })

    # Remove duplicates (hyperedges selected for multiple chromosomes)
    selected_hyperedges <- selected_hyperedges |>
      dplyr::arrange(source_chrom, dplyr::desc(weight)) |>
      dplyr::distinct(hyperedge_idx, .keep_all = TRUE)

    # Filter original data
    res <- hg_tidy |>
      dplyr::filter(hyperedge_idx %in% selected_hyperedges$hyperedge_idx) |>
      dplyr::left_join(
        selected_hyperedges |> dplyr::select(hyperedge_idx, type, source_chrom),
        by = "hyperedge_idx"
      )

    if (append && !is.null(x@select_hypergraph)) {
      x@select_hypergraph <- dplyr::bind_rows(res, x@select_hypergraph) |>
        dplyr::distinct(hyperedge_idx, bin_idx, .keep_all = TRUE) |>
        dplyr::arrange(chrom, bin, hyperedge_idx)
    } else {
      x@select_hypergraph <- res
    }

    x
  }
)
