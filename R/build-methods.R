#' Build hypergraph from multi-way contacts
#'
#' @name build
#' @aliases build,MultiWayContacts-method
#'
#' @description
#' Constructs hypergraph from pairs data by binning, filtering, and building
#' sparse incidence matrix.
#'
#' @param x MultiWayContacts object with loaded pairs.
#' @param bin_size Integer. Bin size in base pairs (default: 1000000L).
#' @param quantile Numeric (0-1). Quantile threshold for filtering bin pairs by
#'   contact frequency (default: 0.85).
#' @param min_multiway Integer. Minimum contacts per hyperedge to retain
#'   (default: 2).
#'
#' @return MultiWayContacts with populated hypergraph slots: hypergraph
#'   (incidence matrix), weights (normalizations), bin_info, hyperedge_reads,
#'   multiways, and bin_size.
#'
#' @details
#' Steps: bin positions → filter by quantile → build incidence matrix →
#' deduplicate patterns → filter by min_multiway.
#'
#' Weight normalizations: raw (frequency), log (log-transformed), by_order
#' (normalized within order), minmax (scaled to \[0,1\]).
#'
#' Higher quantile = stricter filtering. Typical: 0.8-0.95.
#'
#' @examples
#' \dontrun{
#' # Default parameters
#' mc <- MultiWayContacts("sample.pairs.gz") |> import() |> build()
#'
#' # Custom parameters
#' mc <- mc |> build(bin_size = 500000L, quantile = 0.9, min_multiway = 3)
#' }
#'
#' @export
methods::setMethod(
  "build", "MultiWayContacts",
  function(x, bin_size = 1000000L, quantile = 0.85, min_multiway = 2) {
    pairs <- x@pairs
    if (is.null(pairs)) {
      stop("Pairs data is not loaded in the MultiWayContacts object.")
    }
    stopifnot(
      "bin_size must be a single positive integer" =
        is.integer(bin_size) & length(bin_size) == 1 & bin_size > 0
    )

    dat_fil <- .binAndFilterContacts(pairs, bin_size, NULL, quantile)
    incidence_res <- .buildIncidenceMatrix(dat_fil)
    dedup_res <- .deduplicateHyperedges(
      incidence_res$incidence, incidence_res$read_levels
    )
    res <- .finalizeHypergraph(
      incidence = dedup_res$incidence,
      hyperedge_weights = dedup_res$hyperedge_weights,
      hyperedge_read_names = dedup_res$hyperedge_read_names,
      bin_ids_sorted = incidence_res$bin_ids_sorted,
      bin_info = incidence_res$bin_info,
      min_multiway = min_multiway,
      bin_size = bin_size,
      chrom = x@focus
    )

    x@hypergraph <- res$incidence
    x@weights <- res$hyperedge_weights
    x@bin_info <- res$bin_info
    x@hyperedge_reads <- res$hyperedge_reads
    x@multiways <- res$multiways_per_hyperedge
    x@bin_size <- as.integer(bin_size)

    x
  }
)
