#' Set genomic features in ChromatinContacts object
#'
#' @name features-set
#' @aliases features<-,ChromatinContacts,character,GRanges-method
#' @description
#' Assigns or updates genomic features (TADs, loops, compartments, tracks) in
#' a ChromatinContacts object. Features are automatically subsetted to overlap
#' with loaded interactions, ensuring only relevant data is retained.
#'
#' @param x ChromatinContacts object with imported interaction data.
#' @param name Character. Feature slot name:
#'   * `"compartments"`: A/B compartment annotations (GRanges)
#'   * `"TADs"`: Topologically associating domains (GRanges)
#'   * `"multi_contacts"`: Multi-way contact regions (GRanges)
#'   * `"loops"`: Chromatin loop calls (GInteractions)
#'   * `"tracks"`: Genomic signal tracks (GRangesList)
#' @param ... Additional arguments (not currently used).
#' @param value Feature data with appropriate type:
#'   * GRanges: for compartments, TADs, multi_contacts
#'   * GInteractions: for loops
#'   * GRangesList or list: for tracks (will be converted to GRangesList)
#'
#' @return ChromatinContacts object with updated features slot. The modified
#'   object is returned invisibly to support piping.
#'
#' @details
#' ## Requirements
#' Interactions must be imported before assigning features. Use
#' `import()` first.
#'
#' ## Automatic subsetting
#' Features are automatically filtered to retain only elements overlapping the
#' loaded Hi-C interaction data. This prevents visualization errors and reduces
#' memory usage.
#'
#' ## Track naming
#' When providing tracks as GRangesList, ensure they are named for proper
#' labeling in plots. Unnamed tracks receive default names (`track1`, `track2`,
#' etc.).
#'
#' @examples
#' \dontrun{
#' # Load Hi-C data first
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # Add TADs from BED file
#' features(cc, "TADs") <- rtracklayer::import("TADs.bed")
#'
#' # Add chromatin loops
#' features(cc, "loops") <- loops_ginteractions
#'
#' # Add multiple ChIP-seq tracks
#' features(cc, "tracks") <- GRangesList(
#'   H3K27ac = rtracklayer::import("H3K27ac.bw"),
#'   H3K4me3 = rtracklayer::import("H3K4me3.bw")
#' )
#'
#' # Add compartments
#' features(cc, "compartments") <- compartment_gr
#'
#' # Chain operations
#' cc <- ChromatinContacts("sample.cool") |>
#'   import() |>
#'   `features<-`("TADs", value = tads) |>
#'   `features<-`("loops", value = loops)
#' }
#'
#' @seealso [features()], [ChromatinContacts()], [import-ChromatinContacts]
#' @export
methods::setMethod(
  "features<-", c("ChromatinContacts", "character", "GRanges"),
  function(x, name, ..., value) {
    stopifnot(
      "name must be one of 'compartments', 'TADs', 'multi_contacts', 'tracks'" =
        name %in% c("compartments", "TADs", "multi_contacts", "tracks")
    )
    stopifnot(
      "interactions haven't been imported yet" = !is.null(x@interactions)
    )

    value <- .confineGRanges(value, x)

    if (name == "tracks") {
      warning(
        "Assigning a single GRanges to 'tracks' slot. Converting to ",
        "GRangesList with one element."
      )
      value <- GenomicRanges::GRangesList(value)
      names(value) <- "track1"
    }

    methods::slot(x, name) <- value

    x
  }
)

#' @rdname features-set
#' @export
methods::setMethod(
  "features<-", c("ChromatinContacts", "character", "GInteractions"),
  function(x, name, ..., value) {
    stopifnot("name must be 'loops'" = name %in% c("loops"))
    stopifnot(
      "interactions haven't been imported yet" = !is.null(x@interactions)
    )

    methods::slot(x, name) <- .confineGInteractions(value, x)

    x
  }
)

#' @rdname features-set
#' @export
methods::setMethod(
  "features<-", c("ChromatinContacts", "character", "GRangesList"),
  function(x, name, ..., value) {
    stopifnot(
      "name must be 'tracks'" = name %in% c("tracks")
    )
    stopifnot(
      "interactions haven't been imported yet" = !is.null(x@interactions)
    )

    grl <- purrr::map(value, ~ .confineGRanges(.x, x)) |>
      GenomicRanges::GRangesList()
    if (is.null(names(grl))) {
      warning("tracks GRangesList has no names, setting default names.")
      names(grl) <- paste0("track", seq_along(grl))
    }
    methods::slot(x, name) <- grl

    x
  }
)

#' @rdname features-set
#' @export
methods::setMethod(
  "features<-", c("ChromatinContacts", "character", "list"),
  function(x, name, ..., value) {
    stopifnot(
      "name must be 'tracks'" = name %in% c("tracks")
    )
    stopifnot(
      "interactions haven't been imported yet" = !is.null(x@interactions)
    )
    warning(
      "Converting list to GRangesList for 'tracks' slot assignment."
    )

    dots <- list(...)
    which <- if ("which" %in% names(dots)) dots$which else NULL

    grl <- purrr::map(value, ~ .confineGRanges(.x, x)) |>
      GenomicRanges::GRangesList()
    if (is.null(names(grl)) && is.null(which)) {
      warning("tracks GRangesList has no names, setting default names.")
      names(grl) <- paste0("track", seq_along(grl))
    } else if (!is.null(which)) {
      names(grl) <- which
    }
    methods::slot(x, name) <- grl

    x
  }
)
