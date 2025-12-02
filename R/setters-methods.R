#' Set genomic features in ChromatinContacts object
#'
#' @description
#' Assigns genomic features (TADs, compartments, multi_contacts, or tracks) to
#' a `ChromatinContacts` object. The features are automatically confined to
#' overlap with the imported interaction regions.
#'
#' @param x A `ChromatinContacts` object.
#' @param name Character string. Name of the feature slot. Must be one of:
#'   `"compartments"`, `"TADs"`, `"multi_contacts"`, or `"tracks"`.
#' @param ... Additional arguments (not used).
#' @param value A `GRanges` object containing the feature data.
#'
#' @return The modified `ChromatinContacts` object.
#'
#' @details
#' Features are automatically subsetted to only include regions that overlap
#' with the imported interaction data. This ensures consistency when subsetting
#' the `ChromatinContacts` object.
#'
#' **Note:** For tracks, use `GRangesList` instead. If a single `GRanges` is
#' provided for tracks, it will be converted to a `GRangesList` with one element.
#'
#' **Important:** Interactions must be imported before features can be assigned.
#'
#' @examples
#' \dontrun{
#' # Import interactions first
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # Add TADs
#' tads <- rtracklayer::import("TADs.bed")
#' features(cc, "TADs") <- tads
#'
#' # Add compartments
#' features(cc, "compartments") <- compartments_gr
#' }
#'
#' @seealso [features()]
#' @rdname features-set
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

#' Set loop features in ChromatinContacts object
#'
#' @description
#' Assigns chromatin loop data to a `ChromatinContacts` object. Loops are
#' automatically confined to overlap with the imported interaction regions.
#'
#' @param x A `ChromatinContacts` object.
#' @param name Character string. Must be `"loops"`.
#' @param ... Additional arguments (not used).
#' @param value A `GInteractions` object containing loop data.
#'
#' @return The modified `ChromatinContacts` object.
#'
#' @details
#' Loops are stored as `GInteractions` objects representing paired genomic
#' regions. They are automatically subsetted to only include loops that overlap
#' with the imported interaction data.
#'
#' **Important:** Interactions must be imported before loops can be assigned.
#'
#' @examples
#' \dontrun{
#' # Import interactions first
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # Add loops
#' loops <- rtracklayer::import("loops.bedpe") |>
#'   makeGInteractionsFromGRangesPairs()
#' features(cc, "loops") <- loops
#' }
#'
#' @seealso [features()]
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

#' Set track features in ChromatinContacts object
#'
#' @description
#' Assigns genomic signal tracks (e.g., ChIP-seq, ATAC-seq) to a
#' `ChromatinContacts` object. Tracks are automatically confined to overlap with
#' the imported interaction regions.
#'
#' @param x A `ChromatinContacts` object.
#' @param name Character string. Must be `"tracks"`.
#' @param ... Additional arguments. Can include `which` to specify track names.
#' @param value A `GRangesList` object containing track data. Each element
#'   should be a `GRanges` object representing a single track.
#'
#' @return The modified `ChromatinContacts` object.
#'
#' @details
#' Tracks should be provided as a `GRangesList` where each element represents
#' a different signal track. If the list is unnamed, default names
#' (`track1`, `track2`, etc.) will be assigned.
#'
#' **Important:** Interactions must be imported before tracks can be assigned.
#'
#' @examples
#' \dontrun{
#' # Import interactions first
#' cc <- ChromatinContacts("sample.cool") |> import()
#'
#' # Add tracks from BigWig files
#' track1 <- rtracklayer::import("signal1.bw")
#' track2 <- rtracklayer::import("signal2.bw")
#' tracks <- GRangesList(ChIP1 = track1, ChIP2 = track2)
#' features(cc, "tracks") <- tracks_list
#' }
#'
#' @seealso [features()]
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
