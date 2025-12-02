#' Get resolution of ChromatinContacts object
#'
#' @description
#' Retrieves the resolution (bin size) in base pairs from a `ChromatinContacts`
#' object.
#'
#' @param x A `ChromatinContacts` object.
#'
#' @return Integer representing resolution in base pairs, or NULL if not set.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.mcool", resolution = 10000)
#' resolution(cc) # Returns 10000
#' }
#'
#' @export
methods::setMethod("resolution", "ChromatinContacts", function(x) x@resolution)

#' Get focus regions of ChromatinContacts object
#'
#' @description
#' Retrieves the focus regions (genomic regions of interest) from a
#' `ChromatinContacts` object.
#'
#' @param x A `ChromatinContacts` object.
#'
#' @return GInteractions object representing the focus regions, or NULL if
#'   genome-wide.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool", focus = "chr1")
#' focus(cc)
#' }
#'
#' @export
methods::setMethod("focus", "ChromatinContacts", function(x) x@focus)

#' Get sequence information from ChromatinContacts object
#'
#' @description
#' Retrieves the sequence information (chromosome names, lengths) from a
#' `ChromatinContacts` object.
#'
#' @param x A `ChromatinContacts` object.
#'
#' @return A `Seqinfo` object containing sequence/chromosome information.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool")
#' seqinfo(cc)
#' }
#'
#' @export
methods::setMethod("seqinfo", "ChromatinContacts", function(x) {
  if (is.null(x@seqinfo)) {
    .getSeqinfo(x@cooler_path, x@resolution)
  }

  x@seqinfo
})

#' Get specific genomic features from ChromatinContacts object
#'
#' @description
#' Retrieves a specific type of genomic feature (TADs, loops, compartments, or
#' multi-contacts) from a `ChromatinContacts` object.
#'
#' @param x A `ChromatinContacts` object.
#' @param name Character string. Name of the feature to retrieve. Must be one
#'   of: `"compartments"`, `"TADs"`, `"loops"`, or `"multi_contacts"`.
#'
#' @return The requested feature as GRanges (for TADs, compartments,
#'   multi_contacts) or GInteractions (for loops).
#'
#' @examples
#' \dontrun{
#' # Get TADs
#' tads <- features(cc, "TADs")
#'
#' # Get loops
#' loops <- features(cc, "loops")
#' }
#'
#' @export
methods::setMethod(
  "features", c("ChromatinContacts", "character"),
  function(x, name) {
    stopifnot(
      "name must be one of 'compartments', 'TADs', 'loops', 'multi_contacts'" =
        name %in% c("compartments", "TADs", "loops", "multi_contacts")
    )
    methods::slot(x, name)
  }
)

#' Get all genomic features from ChromatinContacts object
#'
#' @description
#' Retrieves all genomic features (TADs, loops, compartments, multi-contacts)
#' as a SimpleList.
#'
#' @param x A `ChromatinContacts` object.
#'
#' @return A `SimpleList` containing all features: compartments, TADs, loops,
#'   and multi_contacts.
#'
#' @examples
#' \dontrun{
#' # Get all features
#' all_features <- features(cc)
#' all_features$TADs
#' all_features$loops
#' }
#'
#' @rdname features
#' @export
methods::setMethod("features", c("ChromatinContacts", "missing"), function(x) {
  S4Vectors::SimpleList(list(
    compartments = x@compartments,
    TADs = x@TADs,
    loops = x@loops,
    multi_contacts = x@multi_contacts
  ))
})

#' Get interaction data from ChromatinContacts object
#'
#' @description
#' Retrieves the imported Hi-C interaction data as a `GInteractions` object.
#'
#' @param x A `ChromatinContacts` object.
#'
#' @return A `GInteractions` object containing the Hi-C interactions, or NULL
#'   if data has not been imported yet.
#'
#' @details
#' The interactions must be loaded first using `import()`. The returned
#' `GInteractions` object contains bin coordinates and metadata columns such as
#' `count` and `balanced`.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool") |> import()
#' gis <- interactions(cc)
#' }
#'
#' @export
methods::setMethod(
  "interactions", "ChromatinContacts", function(x) x@interactions
)
