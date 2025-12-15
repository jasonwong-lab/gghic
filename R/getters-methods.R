#' Get resolution of ChromatinContacts object
#'
#' @name resolution
#' @aliases resolution,ChromatinContacts-method
#' @description
#' Retrieves bin size in base pairs.
#'
#' @param x ChromatinContacts object.
#'
#' @return Integer (resolution in bp) or NULL.
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.mcool", resolution = 10000L)
#' resolution(cc) # Returns 10000
#' }
#'
#' @export
methods::setMethod("resolution", "ChromatinContacts", function(x) x@resolution)

#' Get focus regions
#'
#' @name focus
#' @aliases focus,ChromatinContacts-method
#' @description
#' Retrieves focus regions from ChromatinContacts object.
#'
#' @param x ChromatinContacts object.
#'
#' @return GInteractions (focus regions) or NULL (genome-wide).
#'
#' @examples
#' \dontrun{
#' cc <- ChromatinContacts("sample.cool", focus = "chr1")
#' focus(cc)
#' }
#'
#' @export
methods::setMethod("focus", "ChromatinContacts", function(x) x@focus)

#' Get sequence information
#'
#' @description
#' Retrieves chromosome/sequence information.
#'
#' @param x ChromatinContacts object.
#'
#' @return Seqinfo object.
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

#' Get genomic features
#'
#' @name features
#' @aliases features,ChromatinContacts,character-method
#' @description
#' Retrieves specific or all genomic features.
#'
#' @param x ChromatinContacts object.
#' @param name Character. Feature name: `"compartments"`, `"TADs"`, `"loops"`,
#'   or `"multi_contacts"`. If missing, returns all as SimpleList.
#'
#' @return GRanges (TADs, compartments, multi_contacts), GInteractions (loops),
#'   or SimpleList (all features).
#'
#' @examples
#' \dontrun{
#' tads <- features(cc, "TADs")
#' loops <- features(cc, "loops")
#' all_features <- features(cc)
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

#' Get interaction data
#'
#' @description
#' Retrieves imported Hi-C interactions.
#'
#' @param x ChromatinContacts object.
#'
#' @return GInteractions (Hi-C data) or NULL (not yet imported).
#'
#' @details
#' Must load with `import()` first. GInteractions contains bin coordinates and
#' metadata (count, balanced).
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

#' Get hypergraph data from MultiWayContacts object
#'
#' @description
#' Retrieves the selected or tidied hypergraph data from a `MultiWayContacts`
#' object.
#'
#' @name hypergraphData
#' @aliases hypergraphData,MultiWayContacts-method
#'
#' @param x A `MultiWayContacts` object.
#' @param which Character. Which hypergraph data to retrieve: `"selected"`
#'   (default) returns the selected subset for visualization, `"tidied"`
#'   returns the full tidied long-format data.
#'
#' @return A data.frame containing hypergraph data in long format, or NULL if
#'   the requested data has not been computed yet.
#'
#' @details
#' - `which = "selected"`: Returns data from `select()`, a subset of hyperedges
#'   chosen for visualization. Use this for plotting.
#' - `which = "tidied"`: Returns data from `tidy()`, the full hypergraph in
#'   long format. Use this for analysis.
#'
#' @examples
#' \dontrun{
#' mc <- MultiWayContacts("sample.pairs.gz") |>
#'   import() |>
#'   build(bin_size = 1e6, min_multiway = 3) |>
#'   tidy(max_hyperedges = 100) |>
#'   select(n = 50)
#'
#' # Get selected hyperedges for plotting
#' df <- hypergraphData(mc, "selected")
#'
#' # Get full tidied data for analysis
#' df_full <- hypergraphData(mc, "tidied")
#' }
#'
#' @export
methods::setMethod(
  "hypergraphData", "MultiWayContacts",
  function(x, which = c("selected", "tidied")) {
    which <- match.arg(which)
    if (which == "selected") {
      x@select_hypergraph
    } else {
      x@tidied_hypergraph
    }
  }
)
