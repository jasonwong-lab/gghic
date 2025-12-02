methods::setClassUnion("GInteractionsOrNULL", c("GInteractions", "NULL"))
methods::setClassUnion(
  "GInteractionsOrCharacterOrNULL", c("GInteractions", "character", "NULL")
)
methods::setClassUnion("SeqinfoOrNULL", c("Seqinfo", "NULL"))
methods::setClassUnion("IntegerOrNULL", c("integer", "NULL"))
methods::setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))

#' ChromatinContacts S4 class
#'
#' @description
#' An S4 class to represent chromatin contact data from Hi-C experiments stored
#' in cooler format (.cool or .mcool files). This lightweight object acts as a
#' pointer to Hi-C data without loading the full matrix into memory until
#' explicitly imported.
#'
#' @slot cooler_path Character string. Path to the .cool or .mcool file.
#' @slot resolution Integer or NULL. Resolution in base pairs for .mcool files.
#'   Required for multi-resolution cooler files.
#' @slot focus GInteractions or NULL. Genomic regions of interest to
#'   focus on. Can be specified as a GInteractions object,
#'   or NULL for genome-wide data.
#' @slot interactions GInteractions or NULL. The imported Hi-C interaction data.
#'   Initially NULL until `import()` is called.
#' @slot seqinfo Seqinfo or NULL. Sequence information for the genome assembly.
#' @slot compartments GRanges or NULL. A/B compartment annotations.
#' @slot TADs GRanges or NULL. Topologically associating domains (TADs).
#' @slot loops GInteractions or NULL. Chromatin loop calls.
#' @slot multi_contacts GRanges or NULL. Multi-way contact regions (e.g., from
#'   Pore-C data).
#' @slot tracks GRangesList. Additional genomic signal tracks (e.g., ChIP-seq,
#'   ATAC-seq).
#'
#' @details
#' The `ChromatinContacts` class provides an efficient framework for managing
#' Hi-C data. Key features:
#'
#' * **Memory efficient**: Only loads data when `import()` is called
#' * **Flexible focusing**: Subset regions before or after loading
#' * **Feature integration**: Attach TADs, loops, and tracks that auto-subset
#'   with the interaction data
#' * **Multi-resolution support**: Works with both .cool and .mcool files
#'
#' The `focus` parameter accepts flexible specifications:
#' * Single region: `"chr1:1000000-2000000"`
#' * Whole chromosome: `"chr1"`
#' * Multiple regions with OR: `"chr1 | chr2"` (includes all interactions)
#' * Inter-regional with AND: `"chr1 & chr2"` (inter-chromosomal only)
#' * GInteractions object for programmatic definition
#'
#' @seealso
#' * [ChromatinContacts()] - Constructor function
#' * [import-ChromatinContacts] - Load interaction data
#' * [features()] - Access or set genomic features
#' * [resolution()] - Get resolution
#' * [focus()] - Get focus regions
#' * [gghic()] - Create visualizations
#'
#' @examples
#' \dontrun{
#' # Create object pointing to cooler file
#' cc <- ChromatinContacts(cooler_path = "data.cool")
#'
#' # Create with focus on specific region
#' cc <- ChromatinContacts(
#'   cooler_path = "data.cool",
#'   focus = "chr1:1000000-2000000"
#' )
#'
#' # Multi-resolution cooler with specific resolution
#' cc <- ChromatinContacts(
#'   cooler_path = "data.mcool",
#'   resolution = 10000
#' )
#'
#' # Import the data
#' cc <- import(cc)
#'
#' # Add features
#' features(cc, "TADs") <- tads_granges
#' features(cc, "loops") <- loops_ginteractions
#' features(cc, "tracks") <- tracks_grangeslist
#' }
#'
#' @export
methods::setClass(
  "ChromatinContacts",
  slots = c(
    cooler_path = "character",
    resolution = "IntegerOrNULL",
    focus = "GInteractionsOrCharacterOrNULL",
    interactions = "GInteractionsOrNULL",
    seqinfo = "SeqinfoOrNULL",
    compartments = "GRangesOrNULL",
    TADs = "GRangesOrNULL",
    loops = "GInteractionsOrNULL",
    multi_contacts = "GRangesOrNULL",
    tracks = "GRangesList"
  ),
  prototype = list(
    cooler_path = character(0),
    resolution = NULL,
    focus = NULL,
    interactions = NULL,
    seqinfo = NULL,
    compartments = NULL,
    TADs = NULL,
    loops = NULL,
    multi_contacts = NULL,
    tracks = GenomicRanges::GRangesList()
  )
)

.string2Gr <- function(string, with_coordinates, seq_info) {
  if (with_coordinates) {
    # "chr1:1-2"
    methods::as(string, "GRanges")
  } else {
    # "chr1"
    methods::as(seq_info[string], "GRanges")
  }
}

.checkFocus <- function(focus, seq_info) {
  if (!is.character(focus)) {
    return(focus)
  }

  if (length(focus) == 1) {
    specs <- trimws(strsplit(focus, ",", fixed = TRUE)[[1]])
  } else {
    specs <- trimws(focus)
  }

  gis <- purrr::map(specs, function(spec) {
    has_and <- grepl("&", spec, fixed = TRUE)
    has_or <- grepl("|", spec, fixed = TRUE)

    stopifnot(
      "focus string must use either '&' or '|' operator, not both; use a vector of strings instead" = !(has_and & has_or)
    )

    if (has_and || has_or) {
      operator <- if (has_and) "&" else "|"
      regions <- trimws(strsplit(spec, operator, fixed = TRUE)[[1]])

      grs <- purrr::map(regions, function(region) {
        .string2Gr(region, grepl(":", region), seq_info)
      }) |>
        GenomicRanges::GRangesList() |>
        unlist()
      Seqinfo::seqlevels(grs) <- Seqinfo::seqlevels(seq_info)
      GenomicRanges::seqinfo(grs) <- seq_info

      n_regions <- length(grs)
      pairs <- expand.grid(i = seq_len(n_regions), j = seq_len(n_regions))
      if (operator == "&") {
        pairs <- pairs |>
          dplyr::filter(i < j)
      } else {
        pairs <- pairs |>
          dplyr::filter(i <= j)
      }

      InteractionSet::GInteractions(grs[pairs$i], grs[pairs$j])
    } else {
      gr <- .string2Gr(spec, grepl(":", spec), seq_info)
      InteractionSet::GInteractions(gr, gr)
    }
  })

  focus <- do.call(c, gis)

  focus |>
    InteractionSet::swapAnchors(mode = "order") |>
    sort() |>
    unique()
}

#' Create a ChromatinContacts object
#'
#' @description
#' Constructor function for creating a ChromatinContacts object from a cooler
#' file (.cool or .mcool format).
#'
#' @param cooler_path Character string. Path to a .cool or .mcool file.
#' @param focus Character string, GInteractions object, or NULL. Specifies the
#'   genomic region(s) to focus on. Can be specified as:
#'   * A single region: `"chr1"` or `"chr1:1-1000000"`
#'   * Multiple regions with `&` (inter-chromosomal pairs only):
#'     `"chr1&chr2"` or `"chr1:1-1000000&chr2:1-2000000"`
#'   * Multiple regions with `|` (all combinations including intra):
#'     `"chr1|chr2"` or `"chr1:1-1000000|chr2:1-2000000"`
#'   * Multiple specs as a character vector:
#'       c("chr1&chr2", "chr3:1-2e6|chr4", "chr5")
#'   * (Legacy) Single comma-separated string: "chr1&chr2,chr3|chr4"
#'   * A GInteractions object
#'   * NULL for genome-wide (default)
#' @param resolution Integer. Resolution in base pairs. Required for .mcool
#'   files, ignored for .cool files.
#'
#' @return A ChromatinContacts object.
#'
#' @details
#' The focus parameter allows flexible specification of regions of interest:
#'   * `&` operator: Creates inter-chromosomal pairs only (excludes diagonal)
#'   * `|` operator: Creates all combinations including intra-chromosomal
#'   * Multiple specifications can be combined with commas
#'
#' @examples
#' \dontrun{
#' # Create object from .cool file (genome-wide)
#' cc <- ChromatinContacts("path/to/file.cool")
#'
#' # Create object from .mcool file with specific resolution
#' cc <- ChromatinContacts("path/to/file.mcool", resolution = 5000L)
#'
#' # Focus on a single chromosome
#' cc <- ChromatinContacts("path/to/file.cool", focus = "chr1")
#'
#' # Focus on a specific region
#' cc <- ChromatinContacts("path/to/file.cool", focus = "chr1:1-1000000")
#'
#' # Focus on inter-chromosomal interactions between chr1 and chr2
#' cc <- ChromatinContacts("path/to/file.cool", focus = "chr1&chr2")
#'
#' # Focus on all combinations of chr1 and chr2 (including intra)
#' cc <- ChromatinContacts("path/to/file.cool", focus = "chr1|chr2")
#'
#' # Multiple focus specifications
#' cc <- ChromatinContacts(
#'   "path/to/file.cool",
#'   focus = "chr1&chr2,chr3:1-1000000|chr4"
#' )
#'
#' # Import the interaction data
#' cc <- import(cc)
#' }
#' @seealso [import-ChromatinContacts], [seqinfo()],
#'   [resolution()], [focus()], [features()]
#' @export
ChromatinContacts <- function(cooler_path, focus = NULL, resolution = NULL) {
  stopifnot(
    "cooler_path must be a single character string" =
      is.character(cooler_path) & length(cooler_path) == 1
  )
  stopifnot("cooler_path must exist" = file.exists(cooler_path))
  if (grepl("\\.mcool$", cooler_path)) {
    stopifnot(
      "resolution must be provided for .mcool files" =
        !is.null(resolution)
    )
  } else if (grepl("\\.cool$", cooler_path)) {
    if (!is.null(resolution)) {
      warning("resolution is ignored for .cool files")
      resolution <- NULL
    }
  } else {
    stop("cooler_path must end with .cool or .mcool")
  }
  stopifnot(
    "focus must be a GInteractions object, character string, or NULL" =
      is(focus, "GInteractions") | is.character(focus) | is.null(focus)
  )
  stopifnot(
    "resolution must be a single positive integer" =
      is.integer(resolution) & length(resolution) == 1 & resolution > 0
  )

  info_seq <- .getSeqinfo(cooler_path, resolution)
  focus <- .checkFocus(focus, info_seq)

  methods::new(
    "ChromatinContacts",
    cooler_path = cooler_path,
    resolution = resolution,
    focus = focus,
    interactions = NULL,
    seqinfo = info_seq
  )
}

methods::setValidity("ChromatinContacts", function(object) {
  if (!is.null(object@focus) && length(object@focus) == 0) {
    return("@focus GInteractions must not be empty.")
  }

  TRUE
})
