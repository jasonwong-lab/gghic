.filterFeatures <- function(cc) {
  for (name in c("compartments", "TADs", "multi_contacts")) {
    feature <- methods::slot(cc, name)
    if (!is.null(feature)) {
      methods::slot(cc, name) <- .confineGRanges(feature, cc)
    }
  }

  for (name in c("loops")) {
    feature <- methods::slot(cc, name)
    if (!is.null(feature)) {
      methods::slot(cc, name) <- .confineGInteractions(feature, cc)
    }
  }

  for (name in c("tracks")) {
    feature <- methods::slot(cc, name)
    if (!is.null(feature) && length(feature) > 0) {
      methods::slot(cc, name) <- purrr::map(
        feature, ~ .confineGRanges(.x, cc)
      ) |>
        GenomicRanges::GRangesList()
    }
  }

  if (!is.null(cc@focus)) {
    cc@focus <- .confineGInteractions(cc@focus, cc)
  } else {
    grs <- cc@interactions |>
      Seqinfo::seqinfo() |>
      GenomicRanges::GRanges()
    n_regions <- length(grs)
    pairs <- expand.grid(i = seq_len(n_regions), j = seq_len(n_regions)) |>
      dplyr::filter(i <= j)
    gis <- InteractionSet::GInteractions(grs[pairs$i], grs[pairs$j])
    cc@focus <- .confineGInteractions(gis, cc)
  }

  cc
}

#' Subset ChromatinContacts objects by various criteria
#'
#' @description
#' Flexibly subset ChromatinContacts objects using multiple indexing methods:
#' numeric indices, logical vectors, character strings (genomic regions),
#' GRanges, or GInteractions. All attached genomic features (TADs, loops,
#' tracks) are automatically filtered to match the subset.
#'
#' @name subset-methods
#' @aliases [,ChromatinContacts,numeric,ANY,ANY-method
#'   [,ChromatinContacts,logical,ANY,ANY-method
#'   [,ChromatinContacts,character,ANY,ANY-method
#'   [,ChromatinContacts,GRanges,ANY,ANY-method
#'   [,ChromatinContacts,GInteractions,ANY,ANY-method
#'   subsetByOverlaps,ChromatinContacts,numeric-method
#'   subsetByOverlaps,ChromatinContacts,logical-method
#'   subsetByOverlaps,ChromatinContacts,GRanges-method
#'   subsetByOverlaps,ChromatinContacts,GInteractions-method
#'
#' @param x A `ChromatinContacts` object.
#' @param i Index, logical vector, character (region string), GRanges, or
#'   GInteractions for subsetting.
#' @param j,drop Ignored. Included for S4 compatibility.
#' @param ranges For `subsetByOverlaps`: GRanges, GInteractions, numeric, or
#'   logical vector to subset by.
#' @param type Character string specifying overlap type for GRanges subsetting.
#'   Can be `"within"` or `"any"`. Default is `"any"`.
#' @param ... Additional arguments passed to subsetting methods.
#'
#' @return A subsetted `ChromatinContacts` object. Features (TADs, loops,
#'   tracks) are automatically subsetted to match the new interaction data.
#'
#' @details
#' Subsetting methods:
#' * **Numeric**: Select interactions by index
#' * **Logical**: Select interactions matching TRUE values
#' * **Character**: Specify genomic region(s) as strings (e.g., `"chr1"` or
#'   `"chr1:1000000-2000000"`)
#' * **GRanges**: Select interactions overlapping the ranges
#' * **GInteractions**: Select interactions overlapping the interaction pairs
#'
#' All methods automatically update:
#' * `interactions` slot
#' * `focus` slot
#' * All attached features (TADs, loops, tracks, etc.)
#'
#' @examples
#' \dontrun{
#' # Load Hi-C data with features
#' cc <- ChromatinContacts("sample.cool") |> import()
#' features(cc, "TADs") <- rtracklayer::import("tads.bed")
#' features(cc, "loops") <- loops_gi
#'
#' # Subset by character region string
#' cc_chr1 <- cc["chr1"]  # Entire chromosome
#' cc_region <- cc["chr1:1000000-5000000"]  # Specific region
#' cc_multi <- cc["chr1|chr2"]  # Multiple regions
#'
#' # Subset by numeric index
#' cc_first_1000 <- cc[1:1000]  # First 1000 interactions
#'
#' # Subset by logical vector
#' gis <- interactions(cc)
#' high_counts <- gis$count > 100
#' cc_high <- cc[high_counts]  # Only high-count interactions
#'
#' # Subset by GRanges
#' region_gr <- GRanges("chr1:1000000-2000000")
#' cc_overlap <- cc[region_gr]  # Interactions overlapping region
#'
#' # Subset by GInteractions
#' specific_pairs <- GInteractions(
#'   GRanges("chr1:1000000-1100000"),
#'   GRanges("chr1:2000000-2100000")
#' )
#' cc_pairs <- cc[specific_pairs]
#'
#' # All features are automatically subset
#' length(features(cc, "TADs"))  # Original TAD count
#' length(features(cc_chr1, "TADs"))  # Fewer TADs after subsetting
#'
#' # By GRanges
#' roi <- GRanges("chr1:1000000-5000000")
#' cc_sub <- cc[roi]
#'
#' # By GInteractions
#' query_gi <- GInteractions(
#'   GRanges("chr1:1000000-2000000"),
#'   GRanges("chr1:3000000-4000000")
#' )
#' cc_sub <- cc[query_gi]
#' }
#'
#' @export
methods::setMethod(
  "subsetByOverlaps",
  signature = c("ChromatinContacts", "numeric"),
  function(x, ranges) {
    interactions_sub <- x@interactions[ranges]
    regions_sub <- subsetByOverlaps(
      InteractionSet::regions(x@interactions), interactions_sub
    )
    InteractionSet::replaceRegions(interactions_sub) <- regions_sub
    x@interactions <- interactions_sub |>
      InteractionSet::reduceRegions()

    x <- .filterFeatures(x)

    x
  }
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "subsetByOverlaps",
  signature = c("ChromatinContacts", "logical"),
  function(x, ranges) {
    subsetByOverlaps(x, which(ranges))
  }
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "subsetByOverlaps",
  signature = c("ChromatinContacts", "GRanges"),
  function(x, ranges, type = c("within", "any")) {
    type <- match.arg(type)
    if (!type %in% c("within", "any")) stop("Unsupported type")

    Seqinfo::seqlevels(ranges) <- Seqinfo::seqlevels(x@seqinfo)
    Seqinfo::seqinfo(ranges) <- x@seqinfo
    ranges <- ranges |>
      sort() |>
      unique()

    if (type == "any") {
      overlaps <- IRanges::overlapsAny(x@interactions, ranges, type = "within")
      x@interactions <- x@interactions[overlaps] |>
        InteractionSet::reduceRegions()
    }

    if (type == "within") {
      x@interactions <- subsetByOverlaps(
        x@interactions, InteractionSet::GInteractions(ranges, ranges)
      ) |>
        InteractionSet::reduceRegions()
    }

    x <- .filterFeatures(x)

    x
  }
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "subsetByOverlaps",
  signature = c("ChromatinContacts", "GInteractions"),
  function(x, ranges) {
    Seqinfo::seqlevels(ranges) <- Seqinfo::seqlevels(x@seqinfo)
    Seqinfo::seqinfo(ranges) <- x@seqinfo
    ranges <- InteractionSet::swapAnchors(ranges, mode = "order") |>
      sort() |>
      unique()
    anchors <- x@interactions |>
      InteractionSet::anchors()
    idx <- purrr::map(seq_along(ranges), function(i) {
      which(
        IRanges::overlapsAny(
          anchors[[1]], InteractionSet::anchors(ranges[i], "first"),
          type = "within"
        ) &
          IRanges::overlapsAny(
            anchors[[2]], InteractionSet::anchors(ranges[i], "second"),
            type = "within"
          )
      )
    }) |>
      purrr::flatten_int() |>
      unique()

    x@interactions <- x@interactions[idx] |>
      InteractionSet::reduceRegions()
    x <- .filterFeatures(x)

    x
  }
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "[", signature("ChromatinContacts", "numeric"),
  function(x, i) subsetByOverlaps(x, i)
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "[", signature("ChromatinContacts", "logical"),
  function(x, i) subsetByOverlaps(x, i)
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "[", signature("ChromatinContacts", "GRanges"),
  function(x, i, ..., drop = TRUE) subsetByOverlaps(x, i, ...)
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "[", signature("ChromatinContacts", "GInteractions"),
  function(x, i) subsetByOverlaps(x, i)
)

#' @rdname subset-methods
#' @export
methods::setMethod(
  "[", signature("ChromatinContacts", "character"),
  function(x, i) {
    if (length(i) == 1) {
      specs <- trimws(strsplit(i, ",", fixed = TRUE)[[1]])
    } else {
      specs <- trimws(i)
    }

    interactions <- x@interactions

    gis <- purrr::map(specs, function(spec) {
      has_and <- grepl("&", spec, fixed = TRUE)
      has_or <- grepl("|", spec, fixed = TRUE)
      operator <- if (has_and) "&" else "|"

      regions <- trimws(strsplit(spec, operator, fixed = TRUE)[[1]])
      grs <- purrr::map(regions, function(region) {
        .string2Gr(region, grepl(":", region), x@seqinfo)
      }) |>
        GenomicRanges::GRangesList() |>
        unlist()
      Seqinfo::seqlevels(grs) <- Seqinfo::seqlevels(x@seqinfo)
      Seqinfo::seqinfo(grs) <- x@seqinfo

      if (has_or) {
        n_regions <- length(grs)
        pairs <- expand.grid(i = seq_len(n_regions), j = seq_len(n_regions)) |>
          dplyr::filter(i < j)
        gi <- InteractionSet::GInteractions(grs[pairs$i], grs[pairs$j])
        subsetByOverlaps(interactions, gi)
      } else if (has_and) {
        subsetByOverlaps(interactions, grs, "within")
      } else {
        gi <- InteractionSet::GInteractions(grs, grs)
        subsetByOverlaps(interactions, gi)
      }
    })

    gis <- do.call(c, gis) |>
      InteractionSet::swapAnchors(mode = "order") |>
      sort() |>
      unique()

    x@interactions <- gis |>
      InteractionSet::reduceRegions()
    x <- .filterFeatures(x)

    x
  }
)
