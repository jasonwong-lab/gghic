methods::setClassUnion("GInteractionsOrNULL", c("GInteractions", "NULL"))
methods::setClassUnion(
  "GInteractionsOrCharacterOrNULL", c("GInteractions", "character", "NULL")
)
methods::setClassUnion("CharacterOrNULL", c("character", "NULL"))
methods::setClassUnion("SeqinfoOrNULL", c("Seqinfo", "NULL"))
methods::setClassUnion("IntegerOrNULL", c("integer", "NULL"))
methods::setClassUnion("GRangesOrNULL", c("GRanges", "NULL"))
methods::setClassUnion("DgCMatrixOrNULL", c("dgCMatrix", "NULL"))
methods::setClassUnion("Data.frameOrNULL", c("data.frame", "NULL"))

#' ChromatinContacts S4 class
#'
#' @description
#' S4 class representing chromatin contact data from cooler files (.cool or
#' .mcool). Memory-efficient pointer that loads data only when `import()` is
#' called.
#'
#' @slot cooler_path Character. Path to cooler file.
#' @slot resolution Integer or NULL. Resolution in base pairs (required for
#'   .mcool files).
#' @slot focus GInteractions or NULL. Genomic regions of interest (NULL for
#'   genome-wide).
#' @slot interactions GInteractions or NULL. Imported Hi-C interaction data.
#' @slot seqinfo Seqinfo or NULL. Genome assembly sequence information.
#' @slot compartments GRanges or NULL. A/B compartment annotations.
#' @slot TADs GRanges or NULL. Topologically associating domains.
#' @slot loops GInteractions or NULL. Chromatin loop calls.
#' @slot multi_contacts GRanges or NULL. Multi-way contact regions.
#' @slot tracks GRangesList. Genomic signal tracks (ChIP-seq, ATAC-seq, etc.).
#'
#' @details
#' Key features:
#' * Memory-efficient: loads data only when needed
#' * Flexible focusing: subset regions before or after import
#' * Feature integration: automatically subsets attached features
#' * Multi-resolution: supports both .cool and .mcool files
#'
#' Focus parameter options:
#' * Single region: `"chr1:1000000-2000000"`
#' * Whole chromosome: `"chr1"`
#' * Multiple regions (all combinations): `"chr1 | chr2"`
#' * Inter-chromosomal only: `"chr1 & chr2"`
#' * GInteractions object for programmatic control
#'
#' @seealso
#' [ChromatinContacts()], [import-ChromatinContacts], [features()],
#' [resolution()], [focus()], [gghic()]
#'
#' @examples
#' \dontrun{
#' # Create and import
#' cc <- ChromatinContacts("data.cool") |> import()
#'
#' # Focus on specific region
#' cc <- ChromatinContacts("data.cool", focus = "chr1:1-2e6") |> import()
#'
#' # Multi-resolution cooler
#' cc <- ChromatinContacts("data.mcool", resolution = 10000) |> import()
#'
#' # Add features
#' features(cc, "TADs") <- tads_granges
#' features(cc, "loops") <- loops_ginteractions
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

#' Create a ChromatinContacts object
#'
#' @description
#' Constructor for ChromatinContacts from cooler files (.cool or .mcool).
#'
#' @param cooler_path Character. Path to .cool or .mcool file.
#' @param focus Character, GInteractions, or NULL. Genomic region(s) to focus:
#'   * Single region: `"chr1"` or `"chr1:1-1000000"`
#'   * Inter-chromosomal: `"chr1&chr2"` (excludes diagonal)
#'   * All combinations: `"chr1|chr2"` (includes intra-chromosomal)
#'   * Multiple: `c("chr1&chr2", "chr3|chr4")`
#'   * NULL for genome-wide (default)
#' @param resolution Integer. Resolution in base pairs (required for .mcool).
#'
#' @return A ChromatinContacts object.
#'
#' @examples
#' \dontrun{
#' # Genome-wide .cool file
#' cc <- ChromatinContacts("file.cool")
#'
#' # Multi-resolution .mcool
#' cc <- ChromatinContacts("file.mcool", resolution = 5000L)
#'
#' # Focus on chromosome
#' cc <- ChromatinContacts("file.cool", focus = "chr1")
#'
#' # Inter-chromosomal between chr1 and chr2
#' cc <- ChromatinContacts("file.cool", focus = "chr1&chr2")
#' }
#' @seealso [import-ChromatinContacts], [resolution()], [focus()], [features()]
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


#' MultiWayContacts S4 class
#'
#' @description
#' S4 class for multi-way chromatin contacts from long-read sequencing
#' (Pore-C, Micro-C XL). Represents interactions as hypergraphs where reads
#' connect multiple genomic loci.
#'
#' @slot pairs_path Character. Path to .pairs file.
#' @slot pairs Data.frame or NULL. Loaded pairs data.
#' @slot bin_size Integer or NULL. Genomic bin size in base pairs.
#' @slot focus Character or NULL. Chromosome(s) to analyze.
#' @slot hypergraph dgCMatrix or NULL. Sparse incidence matrix (bins ×
#'   hyperedges).
#' @slot weights Data.frame. Hyperedge weights with normalizations (raw, log,
#'   by_order, minmax).
#' @slot bin_info Data.frame. Genomic bin information.
#' @slot hyperedge_reads SimpleList. Read names per hyperedge.
#' @slot multiways Numeric. Contact count per hyperedge.
#' @slot tidied_hypergraph Data.frame or NULL. Long-format hypergraph.
#' @slot select_hypergraph Data.frame or NULL. Selected hyperedges for plotting.
#'
#' @details
#' Workflow: Create → Import → Build → Tidy → Select → Visualize
#'
#' Key features:
#' * Hypergraph representation: reads as hyperedges connecting bins
#' * Sparse matrix storage for memory efficiency
#' * Multiple weight normalizations
#' * Flexible filtering by order, chromosome, and weight
#'
#' @seealso
#' [MultiWayContacts()], [import-MultiWayContacts], [build()], [tidy()],
#' [select()], [gghypergraph()]
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' mc <- MultiWayContacts("sample.pairs.gz")
#'
#' # Focus on chromosome
#' mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1")
#'
#' # Complete workflow
#' mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1") |>
#'   import() |>
#'   build(bin_size = 1000000L) |>
#'   tidy() |>
#'   select(n_intra = 10, n_inter = 5) |>
#'   gghypergraph()
#' }
#'
#' @export
methods::setClass(
  "MultiWayContacts",
  slots = c(
    pairs_path = "character",
    pairs = "Data.frameOrNULL",
    bin_size = "IntegerOrNULL",
    focus = "CharacterOrNULL",
    hypergraph = "DgCMatrixOrNULL",
    weights = "data.frame",
    bin_info = "data.frame",
    hyperedge_reads = "SimpleList",
    multiways = "numeric",
    tidied_hypergraph = "Data.frameOrNULL",
    select_hypergraph = "Data.frameOrNULL"
  ),
  prototype = list(
    pairs_path = character(0),
    pairs = NULL,
    bin_size = NULL,
    focus = NULL,
    hypergraph = NULL,
    weights = tibble::tibble(),
    bin_info = tibble::tibble(),
    hyperedge_reads = S4Vectors::SimpleList(),
    multiways = numeric(0),
    tidied_hypergraph = NULL,
    select_hypergraph = NULL
  )
)

#' Create a MultiWayContacts object
#'
#' @description
#' Constructor for MultiWayContacts from .pairs files.
#'
#' @param pairs_path Character. Path to .pairs or .pairs.gz file.
#' @param focus Character or NULL. Chromosome(s) to analyze (NULL for
#'   genome-wide).
#' @rdname MultiWayContacts-class
#' @export
MultiWayContacts <- function(pairs_path = NULL, focus = NULL) {
  stopifnot(
    "pairs_path must be a single character string" =
      is.character(pairs_path) & length(pairs_path) == 1
  )
  stopifnot("pairs_path file must exist" = file.exists(pairs_path))
  stopifnot(
    "focus must be a character string or NULL" =
      is.character(focus) | is.null(focus)
  )

  methods::new(
    "MultiWayContacts",
    pairs_path = pairs_path,
    focus = focus
  )
}

methods::setValidity(
  "MultiWayContacts", function(object) {
    if (!is.null(object@focus) && length(object@focus) == 0) {
      return("@focus must not be empty.")
    }
    if (!is.null(object@pairs) && nrow(object@pairs) == 0) {
      return("@pairs must not be empty.")
    }
    if (!is.null(object@hypergraph) && ncol(object@hypergraph) == 0) {
      return("@hypergraph must not be empty.")
    }

    TRUE
  }
)
