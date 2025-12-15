#' @useDynLib gghic, .registration = TRUE
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom InteractionSet GInteractions
#' @importClassesFrom Matrix dgCMatrix
#' @importClassesFrom S4Vectors SimpleList
#' @importClassesFrom Seqinfo Seqinfo
#' @importFrom dplyr bind_rows case_when distinct filter first group_by
#'   left_join mutate reframe rename row_number select slice slice_max
#'   slice_min summarise ungroup
#' @importFrom GenomicRanges end findOverlaps gaps GRanges GRangesList
#'   pintersect reduce seqinfo split start
#' @importFrom ggplot2 aes coord_fixed draw_key_blank draw_key_path
#'   draw_key_polygon element_blank element_line Geom ggplot ggproto
#'   ggproto_parent layer scale_fill_gradientn scale_x_continuous Stat theme
#'   theme_bw %+replace%
#' @importFrom grDevices colorRampPalette
#' @importFrom grid arrow gList gpar is.unit nullGrob pathGrob pointsGrob
#'   polygonGrob polylineGrob segmentsGrob textGrob unit
#' @importFrom InteractionSet anchors GInteractions interactions pairdist
#'   regions replaceRegions<- swapAnchors
#' @importFrom IRanges IRanges overlapsAny subsetByOverlaps
#' @importFrom Matrix colSums nnzero rowSums sparseMatrix summary
#' @importFrom methods as is show
#' @importFrom purrr map map2 map_df pmap
#' @importFrom rhdf5 h5read
#' @importFrom rlang sym := !!
#' @importFrom rtracklayer browserSession genome getTable import ucscTableQuery
#' @importFrom S4Vectors mcols mcols<- SimpleList subjectHits
#' @importFrom scales oob_squish percent unit_format
#' @importFrom Seqinfo seqinfo seqlevels seqnames
#' @importFrom stats setNames
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr drop_na replace_na
#' @importFrom utils globalVariables
## usethis namespace: end
NULL
