#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom AnnotationDbi as.list loadDb saveDb
#' @importFrom dplyr bind_rows case_when distinct filter first group_by
#'   left_join mutate reframe rename row_number select slice slice_max
#'   slice_min summarise ungroup
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges end findOverlaps gaps GRanges GRangesList
#'   pintersect reduce split start
#' @importFrom ggplot2 aes coord_fixed draw_key_blank draw_key_path
#'   draw_key_polygon element_blank element_line Geom ggplot ggproto layer
#'   scale_fill_gradientn scale_x_continuous Stat theme theme_bw %+replace%
#' @importFrom glue glue
#' @importFrom grDevices colorRampPalette
#' @importFrom grid arrow gList gpar is.unit nullGrob pathGrob pointsGrob
#'   polygonGrob polylineGrob segmentsGrob textGrob unit
#' @importFrom Gviz GeneRegionTrack
#' @importFrom InteractionSet anchors interactions pairdist
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom methods is
#' @importFrom purrr map map2 map_df pmap
#' @importFrom rappdirs user_cache_dir
#' @importFrom reticulate import_from_path
#' @importFrom rlang sym := !!
#' @importFrom rtracklayer browserSession genome getTable import ucscTableQuery
#' @importFrom S4Vectors mcols mcols<- subjectHits
#' @importFrom scales oob_squish unit_format
#' @importFrom stats setNames
#' @importFrom tibble as_tibble is_tibble tibble
#' @importFrom tidyr drop_na replace_na
#' @importFrom txdbmaker makeTxDb makeTxDbFromGFF
#' @importFrom utils globalVariables
#' @importFrom vroom cols vroom
## usethis namespace: end
NULL
