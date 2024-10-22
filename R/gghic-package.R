#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom AnnotationDbi as.list loadDb saveDb
#' @importFrom biovizBase getBioColor
#' @importFrom dplyr bind_rows case_when distinct filter first group_by left_join mutate row_number select slice slice_max slice_min summarise ungroup
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicRanges findOverlaps gaps GRanges GRangesList pintersect reduce split
#' @importFrom ggplot2 aes coord_fixed draw_key_blank draw_key_polygon element_blank element_line Geom ggplot ggproto layer scale_fill_gradientn scale_x_continuous Stat theme theme_bw %+replace%
#' @importFrom ggrastr rasterise
#' @importFrom glue glue
#' @importFrom grid arrow gList gpar nullGrob polygonGrob polylineGrob segmentsGrob textGrob unit
#' @importFrom Gviz GeneRegionTrack
#' @importFrom InteractionSet interactions pairdist
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom purrr map map2 map_df pmap
#' @importFrom rappdirs user_cache_dir
#' @importFrom reticulate import_from_path
#' @importFrom rlang sym := !!
#' @importFrom rtracklayer browserSession genome getTable import ucscTableQuery
#' @importFrom S4Vectors mcols subjectHits
#' @importFrom scales oob_squish unit_format
#' @importFrom tibble as_tibble is_tibble tibble
#' @importFrom tidyr drop_na replace_na
#' @importFrom txdbmaker makeTxDb makeTxDbFromGFF
#' @importFrom vroom cols vroom
## usethis namespace: end
NULL
