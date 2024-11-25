name_pkg <- "gghic"

ensure_dir <- function(paths) {
  for (path in paths) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

calculate_add_lengths <- function(data) {
  chroms_add <- data |>
    dplyr::group_by(seqnames1) |>
    dplyr::summarise(end1 = max(end1), start1 = min(start1)) |>
    dplyr::mutate(length = end1 - start1 + 1) |>
    dplyr::ungroup() |>
    dplyr::select(seqnames1, length) |>
    dplyr::distinct() |>
    dplyr::mutate(cum_len = cumsum(c(0, length)[seq_len(dplyr::n())])) |>
    dplyr::select(seqnames1, cum_len) |>
    dplyr::rename(seqname = seqnames1)

  chroms_add
}

calculate_subtract_lengths <- function(data) {
  chroms_sub <- data |>
    dplyr::group_by(seqnames1) |>
    dplyr::slice_min(order_by = start1, n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(seqnames1, start1) |>
    dplyr::distinct() |>
    dplyr::rename(seqname = seqnames1, orignal_start = start1)

  chroms_sub
}

add_columns <- function(data, chroms_add, cols) {
  `:=` <- rlang::`:=`
  `!!` <- rlang::`!!`
  for (col in names(cols)) {
    data <- data |>
      dplyr::mutate(
        !!rlang::sym(cols[col]) := !!rlang::sym(col) + cum_len
      )
  }
  data
}

sub_columns <- function(data, chroms_sub, cols) {
  for (col in names(cols)) {
    data <- data |>
      dplyr::mutate(
        !!rlang::sym(cols[col]) := !!rlang::sym(col) - orignal_start + 1
      )
  }
  data
}

adjust_coordinates <- function(data, columns_list, second = TRUE) {
  chroms_add <- data |>
    calculate_add_lengths()
  chroms_sub <- data |>
    calculate_subtract_lengths()

  data <- data |>
    dplyr::left_join(chroms_add, by = c("seqnames1" = "seqname")) |>
    add_columns(chroms_add, columns_list[[1]]) |>
    dplyr::select(-cum_len)

  if (second) {
    data <- data |>
      dplyr::left_join(chroms_add, by = c("seqnames2" = "seqname")) |>
      add_columns(chroms_add, columns_list[[2]]) |>
      dplyr::select(-cum_len)
  }

  data <- data |>
    dplyr::left_join(chroms_sub, by = c("seqnames1" = "seqname")) |>
    sub_columns(chroms_sub, columns_list[[1]]) |>
    dplyr::select(-orignal_start)

  if (second) {
    data <- data |>
      dplyr::left_join(chroms_sub, by = c("seqnames2" = "seqname")) |>
      sub_columns(chroms_sub, columns_list[[2]]) |>
      dplyr::select(-orignal_start)
  }

  data
}

adjust_coordinates2 <- function(data, chroms_add, chroms_sub, columns) {
  data <- data |>
    dplyr::left_join(chroms_add, by = "seqname") |>
    add_columns(chroms_add, columns) |>
    dplyr::select(-cum_len)

  data <- data |>
    dplyr::left_join(chroms_sub, by = "seqname") |>
    sub_columns(chroms_sub, columns) |>
    dplyr::select(-orignal_start)

  data
}

tbl2gis <- function(data) {
  anchor1 <- GenomicRanges::GRanges(
    seqnames = data$seqnames1,
    ranges = IRanges::IRanges(data$start1, data$end1)
  )
  anchor2 <- GenomicRanges::GRanges(
    seqnames = data$seqnames2,
    ranges = IRanges::IRanges(data$start2, data$end2)
  )

  InteractionSet::GInteractions(anchor1, anchor2)
}

# *--------------------------------------------------------------------------* #
# * Use the cache directory to store the cytoBand, TxDb, txdump, and tx2gene.  #
# *--------------------------------------------------------------------------* #
dir_cache <- rappdirs::user_cache_dir(appname = name_pkg)

stop_if_null <- function(object, message) if (is.null(object)) stop(message)

ensure_cytoband <- function(genome) {
  path_cytoband <- glue::glue("{dir_cache}/cytoBand.{genome}.rds")
  if (!file.exists(path_cytoband)) {
    session <- rtracklayer::browserSession()
    rtracklayer::genome(session) <- genome
    query <- rtracklayer::ucscTableQuery(session, table = "cytoBandIdeo")
    bands_all <- rtracklayer::getTable(query)
    saveRDS(bands_all, path_cytoband)
  } else {
    bands_all <- readRDS(path_cytoband)
  }
  stop_if_null(bands_all, "cytoBand could not be created or loaded.")
  bands_all
}

ensure_txdb <- function(txdb, gtf_path) {
  if (is.null(txdb) && !is.null(gtf_path)) {
    path_txdb <- glue::glue("{dir_cache}/{basename(gtf_path)}.sqlite")
    if (!file.exists(path_txdb)) {
      txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
      AnnotationDbi::saveDb(txdb, path_txdb)
    } else {
      txdb <- suppressMessages(AnnotationDbi::loadDb(path_txdb))
    }
  }
  stop_if_null(txdb, "TxDb could not be created or loaded.")
  txdb
}

ensure_txdump <- function(txdb, gtf_path) {
  path_txdump <- glue::glue("{dir_cache}/{basename(gtf_path)}.txdump.rds")
  if (!file.exists(path_txdump)) {
    txdump <- AnnotationDbi::as.list(txdb)
    saveRDS(txdump, path_txdump)
  } else {
    txdump <- readRDS(path_txdump)
  }
  txdump
}

ensure_tx2gene <- function(tx2gene, gtf_path) {
  if (is.null(tx2gene) && !is.null(gtf_path)) {
    path_tx2gene <- glue::glue("{dir_cache}/{basename(gtf_path)}.tx2gene.rds")
    tsv <- glue::glue("{dir_cache}/{basename(gtf_path)}.tx2gene.tsv")
    if (!file.exists(path_tx2gene)) {
      module_py <- reticulate::import_from_path(
        module = "tx2gene", path = system.file("python", package = name_pkg)
      )
      module_py$generate_tx2gene(gtf_path, tsv)
      vroom::vroom(
        tsv, col_names = c(
          "chrom", "start", "end", "gene_id", "gene_symbol", "tx_id", "tx_name",
          "gene_type", "tx_type"
        ),
        col_types = vroom::cols()
      ) |>
        saveRDS(path_tx2gene)
      file.remove(tsv)
    }
    tx2gene <- readRDS(path_tx2gene)
  }
  stop_if_null(tx2gene, "tx2gene data could not be created or loaded.")
  tx2gene
}

#' download_example_files
#'
#' @description Download example files under the cache directory.
#' @return NULL
#' @examples
#' \dontrun{
#' library(gghic)
#' download_example_files()
#' }
#' @export download_example_files
#' @aliases download_example_files
download_example_files <- function() {
  files <- list(
    "chr4_11-100kb.cool" = "cooler/chr4_11-100kb.cool",
    "chr4_11-5kb.cool" = "cooler/chr4_11-5kb.cool",
    "track1.bigWig" = "bigwig/track1.bigWig",
    "track2.bigWig" = "bigwig/track2.bigWig",
    "gencode-chr4_11.gtf.gz" = "gtf/gencode-chr4_11.gtf.gz",
    "TADs_500kb-chr4_11.tsv" = "tad/TADs_500kb-chr4_11.tsv",
    "loops-chr4_11.txt" = "loop/loops-chr4_11.txt"
  )

  url_base <- paste0(
    "https://raw.githubusercontent.com/",
    "mhjiang97/gghic-data/refs/heads/master/"
  )

  for (file_name in names(files)) {
    file_path <- file.path(dir_cache, file_name)
    if (!file.exists(file_path)) {
      utils::download.file(paste0(url_base, files[[file_name]]), file_path)
    }
  }
}
