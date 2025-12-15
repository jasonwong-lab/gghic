.getPkgName <- function() "gghic"

.ensureDir <- function(paths) {
  for (path in paths) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

.calculateAddLengths <- function(data) {
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

.calculateAddLengths1d <- function(data) {
  chroms_add <- data |>
    dplyr::group_by(seqnames) |>
    dplyr::summarise(end = max(end), start = min(start)) |>
    dplyr::mutate(length = end - start + 1) |>
    dplyr::ungroup() |>
    dplyr::select(seqnames, length) |>
    dplyr::distinct() |>
    dplyr::mutate(cum_len = cumsum(c(0, length)[seq_len(dplyr::n())])) |>
    dplyr::select(seqnames, cum_len) |>
    dplyr::rename(seqname = seqnames)

  chroms_add
}

.calculateSubtractLengths <- function(data) {
  chroms_sub <- data |>
    dplyr::group_by(seqnames1) |>
    dplyr::slice_min(order_by = start1, n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(seqnames1, start1) |>
    dplyr::distinct() |>
    dplyr::rename(seqname = seqnames1, orignal_start = start1)

  chroms_sub
}

.calculateSubtractLengths1d <- function(data) {
  chroms_sub <- data |>
    dplyr::group_by(seqnames) |>
    dplyr::slice_min(order_by = start, n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(seqnames, start) |>
    dplyr::distinct() |>
    dplyr::rename(seqname = seqnames, orignal_start = start)

  chroms_sub
}

.addColumns <- function(data, chroms_add, cols) {
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

.subColumns <- function(data, chroms_sub, cols) {
  for (col in names(cols)) {
    data <- data |>
      dplyr::mutate(
        !!rlang::sym(cols[col]) := !!rlang::sym(col) - orignal_start + 1
      )
  }
  data
}

.adjustCoordinates <- function(data, columns_list, second = TRUE) {
  chroms_add <- data |>
    .calculateAddLengths()
  chroms_sub <- data |>
    .calculateSubtractLengths()

  data <- data |>
    dplyr::left_join(chroms_add, by = c("seqnames1" = "seqname")) |>
    .addColumns(chroms_add, columns_list[[1]]) |>
    dplyr::select(-cum_len)

  if (second) {
    data <- data |>
      dplyr::left_join(chroms_add, by = c("seqnames2" = "seqname")) |>
      .addColumns(chroms_add, columns_list[[2]]) |>
      dplyr::select(-cum_len)
  }

  data <- data |>
    dplyr::left_join(chroms_sub, by = c("seqnames1" = "seqname")) |>
    .subColumns(chroms_sub, columns_list[[1]]) |>
    dplyr::select(-orignal_start)

  if (second) {
    data <- data |>
      dplyr::left_join(chroms_sub, by = c("seqnames2" = "seqname")) |>
      .subColumns(chroms_sub, columns_list[[2]]) |>
      dplyr::select(-orignal_start)
  }

  data
}

.adjustCoordinates2 <- function(data, chroms_add, chroms_sub, columns) {
  data <- data |>
    dplyr::left_join(chroms_add, by = "seqname") |>
    .addColumns(chroms_add, columns) |>
    dplyr::select(-cum_len)

  data <- data |>
    dplyr::left_join(chroms_sub, by = "seqname") |>
    .subColumns(chroms_sub, columns) |>
    dplyr::select(-orignal_start)

  data
}

.tbl2Gis <- function(data) {
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

.grs2Gis <- function(grs, res) {
  starts_1 <- start(grs)
  ends_2 <- end(grs)
  starts_2 <- starts_1 + res - 1
  ends_1 <- ends_2 - res + 1

  anchor1 <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(grs),
    ranges = IRanges::IRanges(starts_1, ends_1)
  )
  anchor2 <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(grs),
    ranges = IRanges::IRanges(starts_2, ends_2)
  )

  InteractionSet::GInteractions(anchor1, anchor2)
}

#' Generate axis breaks and labels for multi-chromosome plots
#'
#' @description
#' Calculates appropriate breaks and labels for x-axis when visualizing
#' multiple chromosomes simultaneously. Adjusts break density based on
#' chromosome length.
#'
#' @param data A data frame or tibble with columns `seqnames`, `start`, and
#'   `end` representing genomic coordinates.
#'
#' @return A list with two elements:
#'   * `breaks`: Numeric vector of x-axis break positions
#'   * `labels`: Character vector of formatted labels
#'
#' @details
#' The function:
#' 1. Groups data by chromosome
#' 2. Calculates chromosome lengths
#' 3. Assigns more breaks to longer chromosomes (2-10 breaks per chromosome)
#' 4. Formats labels in megabases (M)
#' 5. Adds newlines before chromosome transitions for clarity
#'
#' Useful for multi-chromosome Hi-C visualizations where chromosomes are
#' displayed side-by-side.
#'
#' @examples
#' \dontrun{
#' # Get breaks and labels for multi-chromosome plot
#' df <- scaleData(cc, "balanced", log10)
#' breaks_labels <- getBreaksLabels(df)
#'
#' ggplot(df) +
#'   geom_hic() +
#'   scale_x_continuous(
#'     breaks = breaks_labels$breaks,
#'     labels = breaks_labels$labels
#'   )
#' }
#'
#' @export
#' @aliases get_breaks_labels
getBreaksLabels <- function(data) {
  labels <- scales::unit_format(unit = "M", scale = 1e-6, accuracy = 0.01)

  lengths <- data |>
    dplyr::group_by(seqnames) |>
    dplyr::summarise(length = max(end) - min(start)) |>
    dplyr::mutate(
      n_breaks = pmax(2, round(2 + 8 * (length / max(length))))
    )

  tmp <- data |>
    dplyr::left_join(lengths, by = "seqnames") |>
    dplyr::group_by(seqnames) |>
    dplyr::reframe(
      bin = seq(
        from = min(end), to = max(end), length.out = dplyr::first(n_breaks)
      )
    ) |>
    dplyr::rename(seqname = seqnames)

  chroms_add <- data |>
    .calculateAddLengths1d()
  chroms_sub <- data |>
    .calculateSubtractLengths1d()

  tmp2 <- tmp |>
    .adjustCoordinates2(chroms_add, chroms_sub, c(bin = "bin"))

  breaks <- tmp2$bin
  labels <- labels(tmp$bin)

  indices <- match(unique(tmp$seqname), tmp$seqname)[-1]
  labels[indices] <- paste0("\n", labels[indices])

  list(breaks = breaks, labels = labels)
}

.getMemLim <- function() {
  tryCatch(
    {
      if (Sys.info()["sysname"] == "Linux") {
        # Linux: read from /proc/meminfo
        meminfo <- readLines("/proc/meminfo")
        mem_total_line <- grep("^MemTotal:", meminfo, value = TRUE)
        as.numeric(sub("MemTotal:\\s+(\\d+).*", "\\1", mem_total_line)) / 1024
      } else if (Sys.info()["sysname"] == "Darwin") {
        # macOS: use sysctl
        mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
        mem_bytes / (1024^2)
      } else if (Sys.info()["sysname"] == "Windows") {
        # Windows: use memory.size()
        memory.limit()
      } else {
        # Unknown OS - use conservative default
        warning("Unknown OS, using default 16 GB", call. = FALSE)
        16000
      }
    },
    error = function(e) {
      warning(
        "Could not detect system memory, using default 16 GB",
        call. = FALSE
      )
      16000
    }
  )
}

# *--------------------------------------------------------------------------* #
# * geom_annotation                                                          * #
# *--------------------------------------------------------------------------* #
.extractTrs <- function(gr, gene, txdump, include_ncrna, gene_info) {
  if (!requireNamespace("txdbmaker", quietly = TRUE) ||
    !requireNamespace("Gviz", quietly = TRUE)) {
    stop(
      "Packages 'txdbmaker' and 'Gviz' are required for gene annotation. ",
      "Install them with: BiocManager::install(c('txdbmaker', 'Gviz'))",
      call. = FALSE
    )
  }

  `!!` <- rlang::`!!`

  chrom <- as.character(gr@seqnames)
  start <- gr@ranges@start
  end <- start + gr@ranges@width - 1

  r <- range(
    gene[
      S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(gr, gene, ignore.strand = TRUE)
      )
    ],
    ignore.strand = TRUE
  )

  .txdump <- txdump
  .txdump$chrominfo <- .txdump$chrominfo[
    .txdump$chrominfo$chrom == chrom, ,
    drop = FALSE
  ]

  if (is.null(gene_info)) {
    .txdump$transcripts <- .txdump$transcripts[
      .txdump$transcripts$tx_chrom == chrom &
        (
          .txdump$transcripts$tx_start < GenomicRanges::end(r) |
            .txdump$transcripts$tx_end > GenomicRanges::start(r)
        ), ,
      drop = FALSE
    ]

    transcripts <- .txdump$transcripts
  }
  if (!is.null(gene_info)) {
    gene_info <- gene_info |>
      dplyr::filter(
        chrom == !!chrom, start < GenomicRanges::end(r) |
          end > GenomicRanges::start(r)
      )
    if (nrow(gene_info) == 0) {
      return(GenomicRanges::GRangesList())
    }

    transcripts <- .txdump$transcripts[
      .txdump$transcripts$tx_name %in% gene_info$tx_id, ,
      drop = FALSE
    ]

    .txdump$transcripts <- .txdump$transcripts[
      .txdump$transcripts$tx_name %in% gene_info$tx_id, ,
      drop = FALSE
    ]
  }

  transcripts_with_genes <- transcripts |>
    dplyr::left_join(.txdump$genes, by = "tx_id") |>
    dplyr::filter(!is.na(gene_id))

  longest_transcripts <- transcripts_with_genes |>
    dplyr::mutate(length = tx_end - tx_start) |>
    dplyr::group_by(gene_id) |>
    dplyr::filter(length == max(length)) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  .txdump$transcripts <- .txdump$transcripts[
    .txdump$transcripts$tx_id %in% longest_transcripts$tx_id, ,
    drop = FALSE
  ]

  .txdump$genes <- .txdump$genes[
    .txdump$genes$tx_id %in% .txdump$transcripts$tx_id, ,
    drop = FALSE
  ]
  .txdump$splicings <- .txdump$splicings[
    .txdump$splicings$tx_id %in% .txdump$transcripts$tx_id, ,
    drop = FALSE
  ]

  txdb_subset <- do.call(txdbmaker::makeTxDb, .txdump)

  exons_gviz <- Gviz::GeneRegionTrack(
    txdb_subset,
    chromosome = chrom, start = GenomicRanges::start(r),
    end = GenomicRanges::end(r),
    strand = "*"
  )
  exons <- exons_gviz@range

  if (!include_ncrna) {
    ncrnas <- unique(exons$transcript[exons$feature == "ncRNA"])
    exons <- exons[!exons$transcript %in% ncrnas]
  }

  trs <- GenomicRanges::split(exons, as.character(exons$transcript))

  trs_exon_intron <- purrr::map(trs, function(tr) {
    introns <- GenomicRanges::gaps(tr)[-1]
    if (length(introns) > 0) {
      introns$gene <- unique(tr$gene)
      introns$feature <- "intron"
    }
    tmp <- GenomicRanges::pintersect(c(tr, introns), gr, ignore.strand = TRUE)
    tmp[width(tmp) > 0]
  })

  GenomicRanges::GRangesList(trs_exon_intron)
}

.decideTrackLines <- function(genes_span, maxgap) {
  names_gene <- names_gene_copy <- names(genes_span)
  dat_line <- tibble::tibble(gene_id = names_gene, line = 0)

  while (length(names_gene_copy) > 1) {
    n1 <- names_gene_copy[1]
    n1o <- setdiff(names_gene_copy, n1)
    genes_n1 <- genes_span[[n1]]
    genes_n1o <- GenomicRanges::GRangesList(genes_span[n1o]) |>
      unlist()
    genes_n1o_nonol <- genes_n1o[
      !IRanges::overlapsAny(
        genes_n1o, genes_n1,
        maxgap = maxgap, ignore.strand = TRUE
      )
    ]
    n1o_nonol <- names(genes_n1o_nonol)
    if (length(n1o_nonol) == 0) {
      dat_line$line[dat_line$gene_id == n1] <- max(dat_line$line) + 1
      names_gene_copy <- setdiff(names_gene_copy, n1)
      next
    }
    dat_line$line[dat_line$gene_id %in% c(n1, n1o_nonol)] <-
      max(dat_line$line) + 1
    names_gene_copy <- setdiff(names_gene_copy, n1o_nonol)
  }
  if (sum(dat_line$line == 0)) {
    dat_line$line[dat_line$line == 0] <- max(dat_line$line) + 1
  }

  dat_line
}

.retriveGenes <- function(
  data, txdb, tx2gene, gtf_path, maxgap, include_ncrna, gene_symbols,
  chrom_prefix
) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "Package 'GenomicFeatures' is required for gene annotation. ",
      "Install it with: BiocManager::install('GenomicFeatures')",
      call. = FALSE
    )
  }

  txdb <- .ensureTxdb(txdb, gtf_path)
  tx2gene <- .ensureTx2gene(tx2gene, gtf_path)
  txdump <- .ensureTxdump(txdb, gtf_path)

  info_genes <- NULL
  if (!is.null(gene_symbols)) {
    info_genes <- tx2gene |>
      dplyr::filter(gene_symbol %in% gene_symbols)
  }

  grs <- purrr::map(unique(data$seqnames1), function(.x) {
    tmp <- .x
    if (!chrom_prefix) tmp <- paste0("chr", tmp)
    GenomicRanges::GRanges(
      seqnames = tmp,
      ranges = IRanges::IRanges(
        start = min(data$start1[data$seqnames1 == .x]),
        end = max(data$end1[data$seqnames1 == .x])
      )
    )
  })

  genes <- purrr::map(grs, function(.x) {
    GenomicFeatures::genes(
      txdb,
      columns = c("exon_id"), filter = list(tx_chrom = .x@seqnames)
    )
  })

  trs <- purrr::map2(
    grs, genes, .extractTrs,
    txdump = txdump, include_ncrna = include_ncrna, gene_info = info_genes
  )
  trs <- purrr::keep(trs, ~ length(.x) > 0)

  dat <- purrr::pmap(
    list(trs), function(.x) {
      genes <- GenomicRanges::split(unlist(.x), unlist(.x)$gene)

      genes_reduced <- purrr::map(
        genes, function(g) {
          features <- GenomicRanges::split(g, g$feature)
          features_reduced <- purrr::map(
            features, function(f) {
              if (f$feature[1] == "intron") {
                tmp <- f
                S4Vectors::mcols(tmp) <- NULL
                names(tmp) <- NULL
              } else {
                tmp <- GenomicRanges::reduce(f)
              }
              tmp$feature <- unique(f$feature)
              tmp$gene_id <- unique(g$gene)
              tmp
            }
          )
        }
      )

      genes_span <- purrr::map(
        genes_reduced, function(.x) {
          tmp <- unlist(GenomicRanges::GRangesList(.x))
          GenomicRanges::GRanges(
            seqnames = tmp@seqnames[1],
            ranges = IRanges::IRanges(
              start = min(GenomicRanges::start(tmp)),
              end = max(GenomicRanges::end(tmp))
            ),
            strand = tmp@strand[1]
          )
        }
      )

      dat_line <- .decideTrackLines(genes_span, maxgap)

      dat_gene <- purrr::map_df(
        genes_reduced, function(.x) {
          tmp <- GenomicRanges::GRangesList(.x) |>
            unlist()
          tibble::as_tibble(stats::setNames(tmp, NULL))
        }
      ) |>
        dplyr::left_join(dat_line, by = "gene_id") |>
        dplyr::left_join(
          dplyr::select(
            dplyr::distinct(tx2gene, gene_id, .keep_all = TRUE),
            gene_id, gene_symbol
          ),
          by = "gene_id"
        )

      dat_gene
    }
  )

  dplyr::bind_rows(dat)
}

# *--------------------------------------------------------------------------* #
# * geom_hic, geom_hic_under                                                 * #
# *--------------------------------------------------------------------------* #
.calculateHicCoordinates <- function(data, lower = FALSE) {
  n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))

  if (n_sn == 1 || (n_sn == 2 && all(data$seqnames1 != data$seqnames2))) {
    dat <- data
  } else {
    dat <- data |>
      .adjustCoordinates(
        list(
          c(start1 = "start1", end1 = "end1"),
          c(start2 = "start2", end2 = "end2")
        )
      )
  }

  dat <- dat |>
    dplyr::mutate(
      x = (end1 + start2) / 2,
      xmin = (start1 + start2) / 2,
      xmax = (start1 + end2) / 2,
      xend = (end1 + end2) / 2,
      y = (x - end1),
      ymin = (xmin - start1),
      ymax = (xmax - start1),
      yend = (xend - end1)
    )

  if (lower) {
    dat <- dat |>
      dplyr::mutate(
        y = -1 * y, ymin = -1 * ymin, ymax = -1 * ymax, yend = -1 * yend
      )
  }

  dat
}

# *--------------------------------------------------------------------------* #
# * geom_ideogram                                                            * #
# *--------------------------------------------------------------------------* #
.retriveCytoband <- function(data, genome, chrom_prefix) {
  bands_all <- .ensureCytoband(genome = genome)

  if (!chrom_prefix) {
    bands_all <- bands_all |>
      dplyr::mutate(chrom = gsub("chr", "", chrom))
  }

  chroms <- unique(c(data$seqnames1, data$seqnames2)) |>
    as.character()
  bands <- bands_all |>
    dplyr::filter(chrom %in% chroms)

  bnames <- bands$name
  indices <- is.na(bnames)
  if (any(indices)) {
    bnames[indices] <- paste("band_na", seq_len(sum(indices)), sep = "_")
  }
  if (any(bnames == "")) {
    bnames[bnames == ""] <- paste("band_null", which(bnames == ""), sep = "_")
  }

  cols_all <- .getCytobandColors()
  cols <- c(
    cols_all[c("gneg", "stalk", "acen")],
    gpos = unname(cols_all["gpos100"]),
    gvar = unname(cols_all["gpos100"])
  )
  gpcols <- unique(grep("gpos", bands$gieStain, value = TRUE))
  crmp <- grDevices::colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
  posCols <- stats::setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
  cols <- c(cols, posCols) |>
    as.data.frame() |>
    tibble::rownames_to_column("type") |>
    dplyr::rename(band_fill = `c(cols, posCols)`)

  dat <- bands |>
    dplyr::mutate(chromStart = chromStart + 1, name = bnames) |>
    dplyr::rename(
      type = gieStain, seqname = chrom, start = chromStart, end = chromEnd
    ) |>
    dplyr::left_join(cols, by = "type")

  dat
}

# *--------------------------------------------------------------------------* #
# * GRanges, GInteractions                                                   * #
# *--------------------------------------------------------------------------* #
.confineGRanges <- function(grs, cc) {
  grs <- grs[IRanges::overlapsAny(
    grs, InteractionSet::regions(cc@interactions),
    ignore.strand = TRUE
  )]
  Seqinfo::seqlevels(grs) <- Seqinfo::seqlevels(cc@seqinfo)
  Seqinfo::seqinfo(grs) <- cc@seqinfo

  grs <- grs |>
    sort() |>
    unique()

  grs
}

.confineGInteractions <- function(gis, cc) {
  gis <- gis[IRanges::overlapsAny(
    gis, cc@interactions, # InteractionSet::regions(cc@interactions),
    ignore.strand = TRUE,
    use.region = "both"
  )]
  Seqinfo::seqlevels(gis) <- Seqinfo::seqlevels(cc@seqinfo)
  Seqinfo::seqinfo(gis) <- cc@seqinfo

  gis <- gis |>
    InteractionSet::reduceRegions() |>
    InteractionSet::swapAnchors(mode = "order") |>
    sort() |>
    unique()

  gis
}

# *--------------------------------------------------------------------------* #
# * Use the cache directory to store the cytoBand, TxDb, txdump, and tx2gene.  #
# *--------------------------------------------------------------------------* #
.getCacheDir <- function() {
  # Determine cache directory based on OS
  if (.Platform$OS.type == "windows") {
    # Windows: use LOCALAPPDATA or USERPROFILE
    cache_home <- Sys.getenv("LOCALAPPDATA", "")
    if (cache_home == "") {
      cache_home <- Sys.getenv("USERPROFILE", "")
      cache_home <- file.path(cache_home, "AppData", "Local")
    }
  } else if (Sys.info()["sysname"] == "Darwin") {
    # macOS: use ~/Library/Caches
    cache_home <- "~/Library/Caches"
  } else {
    # Linux/Unix: use XDG_CACHE_HOME or ~/.cache
    cache_home <- Sys.getenv("XDG_CACHE_HOME", "")
    if (cache_home == "") {
      cache_home <- "~/.cache"
    }
  }

  cache_dir <- file.path(cache_home, .getPkgName())
  cache_dir <- path.expand(cache_dir)

  # Create directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cache_dir
}

.stopIfNull <- function(object, message) if (is.null(object)) stop(message)

.ensureCytoband <- function(genome) {
  dir_cache <- .getCacheDir()
  path_cytoband <- file.path(dir_cache, paste0("cytoBand.", genome, ".rds"))
  if (!file.exists(path_cytoband)) {
    session <- rtracklayer::browserSession()
    rtracklayer::genome(session) <- genome
    query <- rtracklayer::ucscTableQuery(session, table = "cytoBandIdeo")
    bands_all <- rtracklayer::getTable(query)
    saveRDS(bands_all, path_cytoband)
  } else {
    bands_all <- readRDS(path_cytoband)
  }
  .stopIfNull(bands_all, "cytoBand could not be created or loaded.")
  bands_all
}

.ensureTxdb <- function(txdb, gtf_path) {
  if (!requireNamespace("txdbmaker", quietly = TRUE) ||
    !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Packages 'txdbmaker' and 'AnnotationDbi' are required for gene annotation. ",
      "Install them with: BiocManager::install(c('txdbmaker', 'AnnotationDbi'))",
      call. = FALSE
    )
  }

  dir_cache <- .getCacheDir()
  if (is.null(txdb) && !is.null(gtf_path)) {
    path_txdb <- file.path(dir_cache, paste0(basename(gtf_path), ".sqlite"))
    if (!file.exists(path_txdb)) {
      txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
      AnnotationDbi::saveDb(txdb, path_txdb)
    } else {
      txdb <- suppressMessages(AnnotationDbi::loadDb(path_txdb))
    }
  }
  .stopIfNull(txdb, "TxDb could not be created or loaded.")
  txdb
}

.ensureTxdump <- function(txdb, gtf_path) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Package 'AnnotationDbi' is required for gene annotation. ",
      "Install it with: BiocManager::install('AnnotationDbi')",
      call. = FALSE
    )
  }

  dir_cache <- .getCacheDir()
  path_txdump <- file.path(dir_cache, paste0(basename(gtf_path), ".txdump.rds"))
  if (!file.exists(path_txdump)) {
    txdump <- AnnotationDbi::as.list(txdb)
    saveRDS(txdump, path_txdump)
  } else {
    txdump <- readRDS(path_txdump)
  }
  txdump
}

.generateTx2gene <- function(gtf_path) {
  grs <- rtracklayer::import(gtf_path)
  grs <- grs[grs$type == "transcript"]

  tibble::tibble(
    chrom = as.character(GenomicRanges::seqnames(grs)),
    gene_id = grs$gene_id,
    gene_symbol = grs$gene_name,
    tx_id = grs$transcript_id,
    gene_type = grs$gene_type
  )
}

.ensureTx2gene <- function(tx2gene, gtf_path) {
  dir_cache <- .getCacheDir()
  name_pkg <- .getPkgName()
  if (is.null(tx2gene) && !is.null(gtf_path)) {
    path_tx2gene <- file.path(
      dir_cache, paste0(basename(gtf_path), ".tx2gene.rds")
    )
    tsv <- file.path(dir_cache, paste0(basename(gtf_path), ".tx2gene.tsv"))
    if (!file.exists(path_tx2gene)) {
      .generateTx2gene(gtf_path) |>
        saveRDS(path_tx2gene)
    }
    tx2gene <- readRDS(path_tx2gene)
  }
  .stopIfNull(tx2gene, "tx2gene data could not be created or loaded.")
  tx2gene
}

.getBreaksLabels <- function(data) {
  n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
  fmt_label <- scales::unit_format(unit = "M", scale = 1e-6)

  if (n_sn == 1 || (n_sn == 2 && all(data$seqnames1 != data$seqnames2))) {
    return(list(breaks = ggplot2::waiver(), labels = fmt_label))
  }

  # Calculate lengths and assign breaks scaled to chromosome length
  lengths <- data |>
    dplyr::group_by(seqnames1) |>
    dplyr::summarise(length = max(end1) - min(end1)) |>
    dplyr::mutate(
      n_breaks = pmax(2, pmin(5, round(2 + 3 * (length / max(length)))))
    )

  tmp <- data |>
    dplyr::left_join(lengths, by = "seqnames1") |>
    dplyr::group_by(seqnames1) |>
    dplyr::reframe(
      bin = seq(from = min(end1), to = max(end1), length.out = dplyr::first(n_breaks))
    ) |>
    dplyr::rename(seqname = seqnames1)

  chroms_add <- data |> .calculateAddLengths()
  chroms_sub <- data |> .calculateSubtractLengths()

  tmp_adjusted <- tmp |>
    .adjustCoordinates2(chroms_add, chroms_sub, c(bin = "bin"))

  breaks <- tmp_adjusted$bin
  labels <- fmt_label(tmp$bin)

  indices <- which(!duplicated(tmp$seqname))[-1]
  labels[indices] <- paste0("\n", labels[indices])

  list(breaks = breaks, labels = labels)
}

# *--------------------------------------------------------------------------* #
# * ChromatinContacts                                                       *  #
# *--------------------------------------------------------------------------* #
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

# *--------------------------------------------------------------------------* #
# * Cooler                                                                  *  #
# *--------------------------------------------------------------------------* #
.queryCool <- function(file, path, resolution = NULL, idx = NULL, ...) {
  path <- ifelse(
    is.null(resolution),
    paste0("/", path), paste0("/resolutions/", resolution, "/", path)
  )
  rhdf5::h5read(
    file,
    name = path, index = list(idx), bit64conversion = "double", ...
  ) |>
    as.vector()
}

.getSeqinfo <- function(file, resolution = NULL) {
  chroms <- .queryCool(file, "chroms", resolution)
  Seqinfo::Seqinfo(
    seqnames = as.vector(chroms$name),
    seqlengths = as.vector(chroms$length)
  )
}

.getBins <- function(file, resolution = NULL) {
  bins <- .queryCool(file, "bins", resolution)

  grs <- GenomicRanges::GRanges(
    bins$chr,
    IRanges::IRanges(bins$start + 1, bins$end),
    bin_id = seq_along(bins$chr) - 1,
    seqinfo = .getSeqinfo(file, resolution)
  )
  names(grs) <- paste(
    GenomicRanges::seqnames(grs),
    GenomicRanges::start(grs),
    GenomicRanges::end(grs),
    sep = "_"
  )

  if ("weight" %in% names(bins)) {
    grs$weight <- as.numeric(bins$weight)
  }

  grs
}

.getPixels <- function(path, resolution, chunks) {
  tibble::tibble(
    bin1_id = .queryCool(path, "pixels/bin1_id", resolution, chunks),
    bin2_id = .queryCool(path, "pixels/bin2_id", resolution, chunks),
    count = .queryCool(path, "pixels/count", resolution, chunks)
  )
}

.getChunks <- function(
  path, resolution, bins, anchor, max_offset, same_region = FALSE
) {
  idx <- which(IRanges::overlapsAny(bins, anchor, ignore.strand = TRUE))
  offsets <- .queryCool(path, "indexes/bin1_offset", resolution, idx)
  from <- min(offsets) + 1
  to <- ifelse(
    (max(offsets) + 1) > max_offset,
    max(offsets), max(offsets) + 1
  )
  if (!same_region && (max(offsets) + 1) < max_offset) to <- to - 1
  list(
    chunks = seq(from, to, by = 1),
    idx = idx
  )
}

# *--------------------------------------------------------------------------* #
# * C                                                                        * #
# *--------------------------------------------------------------------------* #
.readPairsChromC <- function(pairs_file, chrom, inter_chrom = FALSE) {
  pairs_file <- normalizePath(pairs_file, mustWork = TRUE)

  tryCatch(
    {
      result <- .Call(
        "read_pairs_chrom_c",
        pairs_file,
        chrom,
        inter_chrom,
        PACKAGE = "gghic"
      )
      return(tibble::as_tibble(result))
    },
    error = function(e) {
      stop("C implementation failed: ", e$message)
    }
  )
}

.readPairsAllIntraC <- function(pairs_file) {
  pairs_file <- normalizePath(pairs_file, mustWork = TRUE)

  tryCatch(
    {
      result <- .Call(
        "read_pairs_all_intra_c",
        pairs_file,
        PACKAGE = "gghic"
      )
      return(tibble::as_tibble(result))
    },
    error = function(e) {
      stop("C implementation failed: ", e$message)
    }
  )
}

.read_pairs_all_c <- function(pairs_file) {
  pairs_file <- normalizePath(pairs_file, mustWork = TRUE)

  tryCatch(
    {
      result <- .Call(
        "read_pairs_all_c",
        pairs_file,
        PACKAGE = "gghic"
      )
      return(tibble::as_tibble(result))
    },
    error = function(e) {
      stop("C implementation failed: ", e$message)
    }
  )
}

# *--------------------------------------------------------------------------* #
# * buildHypergraph                                                         *  #
# *--------------------------------------------------------------------------* #
# Load and filter pairs data from file or data frame
.loadPairsData <- function(pairs, pairs_file, chrom, inter_chrom) {
  if (!is.null(pairs_file)) {
    if (!is.null(chrom)) {
      message("Reading pairs from file using C implementation...")
      pairs <- .readPairsChromC(pairs_file, chrom, inter_chrom)
    } else {
      if (inter_chrom) {
        message(
          paste0(
            "Reading all contacts (intra- and inter-chromosomal, ",
            "this may take a while)..."
          )
        )
        pairs <- .read_pairs_all_c(pairs_file)
      } else {
        message(
          "Reading all intra-chromosomal contacts (this may take a while)..."
        )
        pairs <- .readPairsAllIntraC(pairs_file)
      }
    }
  } else {
    # Filter in-memory data to chromosome(s)
    if (!is.null(chrom)) {
      if (length(chrom) == 1) {
        if (inter_chrom) {
          pairs <- pairs |>
            dplyr::filter(chrom1 == chrom | chrom2 == chrom)
        } else {
          pairs <- pairs |>
            dplyr::filter(chrom1 == chrom & chrom2 == chrom)
        }
      } else {
        if (inter_chrom) {
          pairs <- pairs |>
            dplyr::filter(chrom1 %in% chrom & chrom2 %in% chrom)
        } else {
          pairs <- pairs |>
            dplyr::filter(
              chrom1 %in% chrom & chrom2 %in% chrom & chrom1 == chrom2
            )
        }
      }
    } else {
      if (!inter_chrom) {
        pairs <- pairs |>
          dplyr::filter(chrom1 == chrom2)
      }
    }
  }

  if (nrow(pairs) == 0) {
    contact_type <- if (inter_chrom) "contacts" else "intra-chromosomal contacts"
    if (!is.null(chrom)) {
      stop("No ", contact_type, " found for ", paste(chrom, collapse = ", "))
    } else {
      stop("No ", contact_type, " found")
    }
  }

  chrom_label <- if (!is.null(chrom)) {
    if (length(chrom) == 1) chrom else paste0(length(chrom), " chromosomes")
  } else {
    "genome-wide"
  }

  message(sprintf(
    "Processing %s contacts (%s)",
    format(nrow(pairs), big.mark = ","), chrom_label
  ))

  # Remove reads with inter-chromosomal contacts if needed
  if (!inter_chrom) {
    read_chroms <- pairs |>
      dplyr::select(read_name, chrom1, chrom2) |>
      tidyr::pivot_longer(cols = c(chrom1, chrom2), values_to = "chrom") |>
      dplyr::distinct(read_name, chrom) |>
      dplyr::count(read_name, name = "n_chroms")

    reads_to_keep <- read_chroms |>
      dplyr::filter(n_chroms == 1) |>
      dplyr::pull(read_name)

    n_reads_before <- length(unique(pairs$read_name))
    pairs <- pairs |>
      dplyr::filter(read_name %in% reads_to_keep)
    n_reads_after <- length(unique(pairs$read_name))

    if (n_reads_before > n_reads_after) {
      message(sprintf(
        "Removed %d reads with inter-chromosomal contacts (%d reads remaining)",
        n_reads_before - n_reads_after, n_reads_after
      ))
    }

    if (nrow(pairs) == 0) {
      stop("No reads remain after filtering inter-chromosomal contacts")
    }
  }

  pairs
}

# Bin genomic positions and filter by contact frequency
.binAndFilterContacts <- function(pairs, bin_size, min_contacts, quantile) {
  # Bin positions
  dat <- pairs |>
    dplyr::mutate(
      bin1 = ceiling(pos1 / bin_size),
      bin2 = ceiling(pos2 / bin_size),
      bin1_id = paste0(chrom1, ":", bin1),
      bin2_id = paste0(chrom2, ":", bin2)
    ) |>
    dplyr::select(read_name, chrom1, bin1, bin1_id, chrom2, bin2, bin2_id)

  # Remove duplicate pairwise contacts within each read
  n_contacts_before <- nrow(dat)
  dat <- dat |>
    dplyr::distinct(read_name, bin1_id, bin2_id, .keep_all = TRUE)
  n_contacts_after <- nrow(dat)

  if (n_contacts_before > n_contacts_after) {
    message(sprintf(
      "Removed %s duplicate pairwise contacts within reads (%s remaining)",
      format(n_contacts_before - n_contacts_after, big.mark = ","),
      format(n_contacts_after, big.mark = ",")
    ))
  }

  # Count and filter by frequency
  counts_pw <- dat |>
    dplyr::count(chrom1, bin1_id, chrom2, bin2_id, name = "count")

  threshold <- if (!is.null(min_contacts)) {
    min_contacts
  } else {
    quantile(counts_pw$count, quantile)
  }

  message(sprintf(
    "Filtering bin pairs with >= %.0f contacts (%.0f%% quantile)",
    threshold, quantile * 100
  ))

  counts_pw_fil <- counts_pw |>
    dplyr::filter(count >= threshold)

  dat_fil <- dat |>
    dplyr::inner_join(
      counts_pw_fil |> dplyr::select(-count),
      by = c("chrom1", "bin1_id", "chrom2", "bin2_id")
    )

  if (nrow(dat_fil) == 0) {
    stop(
      "No contacts remain after filtering. Try lower quantile or min_contacts"
    )
  }

  message(sprintf(
    "Retained %s contacts from %s reads",
    format(nrow(dat_fil), big.mark = ","),
    format(length(unique(dat_fil$read_name)), big.mark = ",")
  ))

  dat_fil
}

# Build sparse incidence matrix from binned contacts
.buildIncidenceMatrix <- function(dat_fil) {
  # Create and sort bin IDs
  all_bin_ids <- unique(c(dat_fil$bin1_id, dat_fil$bin2_id))

  bin_info <- tibble::tibble(bin_id = all_bin_ids) |>
    dplyr::mutate(
      chrom = sub(":.*", "", bin_id),
      bin_num = as.integer(sub(".*:", "", bin_id)),
      chrom = factor(
        chrom,
        levels = unique(chrom)[
          order(as.numeric(gsub("\\D", "", unique(chrom))), unique(chrom))
        ]
      )
    ) |>
    dplyr::arrange(chrom, bin_num)

  all_bin_ids_sorted <- bin_info$bin_id
  read_levels <- unique(dat_fil$read_name)

  dat_indexed <- dat_fil |>
    dplyr::mutate(
      bin1_idx = match(bin1_id, all_bin_ids_sorted),
      bin2_idx = match(bin2_id, all_bin_ids_sorted),
      read_idx = match(read_name, read_levels)
    )

  # Build sparse incidence matrix
  n_bins <- length(all_bin_ids_sorted)
  n_reads <- length(read_levels)

  mat1 <- Matrix::sparseMatrix(
    i = dat_indexed$bin1_idx,
    j = dat_indexed$read_idx,
    x = rep(1, nrow(dat_indexed)),
    dims = c(n_bins, n_reads)
  )

  mat2 <- Matrix::sparseMatrix(
    i = dat_indexed$bin2_idx,
    j = dat_indexed$read_idx,
    x = rep(1, nrow(dat_indexed)),
    dims = c(n_bins, n_reads)
  )

  incidence <- (mat1 + mat2 > 0) * 1

  list(
    incidence = incidence,
    bin_ids_sorted = all_bin_ids_sorted,
    bin_info = bin_info,
    read_levels = read_levels
  )
}

# Deduplicate hyperedges and count frequencies
.deduplicateHyperedges <- function(incidence, read_levels) {
  n_bins <- nrow(incidence)
  n_reads <- ncol(incidence)

  estimated_size_gb <- (as.numeric(n_bins) * as.numeric(n_reads) * 8) / (1024^3)
  mem_limit_mb <- .getMemLim()
  available_gb <- mem_limit_mb / 1024
  threshold_gb <- available_gb * 0.5

  message(sprintf(
    "Estimated matrix size: %.1f GB, System RAM: %.1f GB (threshold: %.1f GB)",
    estimated_size_gb, available_gb, threshold_gb
  ))

  if (estimated_size_gb < threshold_gb) {
    # Fast path: dense matrix deduplication
    result <- .deduplicateDense(incidence, read_levels)
  } else {
    # Memory-efficient path: hash-based deduplication
    result <- .deduplicateSparse(incidence, read_levels, estimated_size_gb)
  }

  message(sprintf(
    "Removed %d duplicate hyperedges (%d unique patterns)",
    n_reads - length(result$hyperedge_weights), length(result$hyperedge_weights)
  ))

  result
}

# Dense matrix deduplication (fast but memory-intensive)
.deduplicateDense <- function(incidence, read_levels) {
  message("Identifying unique hyperedge patterns...")
  incidence_df <- as.data.frame(as.matrix(Matrix::t(incidence)))
  incidence_df$read_name <- read_levels

  incidence_grouped <- incidence_df |>
    dplyr::group_by(dplyr::across(-read_name)) |>
    dplyr::summarise(
      read_names = list(read_name),
      frequency = dplyr::n(),
      .groups = "drop"
    )

  incidence <- t(as.matrix(
    incidence_grouped |> dplyr::select(-read_names, -frequency)
  )) |>
    Matrix::Matrix(sparse = TRUE)

  list(
    incidence = incidence,
    hyperedge_weights = incidence_grouped$frequency,
    hyperedge_read_names = S4Vectors::SimpleList(incidence_grouped$read_names)
  )
}

# Sparse matrix deduplication (memory-efficient but slower)
.deduplicateSparse <- function(incidence, read_levels, estimated_size_gb) {
  chunck_size <- 10000
  n_reads <- ncol(incidence)

  message(sprintf(
    "Matrix too large for dense conversion (%.1f GB estimated), using hash-based deduplication...",
    estimated_size_gb
  ))

  # Pre-extract sparse matrix structure
  mat_triplet <- Matrix::summary(incidence)
  col_to_rows <- split(mat_triplet$i, mat_triplet$j)

  # Check for parallel processing
  has_furrr <- requireNamespace("furrr", quietly = TRUE)

  if (has_furrr) { # && has_future && !inherits(future::plan(), "sequential")
    message("Using parallel processing for hashing")
    hyperedge_hashes <- furrr::future_map_chr(
      seq_len(n_reads),
      function(j) paste(col_to_rows[[j]], collapse = ","),
      .options = furrr::furrr_options(seed = TRUE, chunk_size = chunck_size),
      .progress = TRUE
    )
  } else {
    if (!has_furrr) {
      warning(
        "Package 'furrr' not installed. Using sequential processing (slower). ",
        "Install furrr for parallel processing: install.packages('furrr')",
        call. = FALSE
      )
    }

    # Sequential processing with progress
    chunk_size <- chunck_size
    n_chunks <- ceiling(n_reads / chunk_size)
    hyperedge_hashes <- character(n_reads)

    for (chunk in seq_len(n_chunks)) {
      start_idx <- (chunk - 1) * chunk_size + 1
      end_idx <- min(chunk * chunk_size, n_reads)
      chunk_indices <- start_idx:end_idx

      hyperedge_hashes[chunk_indices] <- vapply(
        chunk_indices,
        function(j) paste(col_to_rows[[j]], collapse = ","),
        character(1)
      )

      message(sprintf(
        "\r  Processed %s / %s hyperedges",
        format(end_idx, big.mark = ","),
        format(n_reads, big.mark = ",")
      ))
    }
  }

  # Group by hash
  hash_df <- tibble::tibble(
    hash = hyperedge_hashes,
    read_name = read_levels,
    col_idx = seq_len(n_reads)
  )

  hash_grouped <- hash_df |>
    dplyr::group_by(hash) |>
    dplyr::summarise(
      read_names = list(read_name),
      frequency = dplyr::n(),
      first_idx = dplyr::first(col_idx),
      .groups = "drop"
    )

  incidence <- incidence[, hash_grouped$first_idx, drop = FALSE]

  list(
    incidence = incidence,
    hyperedge_weights = hash_grouped$frequency,
    hyperedge_read_names = S4Vectors::SimpleList(hash_grouped$read_names)
  )
}

# Finalize hypergraph: filter by min_multiway and calculate weight normalizations
.finalizeHypergraph <- function(
  incidence, hyperedge_weights, hyperedge_read_names,
  bin_ids_sorted, bin_info, min_multiway, bin_size, chrom
) {
  # Count contacts per hyperedge
  multiways_per_hyperedge <- Matrix::colSums(incidence)

  # Filter by minimum multi-way contacts
  keep_hyperedges <- which(multiways_per_hyperedge >= min_multiway)

  if (length(keep_hyperedges) == 0) {
    stop(sprintf(
      "No hyperedges with >= %d contacts. Try lower min_multiway", min_multiway
    ))
  }

  incidence_fil <- incidence[, keep_hyperedges, drop = FALSE]
  hyperedge_weights_fil <- hyperedge_weights[keep_hyperedges]
  hyperedge_read_names_fil <- hyperedge_read_names[keep_hyperedges]

  # Remove empty bins
  bin_sums <- Matrix::rowSums(incidence_fil)
  keep_bins <- which(bin_sums > 0)
  incidence_final <- incidence_fil[keep_bins, , drop = FALSE]

  message(sprintf(
    "Final hypergraph: %d bins, %d unique hyperedges (min %d-way contacts)",
    nrow(incidence_final), ncol(incidence_final), min_multiway
  ))

  bin_info_final <- bin_info[keep_bins, ]
  contacts_final <- multiways_per_hyperedge[keep_hyperedges]

  # Calculate normalized weights
  hyperedge_weights_df <- tibble::tibble(
    raw = hyperedge_weights_fil,
    log = log(hyperedge_weights_fil + 1)
  )

  hyperedge_weights_df <- hyperedge_weights_df |>
    dplyr::mutate(n_contacts = contacts_final) |>
    dplyr::group_by(n_contacts) |>
    dplyr::mutate(by_order = raw / sum(raw)) |>
    dplyr::ungroup() |>
    dplyr::select(-n_contacts)

  weight_range <- range(hyperedge_weights_fil)
  if (weight_range[1] == weight_range[2]) {
    hyperedge_weights_df$minmax <- 1
  } else {
    hyperedge_weights_df$minmax <- (hyperedge_weights_fil - weight_range[1]) /
      (weight_range[2] - weight_range[1])
  }

  list(
    incidence = incidence_final,
    bin_ids = bin_ids_sorted[keep_bins],
    bin_info = bin_info_final,
    hyperedge_weights = hyperedge_weights_df,
    hyperedge_reads = hyperedge_read_names_fil,
    multiways_per_hyperedge = contacts_final,
    bin_size = bin_size,
    chrom = chrom
  )

  # structure(
  #   list(
  #     incidence = incidence_final,
  #     bin_ids = bin_ids_sorted[keep_bins],
  #     bin_info = bin_info_final,
  #     hyperedge_weights = hyperedge_weights_df,
  #     hyperedge_reads = hyperedge_read_names_fil,
  #     multiways_per_hyperedge = contacts_final,
  #     bin_size = bin_size,
  #     chrom = chrom
  #   ),
  #   class = "hypergraph"
  # )
}

# *--------------------------------------------------------------------------* #
# * Add global variables to bypass check warnings.                             #
# *--------------------------------------------------------------------------* #
globalVariables(
  c(
    ".data", "bin", "bin1", "bin1_id", "bin1_idx", "bin2", "bin2_id",
    "bin2_idx", "bin_id", "bin_idx", "bin_num", "bin_size", "bin_size_kb",
    "bin_size_label", "c(cols, posCols)", "chrom", "chrom1", "chrom2",
    "chromEnd", "chromStart", "col_idx", "colour", "combn", "count", "coverage",
    "cum_len", "end1", "end2", "fill_color", "frequency", "gene_id",
    "gene_symbol", "gieStain", "group", "hash", "head", "hyperedge_idx", "i",
    "j", "median", "memory.limit", "min_bin", "n_breaks", "n_contacts",
    "n_chroms", "n_multiways", "offset", "orignal_start", "pos1", "pos2",
    "read_idx", "read_name", "read_names", "res", "score", "seqnames1",
    "seqnames2", "shape", "source_chrom", "start1", "start2", "tx_end",
    "tx_start", "type", "weight", "width", "x", "xend", "xmax", "xmin", "y",
    "y_composite", "y_end", "y_start", "yend", "ymax", "ymin"
  )
)

# *--------------------------------------------------------------------------* #
# * Public functions                                                         * #
# *--------------------------------------------------------------------------* #
.getCytobandColors <- function() {
  c(
    gneg = "grey100", stalk = "brown3", acen = "brown4",
    gpos = "grey0", gvar = "grey0",
    gpos1 = "#FFFFFF", gpos2 = "#FCFCFC", gpos3 = "#F9F9F9",
    gpos4 = "#F7F7F7", gpos5 = "#F4F4F4", gpos6 = "#F2F2F2",
    gpos7 = "#EFEFEF", gpos8 = "#ECECEC", gpos9 = "#EAEAEA",
    gpos10 = "#E7E7E7", gpos11 = "#E5E5E5", gpos12 = "#E2E2E2",
    gpos13 = "#E0E0E0", gpos14 = "#DDDDDD", gpos15 = "#DADADA",
    gpos16 = "#D8D8D8", gpos17 = "#D5D5D5", gpos18 = "#D3D3D3",
    gpos19 = "#D0D0D0", gpos20 = "#CECECE", gpos21 = "#CBCBCB",
    gpos22 = "#C8C8C8", gpos23 = "#C6C6C6", gpos24 = "#C3C3C3",
    gpos25 = "#C1C1C1", gpos26 = "#BEBEBE", gpos27 = "#BCBCBC",
    gpos28 = "#B9B9B9", gpos29 = "#B6B6B6", gpos30 = "#B4B4B4",
    gpos31 = "#B1B1B1", gpos32 = "#AFAFAF", gpos33 = "#ACACAC",
    gpos34 = "#AAAAAA", gpos35 = "#A7A7A7", gpos36 = "#A4A4A4",
    gpos37 = "#A2A2A2", gpos38 = "#9F9F9F", gpos39 = "#9D9D9D",
    gpos40 = "#9A9A9A", gpos41 = "#979797", gpos42 = "#959595",
    gpos43 = "#929292", gpos44 = "#909090", gpos45 = "#8D8D8D",
    gpos46 = "#8B8B8B", gpos47 = "#888888", gpos48 = "#858585",
    gpos49 = "#838383", gpos50 = "#808080", gpos51 = "#7E7E7E",
    gpos52 = "#7B7B7B", gpos53 = "#797979", gpos54 = "#767676",
    gpos55 = "#737373", gpos56 = "#717171", gpos57 = "#6E6E6E",
    gpos58 = "#6C6C6C", gpos59 = "#696969", gpos60 = "#676767",
    gpos61 = "#646464", gpos62 = "#616161", gpos63 = "#5F5F5F",
    gpos64 = "#5C5C5C", gpos65 = "#5A5A5A", gpos66 = "#575757",
    gpos67 = "#545454", gpos68 = "#525252", gpos69 = "#4F4F4F",
    gpos70 = "#4D4D4D", gpos71 = "#4A4A4A", gpos72 = "#484848",
    gpos73 = "#454545", gpos74 = "#424242", gpos75 = "#404040",
    gpos76 = "#3D3D3D", gpos77 = "#3B3B3B", gpos78 = "#383838",
    gpos79 = "#363636", gpos80 = "#333333", gpos81 = "#303030",
    gpos82 = "#2E2E2E", gpos83 = "#2B2B2B", gpos84 = "#292929",
    gpos85 = "#262626", gpos86 = "#242424", gpos87 = "#212121",
    gpos88 = "#1E1E1E", gpos89 = "#1C1C1C", gpos90 = "#191919",
    gpos91 = "#171717", gpos92 = "#141414", gpos93 = "#121212",
    gpos94 = "#0F0F0F", gpos95 = "#0C0C0C", gpos96 = "#0A0A0A",
    gpos97 = "#070707", gpos98 = "#050505", gpos99 = "#020202",
    gpos100 = "#000000"
  )
}
