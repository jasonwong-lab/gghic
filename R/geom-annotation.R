extract_trs <- function(
  grs, genes, txdump, include_ncrna, gene_info
) {
  `!!` <- rlang::`!!`

  chrom <- as.character(grs@seqnames)
  start <- grs@ranges@start
  end <- start + grs@ranges@width - 1

  r <- range(
    genes[
      S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(grs, genes, ignore.strand = TRUE)
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
        .txdump$transcripts$tx_start < GenomicRanges::end(r) &
        .txdump$transcripts$tx_end > GenomicRanges::start(r), ,
      drop = FALSE
    ]
  }
  if (!is.null(gene_info)) {
    gene_info <- gene_info |>
      dplyr::filter(
        chrom == !!chrom,
        start < GenomicRanges::end(r),
        end > GenomicRanges::start(r)
      )
    if (nrow(gene_info) == 0) {
      return(GenomicRanges::GRangesList())
    }

    .txdump$transcripts <- .txdump$transcripts[
      .txdump$transcripts$tx_name %in% gene_info$tx_id, ,
      drop = FALSE
    ]
  }

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
    chromosome = chrom, start = GenomicRanges::start(r), end = GenomicRanges::end(r),
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
    tmp <- GenomicRanges::pintersect(c(tr, introns), grs, ignore.strand = TRUE)
    tmp[width(tmp) > 0]
  })

  GenomicRanges::GRangesList(trs_exon_intron)
}

decide_track_lines <- function(genes_span, maxgap) {
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
        genes_n1o, genes_n1, maxgap = maxgap, ignore.strand = TRUE
      )
    ]
    n1o_nonol <- names(genes_n1o_nonol)
    if (length(n1o_nonol) == 0) {
      dat_line$line[dat_line$gene_id == n1] <- max(dat_line$line) + 1
      names_gene_copy <- setdiff(names_gene_copy, n1)
      next
    }
    dat_line$line[dat_line$gene_id %in% c(n1, n1o_nonol)] <- max(dat_line$line) + 1
    names_gene_copy <- setdiff(names_gene_copy, n1o_nonol)
  }
  if (sum(dat_line$line == 0)) {
    dat_line$line[dat_line$line == 0] <- max(dat_line$line) + 1
  }

  dat_line
}

retrive_genes <- function(
  data, txdb, tx2gene, gtf_path, maxgap, include_ncrna,
  gene_symbols, chrom_prefix
) {
  txdb <- ensure_txdb(txdb, gtf_path)
  tx2gene <- ensure_tx2gene(tx2gene, gtf_path)
  txdump <- ensure_txdump(txdb, gtf_path)

  info_genes <- NULL
  if (!is.null(gene_symbols)) {
    info_genes <- tx2gene |>
      dplyr::filter(gene_symbol %in% gene_symbols)
  }

  grs <- purrr::map(
    unique(data$seqnames1), function(.x) {
      tmp <- .x
      if (!chrom_prefix) tmp <- paste0("chr", tmp)
      GenomicRanges::GRanges(
        seqnames = tmp,
        ranges = IRanges::IRanges(
          start = min(data$start1[data$seqnames1 == .x]),
          end = max(data$end1[data$seqnames1 == .x])
        )
      )
    }
  )

  genes <- purrr::map(
    grs, function(.x) {
      GenomicFeatures::genes(
        txdb, columns = c("exon_id"), filter = list(tx_chrom = .x@seqnames)
      )
    }
  )

  trs <- purrr::map2(
    grs, genes,
    extract_trs,
    txdump = txdump, include_ncrna = include_ncrna, gene_info = info_genes
  )

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

      dat_line <- decide_track_lines(genes_span, maxgap)

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
          ), by = "gene_id"
        )

      dat_gene
    }
  )

  dplyr::bind_rows(dat)
}

StatAnnotation <- ggplot2::ggproto(
  "StatAnnotation",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "txdb", "tx2gene", "gtf_path", "width_ratio", "spacing_ratio", "maxgap",
    "include_ncrna", "style", "gene_symbols", "chrom_prefix"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_panel = function(
    data, scales,
    txdb, tx2gene, gtf_path, width_ratio, spacing_ratio, maxgap, include_ncrna,
    style, gene_symbols, chrom_prefix
  ) {
    # ======================================================================== #
    #   ^                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, 0)                                                               #
    #   |                      +------------------------------------+          #
    #   | +-----------+        | (x, y)                   (xmax, y) |          #
    #   | |           |<-------|                                    |          #
    #   | +-----------+        | (x, ymin)             (xmax, ymin) |          #
    #   |                      +------------------------------------+          #
    #   |                          [gene name]                                 #
    #   |    , -----------------------------------------------------+          #
    #   |  ,   (x, y)                                     (xmax, y) |          #
    #   | + (xend, yend)                                            |          #
    #   |  `   (x, ymin)                               (xmax, ymin) |          #
    #   |    ` -----------------------------------------------------+          #
    #   |                          [gene name]                                 #
    #   | +--------------------------------------------------- `               #
    #   | | (x, y)                                      (xmax, y) `            #
    #   | |                                            (xend, yend) +          #
    #   | | (x, ymin)                                (xmax, ymin) ,            #
    #   | +---------------------------------------------------- ,              #
    #   |                          [gene name]                                 #
    # ======================================================================== #
    if (is.null(gtf_path) && (is.null(txdb) || is.null(tx2gene))) {
      stop("Either gtf_path or txdb and tx2gene must be provided.")
    }

    name_pkg <- get_pkg_name()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    if (env$n_hic == 1) {
      max_y <- env$max_y
      max_x <- env$max_x
      min_x <- env$min_x
      res <- env$res
      n_sn <- env$n_sn
    } else {
      dat_hic <- data |>
        calculate_hic_coordinates()
      max_y <- max(dat_hic$ymax, na.rm = TRUE)
      max_x <- max(dat_hic$xend, na.rm = TRUE)
      min_x <- min(dat_hic$xmin, na.rm = TRUE)
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
    }
    if (n_sn > 1) {
      if (env$n_hic == 1 || env$n_track == 1) {
        maxs_x <- env$maxs_x
      } else {
        maxs_x <- dat_hic |>
          dplyr::group_by(seqnames2) |>
          dplyr::summarize(maxs_x = max(xend)) |>
          dplyr::pull(maxs_x) |>
          stats::setNames(unique(dat_hic$seqnames2))
      }
    }
    min_y <- ifelse(n_track > n_annotation, env$min_y, 0)

    .height <- max_y * width_ratio
    .height_cds <- .height * 0.8

    genes <- retrive_genes(
      data, txdb, tx2gene, gtf_path, maxgap, include_ncrna,
      gene_symbols, chrom_prefix
    ) |>
      dplyr::rename(seqname = seqnames)

    ys <- genes |>
      dplyr::distinct(line) |>
      dplyr::mutate(
        y = min_y -
          (.height * (dplyr::row_number() - 1)) -
          (dplyr::row_number() * (.height * spacing_ratio)) -
          res / 2
      )

    dat_text <- genes |>
      dplyr::left_join(ys, by = "line") |>
      dplyr::group_by(gene_symbol) |>
      dplyr::summarise(
        gene_symbol = dplyr::first(gene_symbol),
        gene_id = dplyr::first(gene_id),
        line = dplyr::first(line),
        seqname = dplyr::first(seqname),
        strand = dplyr::first(strand),
        start = min(start),
        end = max(end),
        width = end - start,
        x = mean(c(start, end)),
        y = dplyr::first(y),
        ymin = y - .height_cds - (.height_cds / 5),
        feature = "text",
        xmax = end
      )

    if (style[1] == "basic") {
      dat_cds <- genes |>
        dplyr::filter(feature == "CDS") |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds
        )
      dat_intron <- genes |>
        dplyr::filter(feature == "intron") |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds,
          y = 0.5 * (y + ymin)
        )
      dat_others <- genes |>
        dplyr::filter(feature %in% c("utr3", "utr5", "ncRNA")) |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds,
          y = y - (.height_cds / 8),
          ymin = ymin + (.height_cds / 8)
        )
      dat <- dplyr::bind_rows(dat_cds, dat_intron, dat_others, dat_text)
    }
    if (style[1] == "arrow") {
      dat <- genes |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::group_by(gene_symbol) |>
        dplyr::summarise(
          gene_symbol = dplyr::first(gene_symbol),
          gene_id = dplyr::first(gene_id),
          line = dplyr::first(line),
          seqname = dplyr::first(seqname),
          strand = dplyr::first(strand),
          y = dplyr::first(y),
          xend_p = max(end),
          xend_n = min(start),
          xmax_n = max(end),
          x_p = min(start),
          x_n = (xmax_n - xend_n) * 0.1 + xend_n,
          xmax_p = xend_p - ((xend_p - x_p) * 0.1),
          xend = dplyr::case_when(
            strand == "+" ~ xend_p,
            strand == "-" ~ xend_n
          ),
          xmax = dplyr::case_when(
            strand == "+" ~ xmax_p,
            strand == "-" ~ xmax_n
          ),
          x = dplyr::case_when(
            strand == "+" ~ x_p,
            strand == "-" ~ x_n
          ),
          ymin = y - .height_cds,
          yend = (y + ymin) / 2,
          feature = "arrow"
        ) |>
        dplyr::bind_rows(dat_text)
    }

    if (!chrom_prefix) {
      dat <- dat |>
        dplyr::mutate(seqname = stringr::str_remove(seqname, "^chr"))
    }

    if ((n_sn > 1 || (n_sn == 2 && any(data$seqnames1 == data$seqnames2)))) {
      chroms_add <- data |>
        calculate_add_lengths()
      chroms_sub <- data |>
        calculate_subtract_lengths()

      dat <- dat |>
        adjust_coordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))

      if (style[1] == "arrow") {
        dat <- dat |>
          adjust_coordinates2(chroms_add, chroms_sub, c(xend = "xend"))
      }
    }

    env$min_y <- min(dat$ymin)
    env$n_annotation <- env$n_annotation + 1

    if (n_sn > 1) {
      dat_vline <- dat |>
        dplyr::slice(seq_len(length(maxs_x))) |>
        dplyr::mutate(
          xmax = maxs_x,
          ymin = env$min_y - (.height * spacing_ratio),
          y = min_y - (.height * spacing_ratio),
          feature = "vline"
        ) |>
        dplyr::slice(seq_len((dplyr::n() - 1)))

      dat <- dplyr::bind_rows(dat, dat_vline)
    }

    dat
  }
)

GeomAnnotation <- ggplot2::ggproto(
  "GeomAnnotation",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "feature", "gene_symbol", "strand", "seqname"
  ),
  extra_params = c(
    ggplot2::Geom$extra_params, "fontsize", "style", "colour", "fill",
    "draw_boundary", "boundary_colour", "linetype"
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fontsize, style, colour, fill, draw_boundary, boundary_colour, linetype
  ) {
    coords <- coord$transform(data, panel_params)

    coords_text <- coords |>
      dplyr::filter(feature == "text")
    grob_text <- grid::textGrob(
      label = coords_text$gene_symbol,
      x = coords_text$x, y = coords_text$ymin,
      just = c("center", "top"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )

    if (style[1] == "basic") {
      coords_exon <- coords |>
        dplyr::filter(feature %in% c("CDS", "utr3", "utr5", "ncRNA"))
      grob_exon <- grid::polygonGrob(
        x = c(coords_exon$x, coords_exon$x, coords_exon$xmax, coords_exon$xmax),
        y = c(coords_exon$y, coords_exon$ymin, coords_exon$ymin, coords_exon$y),
        id = rep(seq_len(nrow(coords_exon)), 4),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )

      coords_intron <- coords |>
        dplyr::filter(feature == "intron")
      ends <- ifelse(coords_intron$strand == "+", "last", "first")
      lengths_arrow <- rep(1 / 90, nrow(coords_intron))
      lengths_intron <- coords_intron$xmax - coords_intron$x
      lengths_arrow[lengths_arrow > lengths_intron] <- 0
      grob_intron <- grid::segmentsGrob(
        x0 = coords_intron$x, x1 = coords_intron$xmax,
        y0 = coords_intron$y, y1 = coords_intron$y,
        arrow = grid::arrow(
          type = "open",
          length = grid::unit(lengths_arrow, "native"),
          ends = ends
        ),
        gp = grid::gpar(col = colour),
        default.units = "native"
      )

      grids <- grid::gList(grob_exon, grob_intron, grob_text)
    }
    if (style[1] == "arrow") {
      coords_arrow_p <- coords |>
        dplyr::filter(feature == "arrow", strand == "+")
      grob_arrow_p <- grid::nullGrob()
      if (nrow(coords_arrow_p) > 0) {
        grob_arrow_p <- grid::polygonGrob(
          x = c(
            coords_arrow_p$x, coords_arrow_p$x, coords_arrow_p$xmax,
            coords_arrow_p$xend, coords_arrow_p$xmax
          ),
          y = c(
            coords_arrow_p$y, coords_arrow_p$ymin, coords_arrow_p$ymin,
            coords_arrow_p$yend, coords_arrow_p$y
          ),
          id = rep(seq_len(nrow(coords_arrow_p)), 5),
          gp = grid::gpar(col = colour, fill = fill),
          default.units = "native"
        )
      }

      coords_arrow_n <- coords |>
        dplyr::filter(feature == "arrow", strand == "-")
      grob_arrow_n <- grid::nullGrob()
      if (nrow(coords_arrow_n) > 0) {
        grob_arrow_n <- grid::polygonGrob(
          x = c(
            coords_arrow_n$x, coords_arrow_n$xend, coords_arrow_n$x,
            coords_arrow_n$xmax, coords_arrow_n$xmax
          ),
          y = c(
            coords_arrow_n$y, coords_arrow_n$yend, coords_arrow_n$ymin,
            coords_arrow_n$ymin, coords_arrow_n$y
          ),
          id = rep(seq_len(nrow(coords_arrow_n)), 5),
          gp = grid::gpar(col = colour, fill = fill),
          default.units = "native"
        )
      }

      grids <- grid::gList(grob_arrow_p, grob_arrow_n, grob_text)
    }

    grob_vline <- grid::nullGrob()
    if (length(unique(coords$seqname)) > 1 && draw_boundary) {
      coords_vline <- coords |>
        dplyr::filter(feature == "vline")
      grob_vline <- grid::polylineGrob(
        x = c(coords_vline$xmax, coords_vline$xmax),
        y = c(coords_vline$ymin, coords_vline$y),
        id = rep(seq_len(nrow(coords_vline)), 2),
        gp = grid::gpar(col = boundary_colour, lty = linetype),
        default.units = "native"
      )
    }

    grids <- grid::gList(grids, grob_vline)

    grids
  }
)

#' geom_annotation
#'
#' @description A ggplot2 geom for gene model tracks.
#' @inheritParams ggplot2::geom_polygon
#' @inheritParams geom_hic
#' @param txdb The TxDb object. Default is `NULL`.
#' @param tx2gene An optional data frame or tibble that maps transcript
#'   information to gene information. It should include the following columns:
#'   * chrom: Chromosome number or name.
#'   * gene_id: Entrez gene ID.
#'   * gene_symbol: Common symbol or name of the gene.
#'   * tx_id: Entrez transcript ID.
#'   * tx_name: Name of the transcript.
#'   * gene_type: Type or classification of the gene.
#'   * tx_type: Type or classification of the transcript.
#' @param gtf_path The path to the GTF file, which is used to generate `txdb`
#'   and `tx2gene`. Generated files are saved in the cache directory.
#'   Default is `NULL`.
#' @param chrom_prefix Whether the input data has chromosome names
#'   with prefix 'chr' or not. Default is `TRUE`.
#' @param width_ratio The ratio of the width of each gene model track
#'   relative to the height of the Hi-C plot. Default is `1/50`.
#' @param spacing_ratio The ratio of the spacing between two gene model tracks.
#'   Default is `1/3`.
#' @param maxgap The maximum gap between genes to be drawn in the same line.
#'   Default is `-1`.
#' @param include_ncrna Whether to include ncRNA or not. Default is `TRUE`.
#' @param style The style of the gene model track, which can be `"basic"`
#'   or `"arrow"`. Default is `"basic"`.
#' @param gene_symbols A character vector of gene symbols to be included only.
#'   Default is `NULL`.
#' @param fontsize The font size of the gene symbols. Default is `10`.
#' @param colour The color of the gene model track. Default is `"#48CFCB"`.
#' @param fill The fill color of the gene model track. Default is `"#48CFCB"`.
#' @param ... Parameters to be ignored.
#' @details
#' Requires the following aesthetics:
#' * seqnames1
#' * start1
#' * end1
#' * seqnames2
#' * start2
#' * end2
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' library(gghic)
#' library(ggplot2)
#' library(dplyr)
#' library(HiCExperiment)
#' library(InteractionSet)
#' library(scales)
#' library(glue)
#' library(rappdirs)
#'
#' download_example_files()
#' dir_cache_gghic <- user_cache_dir(appname = "gghic")
#'
#' hic <- glue("{dir_cache_gghic}/chr4_11-5kb.cool") |>
#'   CoolFile() |>
#'   import(cf)
#'
#' gis <- interactions(hic)
#' gis$score <- log10(gis$balanced)
#' x <- as_tibble(gis)
#' scores <- x$score[pairdist(gis) != 0 & !is.na(pairdist(gis) != 0)]
#' scores <- scores[!is.na(scores) & !is.infinite(scores)]
#' x$score <- oob_squish(x$score, c(min(scores), max(scores)))
#'
#' p <- x |>
#'   filter(
#'     seqnames1 == "chr11", seqnames2 == "chr11",
#'     center1 > 67000000, center1 < 67100000,
#'     center2 > 67000000, center2 < 67100000
#'   ) |>
#'   ggplot(
#'     aes(
#'       seqnames1 = seqnames1, start1 = start1, end1 = end1,
#'       seqnames2 = seqnames2, start2 = start2, end2 = end2,
#'       fill = score
#'     )
#'   ) +
#'   geom_hic()
#'
#' path_gtf <- glue("{dir_cache_gghic}/gencode-chr4_11.gtf.gz")
#'
#' p + geom_annotation(gtf_path = path_gtf, style = "basic", maxgap = 100000)
#' }
#' @export geom_annotation
#' @aliases geom_annotation
geom_annotation <- function(
  mapping = NULL, data = NULL, stat = StatAnnotation, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...,
  txdb = NULL, tx2gene = NULL, gtf_path = NULL, width_ratio = 1 / 50,
  spacing_ratio = 1 / 3, maxgap = -1, include_ncrna = TRUE,
  style = c("basic", "arrow"), gene_symbols = NULL, chrom_prefix = TRUE,
  fontsize = 10, colour = "#48CFCB", fill = "#48CFCB",
  draw_boundary = TRUE, boundary_colour = "black", linetype = "dashed"
) {
  ggplot2::layer(
    geom = GeomAnnotation, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(
      na.rm = na.rm, ...,
      txdb = txdb, tx2gene = tx2gene, gtf_path = gtf_path,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio, maxgap = maxgap,
      include_ncrna = include_ncrna, style = style, gene_symbols = gene_symbols,
      chrom_prefix = chrom_prefix,
      fontsize = fontsize, colour = colour, fill = fill,
      draw_boundary = draw_boundary, boundary_colour = boundary_colour,
      linetype = linetype
    )
  )
}
