StatAnnotation <- ggplot2::ggproto(
  "StatAnnotation",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  setup_params = function(data, params) params,
  compute_panel = function(
    data, scales,
    txdb = NULL, tx2gene = NULL, gtf_path = NULL, width_ratio = 1 / 50,
    spacing_ratio = 1 / 3, maxgap = -1, include_ncrna = TRUE,
    style = c("basic", "arrow"), gene_symbols = NULL, chrom_prefix = TRUE
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
    style <- match.arg(style)
    if (is.null(gtf_path) && (is.null(txdb) || is.null(tx2gene))) {
      stop("Either gtf_path or txdb and tx2gene must be provided.")
    }

    name_pkg <- .getPkgName()
    env <- get(".env", envir = asNamespace(name_pkg))
    n_annotation <- env$n_annotation
    n_track <- env$n_track
    n_concatemer <- env$n_concatemer
    if (env$n_hic == 1) {
      max_y <- env$max_y
      max_x <- env$max_x
      min_x <- env$min_x
      res <- env$res
      n_sn <- env$n_sn
    } else {
      dat_hic <- data |>
        .calculateHicCoordinates()
      max_y <- max(dat_hic$ymax, na.rm = TRUE)
      max_x <- max(dat_hic$xend, na.rm = TRUE)
      min_x <- min(dat_hic$xmin, na.rm = TRUE)
      res <- data$end1[1] - data$start1[1] + 1
      n_sn <- length(unique(c(data$seqnames1, data$seqnames2)))
    }
    if (n_sn > 1) {
      if (env$n_hic == 1 || env$n_track == 1 || env$n_concatemer == 1) {
        maxs_x <- env$maxs_x
      } else {
        maxs_x <- dat_hic |>
          dplyr::group_by(seqnames2) |>
          dplyr::summarize(maxs_x = max(xend)) |>
          dplyr::pull(maxs_x) |>
          stats::setNames(unique(dat_hic$seqnames2))
      }
    }
    min_y <- ifelse(
      n_track > n_annotation || n_concatemer > n_annotation, env$min_y, 0
    )

    .height <- max_y * width_ratio
    .height_cds <- .height * 0.8

    genes <- .retriveGenes(
      data, txdb, tx2gene, gtf_path, maxgap, include_ncrna, gene_symbols,
      chrom_prefix
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

    if (style == "basic") {
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
    if (style == "arrow") {
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
        .calculateAddLengths()
      chroms_sub <- data |>
        .calculateSubtractLengths()

      dat <- dat |>
        .adjustCoordinates2(chroms_add, chroms_sub, c(x = "x", xmax = "xmax"))

      if (style == "arrow") {
        dat <- dat |>
          .adjustCoordinates2(chroms_add, chroms_sub, c(xend = "xend"))
      }
    }

    env$min_y <- min(dat$ymin) - .height_cds * 2
    env$n_annotation <- env$n_annotation + 1

    if (n_sn > 1) {
      dat_vline <- tibble::tibble(
        xmax = maxs_x[-length(maxs_x)],
        ymin = env$min_y - (.height * spacing_ratio),
        y = min_y - (.height * spacing_ratio),
        feature = "vline",
        x = xmax,
        gene_symbol = "",
        gene_id = "",
        strand = "*",
        seqname = names(maxs_x)[-length(maxs_x)]
      )

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
  default_aes = ggplot2::aes(
    colour = "#48CFCB", fill = "#48CFCB", fontsize = 10, linetype = "dashed"
  ),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord, fontsize = 10, style = c("basic", "arrow"),
    colour = "#48CFCB", fill = "#48CFCB", draw_boundary = TRUE,
    boundary_colour = "black", linetype = "dashed"
  ) {
    style <- match.arg(style)
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

    if (style == "basic") {
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

      intron_mid <- (coords_intron$x + coords_intron$xmax) / 2
      intron_length <- coords_intron$xmax - coords_intron$x
      # 25% of intron or fixed small distance
      arrow_offset <- pmin(intron_length * 0.25, 1 / 10)

      arrow_x0 <- ifelse(
        coords_intron$strand == "+",
        intron_mid - arrow_offset,
        intron_mid + arrow_offset
      )
      arrow_x1 <- ifelse(
        coords_intron$strand == "+",
        intron_mid + arrow_offset,
        intron_mid - arrow_offset
      )

      grob_intron_line <- grid::segmentsGrob(
        x0 = coords_intron$x, x1 = coords_intron$xmax,
        y0 = coords_intron$y, y1 = coords_intron$y,
        gp = grid::gpar(col = colour),
        default.units = "native"
      )

      lengths_arrow <- pmin(arrow_offset, 1 / 120)
      grob_intron_arrow <- grid::segmentsGrob(
        x0 = arrow_x0, x1 = arrow_x1, y0 = coords_intron$y,
        y1 = coords_intron$y, arrow = grid::arrow(
          type = "open", length = grid::unit(lengths_arrow, "native"),
          ends = "last"
        ), gp = grid::gpar(col = colour), default.units = "native"
      )
      grids <- grid::gList(
        grob_exon, grob_intron_line, grob_intron_arrow, grob_text
      )
    }
    if (style == "arrow") {
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
#' # Load Hi-C data
#' cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
#'   import()
#'
#' gtf_file <- "path/to/genes.gtf"
#' gghic(cc) + geom_annotation(gtf_path = gtf_file)
#'
#' # Filter specific genes
#' gghic(cc) +
#'   geom_annotation(
#'     gtf_path = gtf_file,
#'     gene_symbols = c("BRCA1", "TP53", "MYC")
#'   )
#'
#' # Arrow style with custom colors
#' gghic(cc) +
#'   geom_annotation(
#'     gtf_path = gtf_file,
#'     style = "arrow",
#'     colour = "darkblue",
#'     fill = "lightblue",
#'     fontsize = 8
#'   )
#'
#' # Exclude non-coding RNAs
#' gghic(cc) +
#'   geom_annotation(
#'     gtf_path = gtf_file,
#'     include_ncrna = FALSE
#'   )
#' }
#' @export
#' @aliases geom_annotation
geom_annotation <- function(
  mapping = NULL, data = NULL, stat = StatAnnotation, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, txdb = NULL,
  tx2gene = NULL, gtf_path = NULL, width_ratio = 1 / 50, spacing_ratio = 1 / 3,
  maxgap = -1, include_ncrna = TRUE, style = c("basic", "arrow"),
  gene_symbols = NULL, chrom_prefix = TRUE, fontsize = 10,
  colour = "#48CFCB", fill = "#48CFCB", draw_boundary = TRUE,
  boundary_colour = "black", linetype = "dashed", ...
) {
  ggplot2::layer(
    geom = GeomAnnotation, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(
      na.rm = na.rm, txdb = txdb, tx2gene = tx2gene, gtf_path = gtf_path,
      width_ratio = width_ratio, spacing_ratio = spacing_ratio, maxgap = maxgap,
      include_ncrna = include_ncrna, style = style, gene_symbols = gene_symbols,
      chrom_prefix = chrom_prefix, fontsize = fontsize, colour = colour,
      fill = fill, draw_boundary = draw_boundary,
      boundary_colour = boundary_colour, linetype = linetype, ...
    )
  )
}
