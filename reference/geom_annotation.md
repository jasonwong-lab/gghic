# geom_annotation

A ggplot2 geom for gene model tracks.

## Usage

``` r
geom_annotation(
  mapping = NULL,
  data = NULL,
  stat = StatAnnotation,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  txdb = NULL,
  tx2gene = NULL,
  gtf_path = NULL,
  width_ratio = 1/50,
  spacing_ratio = 1/3,
  maxgap = -1,
  include_ncrna = TRUE,
  style = c("basic", "arrow"),
  gene_symbols = NULL,
  chrom_prefix = TRUE,
  fontsize = 10,
  colour = "#48CFCB",
  fill = "#48CFCB",
  draw_boundary = TRUE,
  boundary_colour = "black",
  linetype = "dashed",
  ...
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). If
  specified and `inherit.aes = TRUE` (the default), it is combined with
  the default mapping at the top level of the plot. You must supply
  `mapping` if there is no plot mapping.

- data:

  The data to be displayed in this layer. There are three options:

  If `NULL`, the default, the data is inherited from the plot data as
  specified in the call to
  [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

  A `data.frame`, or other object, will override the plot data. All
  objects will be fortified to produce a data frame. See
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
  for which variables will be created.

  A `function` will be called with a single argument, the plot data. The
  return value must be a `data.frame`, and will be used as the layer
  data. A `function` can be created from a `formula` (e.g.
  `~ head(.x, 10)`).

- stat:

  The statistical transformation to use on the data for this layer. When
  using a `geom_*()` function to construct a layer, the `stat` argument
  can be used to override the default coupling between geoms and stats.
  The `stat` argument accepts the following:

  - A `Stat` ggproto subclass, for example `StatCount`.

  - A string naming the stat. To give the stat as a string, strip the
    function name of the `stat_` prefix. For example, to use
    `stat_count()`, give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- position:

  A position adjustment to use on the data for this layer. This can be
  used in various ways, including to prevent overplotting and improving
  the display. The `position` argument accepts the following:

  - The result of calling a position function, such as
    `position_jitter()`. This method allows for passing extra arguments
    to the position.

  - A string naming the position adjustment. To give the position as a
    string, strip the function name of the `position_` prefix. For
    example, to use `position_jitter()`, give the position as
    `"jitter"`.

  - For more information and other ways to specify the position, see the
    [layer
    position](https://ggplot2.tidyverse.org/reference/layer_positions.html)
    documentation.

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

- show.legend:

  logical. Should this layer be included in the legends? `NA`, the
  default, includes if any aesthetics are mapped. `FALSE` never
  includes, and `TRUE` always includes. It can also be a named logical
  vector to finely select the aesthetics to display. To include legend
  keys for all levels, even when no data exists, use `TRUE`. If `NA`,
  all levels are shown in legend, but unobserved levels are omitted.

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them. This is most useful for helper functions that define both
  data and aesthetics and shouldn't inherit behaviour from the default
  plot specification, e.g.
  [`annotation_borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

- txdb:

  The TxDb object. Default is `NULL`.

- tx2gene:

  An optional data frame or tibble that maps transcript information to
  gene information. It should include the following columns:

  - chrom: Chromosome number or name.

  - gene_id: Entrez gene ID.

  - gene_symbol: Common symbol or name of the gene.

  - tx_id: Entrez transcript ID.

  - tx_name: Name of the transcript.

  - gene_type: Type or classification of the gene.

  - tx_type: Type or classification of the transcript.

- gtf_path:

  The path to the GTF file, which is used to generate `txdb` and
  `tx2gene`. Generated files are saved in the cache directory. Default
  is `NULL`.

- width_ratio:

  The ratio of the width of each gene model track relative to the height
  of the Hi-C plot. Default is `1/50`.

- spacing_ratio:

  The ratio of the spacing between two gene model tracks. Default is
  `1/3`.

- maxgap:

  The maximum gap between genes to be drawn in the same line. Default is
  `-1`.

- include_ncrna:

  Whether to include ncRNA or not. Default is `TRUE`.

- style:

  The style of the gene model track, which can be `"basic"` or
  `"arrow"`. Default is `"basic"`.

- gene_symbols:

  A character vector of gene symbols to be included only. Default is
  `NULL`.

- chrom_prefix:

  Whether the input data has chromosome names with prefix 'chr' or not.
  Default is `TRUE`.

- fontsize:

  The font size of the gene symbols. Default is `10`.

- colour:

  The color of the gene model track. Default is `"#48CFCB"`.

- fill:

  The fill color of the gene model track. Default is `"#48CFCB"`.

- draw_boundary:

  Logical. Draw boundary lines between chromosomes in multi-chromosome
  plots. Default is `TRUE`.

- boundary_colour:

  Character. Color for boundary lines. Default is `"black"`.

- linetype:

  Line type for boundaries. Default is `"dashed"`.

- ...:

  Parameters to be ignored.

## Value

A ggplot object.

## Details

Requires the following aesthetics:

- seqnames1

- start1

- end1

- seqnames2

- start2

- end2

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

gtf_file <- "path/to/genes.gtf"
gghic(cc) + geom_annotation(gtf_path = gtf_file)

# Filter specific genes
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    gene_symbols = c("BRCA1", "TP53", "MYC")
  )

# Arrow style with custom colors
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    style = "arrow",
    colour = "darkblue",
    fill = "lightblue",
    fontsize = 8
  )

# Exclude non-coding RNAs
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    include_ncrna = FALSE
  )
} # }
```
