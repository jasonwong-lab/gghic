# Add gene annotation tracks to Hi-C plots

Displays gene models as annotation tracks below Hi-C heatmaps, showing
gene structure with exons, introns, UTRs, and strand orientation. Gene
annotations are automatically retrieved from GTF files or TxDb objects
and positioned below the Hi-C plot with gene symbols labeled underneath.
The function supports two visualization styles: "basic" (detailed
exon/intron structure) and "arrow" (simplified arrow representation
indicating strand direction).

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

  A TxDb object from the GenomicFeatures package containing transcript
  annotations. Must be provided together with `tx2gene` if `gtf_path` is
  not specified. Default is `NULL`.

- tx2gene:

  A data frame or tibble mapping transcript information to gene
  information. Required if using `txdb` parameter. Should include the
  following columns:

  - `chrom`: Chromosome number or name

  - `gene_id`: Entrez gene ID or unique gene identifier

  - `gene_symbol`: Common gene symbol or name (e.g., "TP53", "BRCA1")

  - `tx_id`: Entrez transcript ID or unique transcript identifier

  - `tx_name`: Name of the transcript

  - `gene_type`: Gene biotype (e.g., "protein_coding", "lncRNA")

  - `tx_type`: Transcript biotype Default is `NULL`.

- gtf_path:

  Character string specifying the path to a GTF/GFF file. The function
  will automatically parse the GTF file to generate `txdb` and `tx2gene`
  objects. Parsed data is cached to speed up subsequent calls. Either
  `gtf_path` or both `txdb` and `tx2gene` must be provided. Default is
  `NULL`.

- width_ratio:

  Numeric value controlling the height of each gene track relative to
  the Hi-C plot height. Smaller values create thinner tracks. Default is
  `1/50` (2% of Hi-C plot height).

- spacing_ratio:

  Numeric value controlling the vertical spacing between gene tracks as
  a proportion of track height. Larger values increase spacing. Default
  is `1/3` (33% of track height).

- maxgap:

  Integer specifying the maximum genomic distance (in bp) between genes
  on the same horizontal line. Genes within this distance will be placed
  on separate lines to prevent overlap. Set to `-1` for automatic
  spacing. Default is `-1`.

- include_ncrna:

  Logical indicating whether to include non-coding RNA genes (lncRNA,
  miRNA, etc.) in the annotation track. Set to `FALSE` to show only
  protein-coding genes. Default is `TRUE`.

- style:

  Character string specifying the visualization style. Options:

  - `"basic"`: Shows detailed gene structure with distinct exons (thick
    boxes), introns (thin lines with directional arrows), and UTRs (thin
    boxes)

  - `"arrow"`: Simplified representation showing each gene as a single
    arrow pointing in the direction of transcription (5' to 3') Default
    is `"basic"`.

- gene_symbols:

  Character vector of specific gene symbols to display. When provided,
  only these genes will be shown in the annotation track. Useful for
  highlighting genes of interest. Default is `NULL` (show all genes).

- chrom_prefix:

  Logical indicating whether chromosome names include the "chr" prefix
  (e.g., "chr1" vs "1"). Should match the naming convention in your Hi-C
  data. Default is `TRUE`.

- fontsize:

  Numeric value specifying the font size for gene symbol labels. Default
  is `10`.

- colour:

  Character string specifying the outline color for gene features.
  Default is `"#48CFCB"` (teal).

- fill:

  Character string specifying the fill color for gene features (exons,
  UTRs, arrows). Default is `"#48CFCB"` (teal).

- draw_boundary:

  Logical indicating whether to draw vertical boundary lines between
  chromosomes in multi-chromosome displays. Default is `TRUE`.

- boundary_colour:

  Character string specifying the color of chromosome boundary lines.
  Default is `"black"`.

- linetype:

  Character string or integer specifying the line type for chromosome
  boundaries. Default is `"dashed"`.

- ...:

  Additional parameters (currently ignored).

## Value

A ggplot2 layer object that can be added to a gghic plot.

## Details

### Required Aesthetics

This geom inherits aesthetics from the Hi-C data and requires:

- `seqnames1`, `seqnames2`: Chromosome names

- `start1`, `start2`: Start positions

- `end1`, `end2`: End positions

### Gene Annotation Sources

Gene annotations can be provided in three ways:

1.  **GTF file** via `gtf_path`: Most convenient, automatically parsed
    and cached

2.  **TxDb + tx2gene**: Pre-processed annotations for custom databases

3.  Downloaded from online repositories (Ensembl, UCSC, GENCODE)

### Visualization Styles

#### Basic Style (`style = "basic"`)

Shows detailed gene structure:

- **Coding exons (CDS)**: Thick rectangles in full color

- **UTRs**: Thinner rectangles showing 5' and 3' untranslated regions

- **Introns**: Thin lines with small directional arrows indicating
  strand

- **ncRNAs**: Thin rectangles for non-coding transcripts (if
  `include_ncrna = TRUE`)

#### Arrow Style (`style = "arrow"`)

Simplified representation:

- Each gene shown as a single arrow shape

- Arrow points from 5' to 3' direction

- More compact, suitable for dense genomic regions

### Gene Placement and Layout

- Genes are automatically arranged on multiple horizontal lines to
  prevent overlap

- Genes closer than `maxgap` are placed on separate lines

- Gene symbols are centered below each gene

- Multiple tracks are vertically stacked with spacing controlled by
  `spacing_ratio`

### Multi-Chromosome Display

When displaying multiple chromosomes:

- Coordinates are automatically adjusted for proper alignment

- Vertical boundary lines separate chromosomes (if
  `draw_boundary = TRUE`)

- Gene labels maintain proper positioning across chromosome boundaries

### Gene Filtering

Use `gene_symbols` parameter to highlight specific genes:

- Provide a character vector of gene names

- Only matching genes will be displayed

- Useful for focusing on candidate genes or pathways

### Performance Considerations

- GTF parsing is cached for faster subsequent loading

- For large genomic regions, consider filtering genes with
  `gene_symbols`

- Setting `include_ncrna = FALSE` reduces the number of features

- Arrow style renders faster than basic style for dense regions

## See also

- [`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)
  for creating the base Hi-C plot

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md)
  for the Hi-C heatmap layer

- [`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md)
  and
  [`geom_track2()`](https://jasonwong-lab.github.io/gghic/reference/geom_track2.md)
  for genomic signal tracks

- [`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md)
  for chromosome ideograms

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with GTF file
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

gtf_file <- "path/to/genes.gtf"
gghic(cc) + geom_annotation(gtf_path = gtf_file)

# Filter specific genes of interest
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

# Exclude non-coding RNAs for cleaner view
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    include_ncrna = FALSE
  )

# Custom track dimensions
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    width_ratio = 1 / 40,    # Taller tracks
    spacing_ratio = 1 / 2     # More spacing
  )

# Using TxDb object directly
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tx2gene <- read.csv("tx2gene_mapping.csv")

gghic(cc) +
  geom_annotation(txdb = txdb, tx2gene = tx2gene)

# Multiple genes with custom styling
candidates <- c("PDGFRA", "KIT", "KDR", "FLT1")
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    gene_symbols = candidates,
    colour = "red",
    fill = "pink",
    fontsize = 12
  )

# Without chromosome prefix (e.g., "1" instead of "chr1")
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    chrom_prefix = FALSE
  )

# Multi-chromosome with boundaries
cc_multi <- ChromatinContacts("path/to/cooler.cool",
                              focus = c("chr4", "chr8")) |>
  import()
gghic(cc_multi) +
  geom_annotation(
    gtf_path = gtf_file,
    draw_boundary = TRUE,
    boundary_colour = "gray50",
    linetype = "dotted"
  )

# Protein-coding genes only with arrow style
gghic(cc) +
  geom_annotation(
    gtf_path = gtf_file,
    include_ncrna = FALSE,
    style = "arrow",
    fill = "#2E86AB"
  )
} # }
```
