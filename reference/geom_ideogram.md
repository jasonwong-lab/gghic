# Add chromosome ideogram annotation to Hi-C plot

Displays chromosome ideograms (cytogenetic band representations) above
the Hi-C contact map, showing the genomic context with Giemsa staining
patterns and highlighting the displayed region. Ideograms provide visual
reference for chromosome structure, centromeres, and cytogenetic
landmarks.

## Usage

``` r
geom_ideogram(
  mapping = NULL,
  data = NULL,
  stat = StatIdeogram,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  genome = "hg19",
  chrom_prefix = TRUE,
  highlight = TRUE,
  show_coord = FALSE,
  width_ratio = 1/30,
  length_ratio = 0.8,
  fontsize = 10,
  colour = "red",
  fill = "#FFE3E680",
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
    [`stat_count()`](https://ggplot2.tidyverse.org/reference/geom_bar.html),
    give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- position:

  A position adjustment to use on the data for this layer. This can be
  used in various ways, including to prevent overplotting and improving
  the display. The `position` argument accepts the following:

  - The result of calling a position function, such as
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html).
    This method allows for passing extra arguments to the position.

  - A string naming the position adjustment. To give the position as a
    string, strip the function name of the `position_` prefix. For
    example, to use
    [`position_jitter()`](https://ggplot2.tidyverse.org/reference/position_jitter.html),
    give the position as `"jitter"`.

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

- genome:

  Character. Genome assembly version for retrieving cytogenetic band
  information from UCSC. Supported genomes include:

  - `"hg38"`: Human (GRCh38/hg38)

  - `"hg19"`: Human (GRCh37/hg19) (default)

  - `"mm10"`: Mouse (GRCm38)

  - `"mm39"`: Mouse (GRCm39)

  - Other UCSC genome assemblies with cytoBand tables

- chrom_prefix:

  Logical. Whether chromosome names in the data include "chr" prefix
  (e.g., "chr1" vs "1"). Set FALSE for Ensembl-style naming (default:
  TRUE).

- highlight:

  Logical. Draw a colored outline around the displayed genomic region on
  the ideogram to emphasize the Hi-C map extent (default: TRUE).

- show_coord:

  Logical. Display genomic coordinates (start-end) next to chromosome
  names on the ideogram (default: FALSE). Useful for showing exact
  region boundaries.

- width_ratio:

  Numeric. Height of ideogram relative to Hi-C plot height (default:
  1/30). Larger values create taller ideograms.

- length_ratio:

  Numeric. Fraction of Hi-C plot width used for ideogram length
  (default: 0.8 = 80%). Controls horizontal scaling to leave space for
  labels.

- fontsize:

  Numeric. Font size in points for chromosome labels and coordinates
  (default: 10).

- colour:

  Character. Border color for the highlighted region box (default:
  `"red"`).

- fill:

  Character. Fill color for the highlighted region box, using RGBA hex
  format for transparency (default: `"#FFE3E680"` = semi-transparent
  light red).

- ...:

  Additional parameters (unused).

## Value

A ggplot2 layer that can be added to a Hi-C plot.

## Details

### Required aesthetics

Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
`start2`, `end2`

### Ideogram structure

Chromosomes are displayed horizontally above the Hi-C map with:

- **Cytogenetic bands**: Giemsa staining patterns (G-bands) shown with
  standard colors (light/dark representing staining intensity)

- **Centromeres**: Typically appear as darker bands near the middle

- **Highlighted region**: Colored box showing the exact genomic region
  displayed in the Hi-C map below

- **Labels**: Chromosome names or coordinates on the right side

### Multi-chromosome display

When visualizing multiple chromosomes, ideograms are stacked vertically
in the same order as they appear in the Hi-C plot.

### Genome assembly selection

Ensure the `genome` parameter matches your data's assembly. Mismatched
assemblies will show incorrect cytogenetic band patterns or fail to
retrieve band data.

### Performance considerations

Ideogram data is fetched from UCSC Genome Browser on first use per
session. Subsequent calls use cached data for improved performance.

## See also

[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md),
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with human genome (hg19)
cc <- ChromatinContacts("sample.cool", focus = "chr4") |> import()
gghic(cc) + geom_ideogram(genome = "hg19")

# Use human GRCh38/hg38 assembly
gghic(cc) + geom_ideogram(genome = "hg38")

# Mouse genome
cc_mouse <- ChromatinContacts("mouse.cool", focus = "chr1") |> import()
gghic(cc_mouse) + geom_ideogram(genome = "mm10")

# Show coordinates instead of just chromosome names
gghic(cc) +
  geom_ideogram(genome = "hg19", show_coord = TRUE)

# Customize highlight colors
gghic(cc) +
  geom_ideogram(
    genome = "hg19",
    highlight = TRUE,
    colour = "blue",
    fill = "#ADD8E680"  # Semi-transparent light blue
  )

# Adjust ideogram size
gghic(cc) +
  geom_ideogram(
    genome = "hg19",
    width_ratio = 1/20,    # Taller ideogram
    length_ratio = 0.9,    # Wider ideogram
    fontsize = 12          # Larger labels
  )

# Multiple chromosomes with ideograms
cc_multi <- ChromatinContacts("sample.cool", focus = "chr1|chr2") |>
  import()
gghic(cc_multi) +
  geom_ideogram(genome = "hg19", show_coord = TRUE)

# Ensembl-style chromosome names (without "chr" prefix)
cc_ensembl <- ChromatinContacts("ensembl.cool", focus = "1") |> import()
gghic(cc_ensembl) +
  geom_ideogram(genome = "hg19", chrom_prefix = FALSE)

# Disable highlighting for cleaner look
gghic(cc) +
  geom_ideogram(genome = "hg19", highlight = FALSE)

# Complete publication-ready plot
gghic(cc, ideogram = FALSE) +  # Disable auto-ideogram
  geom_ideogram(
    genome = "hg19",
    highlight = TRUE,
    colour = "darkred",
    fill = "#FFE3E650",
    fontsize = 11
  ) +
  theme_hic()
} # }
```
