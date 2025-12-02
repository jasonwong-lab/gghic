# geom_hic

Creates a Hi-C contact heatmap layer using diamond/rotated square
geometry. This is the core visualization layer for displaying chromatin
interaction frequencies.

## Usage

``` r
geom_hic(
  mapping = NULL,
  data = NULL,
  stat = StatHic,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  rasterize = TRUE,
  dpi = 300,
  dev = "cairo",
  scale = 1,
  draw_boundary = TRUE,
  boundary_colour = "black",
  linetype = "dashed",
  ...
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).

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

- rasterize:

  Logical. Rasterize the heatmap for better performance with large
  datasets. Default is `TRUE`.

- dpi:

  Numeric. Resolution for rasterization. Default is `300`.

- dev:

  Character. Graphics device for rasterization. Default is `"cairo"`.

- scale:

  Numeric. Scaling factor for rasterization. Default is `1`.

- draw_boundary:

  Logical. Draw boundary lines between chromosomes in multi-chromosome
  plots. Default is `TRUE`.

- boundary_colour:

  Character. Color for boundary lines. Default is `"black"`.

- linetype:

  Line type for boundaries. Default is `"dashed"`.

- ...:

  Additional parameters (currently unused).

## Value

A ggplot2 layer.

## Details

`geom_hic()` transforms rectangular interaction data into diamond-shaped
heatmap tiles. Each interaction is plotted as a rotated square at 45
degrees, creating the characteristic triangular Hi-C heatmap
visualization.

**Required aesthetics:**

- `seqnames1` - Chromosome name for first anchor

- `start1` - Start position for first anchor

- `end1` - End position for first anchor

- `seqnames2` - Chromosome name for second anchor

- `start2` - Start position for second anchor

- `end2` - End position for second anchor

- `fill` - Value to map to color (typically scaled contact frequency)

## See also

[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`theme_hic()`](https://jasonwong-lab.github.io/gghic/reference/theme_hic.md),
[`scaleData()`](https://jasonwong-lab.github.io/gghic/reference/scaleData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load and import Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool") |>
  import()

# Basic Hi-C heatmap using geom_hic
library(ggplot2)
ggplot() +
  geom_hic(
    data = scaleData(cc["chr4"], "balanced", log10), aes(
      seqnames1 = seqnames1, start1 = start1, end1 = end1,
      seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
    )
  ) +
  scale_fill_gradientn(colors = c("white", "red")) +
  theme_hic()

# Multi-chromosome plot
cc_multi <- cc[c("chr4", "chr11")]
ggplot() +
  geom_hic(
    data = scaleData(cc_multi, "balanced", log10), aes(
      seqnames1 = seqnames1, start1 = start1, end1 = end1,
      seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
    ), draw_boundary = TRUE, boundary_colour = "blue"
  ) +
  theme_hic()
} # }
```
