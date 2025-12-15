# Visualize inverted Hi-C heatmap below main contact map

Creates a flipped (inverted) Hi-C contact heatmap positioned below the
main plot. Useful for comparing two datasets or showing the same data
with different color scales/transformations.

## Usage

``` r
geom_hic_under(
  mapping = NULL,
  data = NULL,
  stat = StatHicUnder,
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
  Must include `fill` aesthetic.

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

  Logical. Rasterize for performance (default: TRUE).

- dpi:

  Numeric. Rasterization resolution (default: 300).

- dev:

  Character. Graphics device for rasterization (default: `"cairo"`).

- scale:

  Numeric. Rasterization scaling factor (default: 1).

- draw_boundary:

  Logical. Draw chromosome boundaries for multi-chromosome plots
  (default: TRUE).

- boundary_colour:

  Character. Boundary line color (default: `"black"`).

- linetype:

  Boundary line type (default: `"dashed"`).

- ...:

  Additional parameters (unused).

## Value

A ggplot2 layer for the inverted Hi-C heatmap.

## Details

### Required aesthetics

- `seqnames1`, `start1`, `end1`: First anchor coordinates

- `seqnames2`, `start2`, `end2`: Second anchor coordinates

- `fill`: Color scale values

### Usage

This geom requires a main Hi-C plot created with
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md)
first. The inverted heatmap is automatically positioned below, using
independent color scales via `fill2` aesthetic.

### Common use cases

- Compare two experimental conditions

- Show raw vs. normalized data

- Display same data with different color scales

## See also

[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md),
[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`renameGeomAes()`](https://jasonwong-lab.github.io/gghic/reference/renameGeomAes.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare two datasets
cc1 <- ChromatinContacts("sample1.cool") |> import()
cc2 <- ChromatinContacts("sample2.cool") |> import()

library(ggplot2)
ggplot() +
  geom_hic(data = scaleData(cc1, "balanced", log10),
           aes(seqnames1 = seqnames1, start1 = start1, end1 = end1,
               seqnames2 = seqnames2, start2 = start2, end2 = end2,
               fill = score)) +
  geom_hic_under(data = scaleData(cc2, "balanced", log10),
                 aes(seqnames1 = seqnames1, start1 = start1, end1 = end1,
                     seqnames2 = seqnames2, start2 = start2, end2 = end2,
                     fill2 = score)) +
  scale_fill_viridis_c() +
  scale_fill2_viridis_c(option = "magma")
# Load two Hi-C datasets for comparison
cc1 <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

# Simulate a second dataset (in practice, load a different sample)
cc2 <- cc1

# Compare two Hi-C maps with different color scales
library(ggplot2)
ggplot() +
  geom_hic(
    data = scaleData(cc1, "balanced", log10),
    aes(
      seqnames1 = seqnames1, start1 = start1, end1 = end1,
      seqnames2 = seqnames2, start2 = start2, end2 = end2, fill = score
    )
  ) +
  geom_hic_under(
    data = scaleData(cc2, "balanced", log10),
    aes(
      seqnames1 = seqnames1, start1 = start1, end1 = end1,
      seqnames2 = seqnames2, start2 = start2, end2 = end2, fill2 = score
    )
  ) |>
  renameGeomAes(new_aes = c("fiil" = "fill2")) +
  scale_fill_gradientn(
    aesthetics = "fill2", colors = c("white", "blue"), name = "Sample 2"
  ) +
  theme_hic()
} # }
```
