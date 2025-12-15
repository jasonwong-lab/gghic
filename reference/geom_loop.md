# Visualize chromatin loops on Hi-C heatmap

Adds chromatin loop annotations to Hi-C contact maps. Loops can be
displayed as circles (points) or arcs. Automatically filters loops to
display only those within the plotted region.

## Usage

``` r
geom_loop(
  mapping = NULL,
  data = NULL,
  stat = StatLoop,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  loop_path = NULL,
  loop_gis = NULL,
  is_0_based = FALSE,
  style = "circle",
  n_arc_points = 50,
  colour = "black",
  shape = 21,
  fill = NA,
  size = grid::unit(1/80, "native"),
  stroke = 1,
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

- loop_path:

  Character. Path to loop file in BEDPE-like format with columns:
  chrom1, start1, end1, chrom2, start2, end2. Either `loop_path` or
  `loop_gis` must be provided (default: NULL).

- loop_gis:

  GInteractions object containing loop coordinates. Either `loop_path`
  or `loop_gis` must be provided (default: NULL).

- is_0_based:

  Logical. Whether input coordinates are 0-based (e.g., BED format). Set
  TRUE for BEDPE files (default: FALSE).

- style:

  Character. Visualization style:

  - `"circle"`: display loops as circular points (default)

  - `"arc"`: display loops as curved arcs connecting anchors

- n_arc_points:

  Integer. Number of points used to draw each arc when `style = "arc"`.
  Higher values produce smoother curves (default: 50).

- colour:

  Character. Color for loop markers or arcs (default: `"black"`).

- shape:

  Integer. Point shape when `style = "circle"` (default: 21 = filled
  circle).

- fill:

  Character. Fill color for points when `style = "circle"` (default: NA
  for transparent).

- size:

  Unit or numeric. Size of loop markers. Can be a grid unit or numeric
  value in mm (default: `unit(1/80, "native")`).

- stroke:

  Numeric. Line width for point borders or arc lines (default: 1).

- ...:

  Additional parameters passed to layer (unused).

## Value

A ggplot2 layer that can be added to a gghic plot.

## Details

### Required aesthetics

Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
`start2`, `end2`

### Loop file format

Tab-delimited file with columns:

    chrom1  start1  end1  chrom2  start2  end2
    chr1    1000000 1005000 chr1  2000000 2005000

### Performance tips

For large loop files, pre-filter to region of interest before plotting.
Arc style is more computationally intensive than circle style.

## See also

[`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md)
for TAD boundaries,
[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)
for creating Hi-C plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with loop file
cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()
gghic(cc) + geom_loop(loop_path = "loops.bedpe")

# Arc style visualization
gghic(cc) + 
  geom_loop(loop_path = "loops.bedpe", style = "arc", colour = "red")

# Using GInteractions object
loops <- rtracklayer::import("loops.bedpe")
gghic(cc) + geom_loop(loop_gis = loops, colour = "blue", size = 2)

# 0-based coordinates (BEDPE format)
gghic(cc) + geom_loop(loop_path = "loops.bedpe", is_0_based = TRUE)

# Customized appearance
gghic(cc) +
  geom_loop(
    loop_path = "loops.bedpe",
    colour = "darkred",
    fill = "red",
    size = 3,
    stroke = 1.5,
    shape = 21
  )
} # }
```
