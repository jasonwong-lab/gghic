# geom_loop

A ggplot2 geom for drawing chromatin loops on the heatmap.

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

  A path to the loop file. Default is `NULL`.

- loop_gis:

  An InteractionSet object of loops. Default is `NULL`.

- is_0_based:

  Whether the loop file is 0-based or not. Default is `FALSE`.

- style:

  The style for drawing loops: `"circle"` for points/circles or `"arc"`
  for arcs under the Hi-C heatmap. Default is `"circle"`.

- n_arc_points:

  Number of points used to draw each arc (only used when
  `style = "arc"`). Default is `50`.

- colour:

  The color of the loops. Default is `"black"`.

- shape:

  The shape of the loops (only used when `style = "circle"`). Default is
  `21`.

- fill:

  The fill color of the loops (only used when `style = "circle"`).
  Default is `NA`.

- size:

  The size of the loops (only used when `style = "circle"`). Default is
  `grid::unit(1 / 80, "native")`.

- stroke:

  The line width of the loops. Default is `1`.

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

# Add loops from file
loop_file <- "path/to/loops.bedpe"
gghic(cc) + geom_loop(loop_path = loop_file)

# Draw loops as arcs
gghic(cc) + geom_loop(loop_path = loop_file, style = "arc")
} # }
```
