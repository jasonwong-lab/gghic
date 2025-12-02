# geom_concatemer

A ggplot2 geom for concatemers containing multi-way contacts.

## Usage

``` r
geom_concatemer(
  mapping = NULL,
  data = NULL,
  stat = StatConcatemer,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  width_ratio = 1/100,
  spacing_ratio = 1/5,
  concatemer_grs = NULL,
  concatemer_path = NULL,
  group_identifier = NULL,
  fill = "black",
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

- width_ratio:

  The ratio of the width of each concatemer track relative to the height
  of the Hi-C plot. Default is `1/100`.

- spacing_ratio:

  The ratio of the spacing between two tracks. Default is `1/5`.

- concatemer_grs:

  The GRanges object of the concatemer tracks. Default is `NULL`.

- concatemer_path:

  The paths to the concatemer files. Default is `NULL`.

- group_identifier:

  A character indicating the column name in the concatemer GRanges
  object that identifies the group of concatemers. Default is `NULL`.

- fill:

  The fill color of the gene model track. Default is `"black"`.

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

# Load concatemer data (multi-way contacts from Pore-C)
concatemers <- rtracklayer::import("path/to/concatemers.bed")

# Add concatemer tracks to Hi-C plot
gghic(cc) +
  geom_concatemer(
    concatemer_grs = concatemers,
    group_identifier = "read_name",
    fill = "blue"
  )
} # }
```
