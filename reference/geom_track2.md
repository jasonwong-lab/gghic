# geom_track2

A second version of
[`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md).

## Usage

``` r
geom_track2(
  mapping = NULL,
  data = NULL,
  stat = StatTrack2,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  width_ratio = 1/20,
  spacing_ratio = 0.5,
  data_range = c("auto", "maximum"),
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

- width_ratio:

  The ratio of the width of each track relative to the height of the
  Hi-C plot. Default is `1/20`.

- spacing_ratio:

  The ratio of the spacing between two tracks. Default is `0.5`.

- data_range:

  The range of the x axis. It can be `"auto"`, `"maximum"`, or a number
  (vector). Default is `"auto"`.

- rasterize:

  Whether to rasterize the plot or not. Default is `TRUE`.

- dpi:

  The resolution of the rasterised plot. Default is `300`.

- dev:

  The device to rasterise the plot. Default is `"cairo"`.

- scale:

  The scale of the rasterised plot. Default is `1`.

- draw_boundary:

  Whether to draw the boundary line or not when plotting multiple
  chromosomes. Default is `TRUE`.

- boundary_colour:

  The color of the boundary line. Default is `"black"`.

- linetype:

  The line type of the boundary line. Default is `"dashed"`.

- ...:

  Parameters to be ignored.

## Value

A ggplot object.

## Details

Requires the following aesthetics:

- seqnames

- start

- end

- score

- name

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

# Load track data from BigWig
track1 <- rtracklayer::import("path/to/track1.bw")
track1$name <- "ChIP-seq"

# Add track using data frame with aesthetics
library(ggplot2)
gghic(cc) +
  geom_track2(
    data = as.data.frame(track1),
    aes(
      seqnames = seqnames,
      start = start,
      end = end,
      score = score,
      name = name
    )
  )

# Multiple tracks with different colors
track2 <- track1
track2$name <- "ATAC-seq"
tracks_df <- rbind(
  as.data.frame(track1),
  as.data.frame(track2)
)

gghic(cc) +
  geom_track2(
    data = tracks_df, aes(
      seqnames = seqnames, start = start, end = end, score = score,
      name = name, fill = name
    ), width_ratio = 1 / 25
  ) +
  scale_fill_manual(values = c("ChIP-seq" = "blue", "ATAC-seq" = "red"))
} # }
```
