# Visualize TAD boundaries on Hi-C heatmap

Adds topologically associating domain (TAD) boundary annotations to Hi-C
contact maps. TAD boundaries are displayed as triangular outlines
marking the diagonal extent of each domain.

## Usage

``` r
geom_tad(
  mapping = NULL,
  data = NULL,
  stat = StatTad,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  tad_path = NULL,
  tad_gis = NULL,
  is_0_based = FALSE,
  colour = "grey",
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

- tad_path:

  Character. Path to BED-like file with TAD coordinates (columns: chrom,
  start, end). Either `tad_path` or `tad_gis` must be provided (default:
  NULL).

- tad_gis:

  GInteractions object containing TAD coordinates as interactions from
  start to end positions. Either `tad_path` or `tad_gis` must be
  provided (default: NULL).

- is_0_based:

  Logical. Whether input coordinates are 0-based (BED format). Set TRUE
  for BED files (default: FALSE).

- colour:

  Character. Color for TAD boundary lines (default: `"grey"`).

- ...:

  Additional parameters passed to layer (unused).

## Value

A ggplot2 layer that can be added to a gghic plot.

## Details

### Required aesthetics

Inherits from Hi-C data: `seqnames1`, `start1`, `end1`, `seqnames2`,
`start2`, `end2`

### TAD file format

Tab-delimited BED-like file with columns:

    chrom  start     end
    chr1   0         1000000
    chr1   1000000   2500000

### Interpretation

TAD boundaries are drawn as V-shaped lines on the Hi-C diagonal. Each
TAD is represented by two lines forming a triangle that encloses the
domain. Darker regions inside triangles represent increased intra-TAD
interactions.

## See also

[`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md)
for chromatin loops,
[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)
for creating Hi-C plots

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with TAD file
cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()
gghic(cc) + geom_tad(tad_path = "tads.bed")

# Custom color
gghic(cc) + geom_tad(tad_path = "tads.bed", colour = "blue")

# Using GInteractions object
tads <- rtracklayer::import("tads.bed") |> as("GInteractions")
gghic(cc) + geom_tad(tad_gis = tads)

# 0-based coordinates (BED format)
gghic(cc) + geom_tad(tad_path = "tads.bed", is_0_based = TRUE)

# Combined with loops
gghic(cc) +
  geom_tad(tad_path = "tads.bed", colour = "grey") +
  geom_loop(loop_path = "loops.bedpe", colour = "red")
} # }
```
