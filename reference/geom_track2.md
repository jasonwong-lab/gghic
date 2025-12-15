# Add genomic signal tracks with direct data input (a second version)

Another version of
[`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md)
that accepts pre-loaded genomic signal data as a data frame or tibble
with explicit aesthetic mappings. Displays continuous genomic signals
(ChIP-seq, ATAC-seq, RNA-seq, etc.) as filled area tracks below Hi-C
heatmaps. Unlike
[`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md),
this function uses standard ggplot2 aesthetic mappings, making it more
flexible for custom styling and integration with ggplot2 scales and
themes.

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

  Numeric value controlling the height of each track relative to the
  Hi-C plot height. Smaller values create shorter tracks. Default is
  `1/20` (5% of Hi-C plot height).

- spacing_ratio:

  Numeric value controlling the vertical spacing between tracks as a
  proportion of track height. Default is `0.5` (50% of track height).

- data_range:

  Character string or numeric value(s) controlling y-axis scaling:

  - `"auto"`: Each track scaled independently to its own maximum value

  - `"maximum"`: All tracks share the same y-axis scale (global maximum)

  - Numeric value: Fixed maximum for all tracks (e.g., `100`)

  - Numeric vector: Individual maximum for each track (length must match
    number of tracks) Default is `"auto"`.

- rasterize:

  Logical indicating whether to rasterize the track polygons for faster
  rendering and smaller file sizes. Recommended for high-resolution
  data. Default is `TRUE`.

- dpi:

  Numeric value specifying the resolution (dots per inch) for rasterized
  tracks. Higher values increase quality but also file size. Default is
  `300`.

- dev:

  Character string specifying the graphics device for rasterization.
  Options include `"cairo"`, `"ragg"`, or other devices. Default is
  `"cairo"`.

- scale:

  Numeric scaling factor for rasterized output. Values \> 1 increase
  resolution. Default is `1`.

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

This geom requires the following aesthetics to be mapped:

- `seqnames`: Chromosome name for each genomic bin

- `start`: Start position of each bin

- `end`: End position of each bin

- `score`: Signal intensity value (numeric)

- `name`: Track identifier/label (groups bins into tracks)

### Comparison with geom_track()

|               |                            |                            |
|---------------|----------------------------|----------------------------|
| Feature       | geom_track()               | geom_track2()              |
| Input         | File paths or GRanges list | Data frame with aesthetics |
| Flexibility   | Simple, automatic          | Full ggplot2 integration   |
| Color control | `fill` parameter           | `aes(fill = ...)` + scales |
| Data loading  | Automatic from files       | Manual pre-loading         |
| Best for      | Quick visualization        | Custom styling & scales    |

### Input Data Format

Data should be a data frame or tibble where each row represents one
genomic bin:

    seqnames  start    end  score      name
    chr1      10000  10100   5.2    "ChIP-seq"
    chr1      10100  10200   6.1    "ChIP-seq"
    chr1      10200  10300   4.8    "ChIP-seq"

### Track Visualization

Each track displays:

- **Signal area**: Filled polygon showing signal intensity

- **Y-axis**: Left side with minimum (0) and maximum value labels

- **Track label**: Name (from `name` aesthetic) displayed on the left

- **Baseline**: Zero line at bottom of each track

### Color and Fill Aesthetics

This function fully integrates with ggplot2's aesthetic system:

- Map `fill` aesthetic to track names for automatic coloring

- Use
  [`scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
  [`scale_fill_brewer()`](https://ggplot2.tidyverse.org/reference/scale_brewer.html),
  etc. for custom colors

- Supports all standard ggplot2 color specifications

- Can combine with other aesthetic mappings (alpha, etc.)

### Y-Axis Scaling

Control track scaling via `data_range`:

- **Independent** (`"auto"`): Compare patterns within each track

- **Shared** (`"maximum"`): Compare absolute signal between tracks

- **Fixed** (numeric): Consistent scaling across multiple plots

### Performance

- Rasterization recommended for dense genomic data

- Adjust `dpi` to balance quality and file size

- Pre-filter data to displayed region for faster rendering

### Multi-Chromosome Support

When data includes multiple chromosomes:

- Coordinates automatically adjusted for proper alignment

- Vertical boundaries separate chromosomes (if `draw_boundary = TRUE`)

- Works seamlessly with Hi-C multi-chromosome displays

## See also

- [`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md)
  for file-based input with automatic data loading

- [`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)
  for creating the base Hi-C plot

- [`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md)
  for gene annotation tracks

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md)
  for the Hi-C heatmap layer

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

# Load and prepare track data
library(rtracklayer)
track1 <- import("path/to/H3K27ac.bw")
track1$name <- "H3K27ac"
track_df <- as.data.frame(track1)

# Basic usage with aesthetic mappings
library(ggplot2)
gghic(cc) +
  geom_track2(
    data = track_df,
    aes(
      seqnames = seqnames,
      start = start,
      end = end,
      score = score,
      name = name
    )
  )

# Multiple tracks with color mapping
track2 <- import("path/to/ATAC.bw")
track2$name <- "ATAC-seq"
tracks_df <- rbind(
  as.data.frame(track1),
  as.data.frame(track2)
)

gghic(cc) +
  geom_track2(
    data = tracks_df,
    aes(
      seqnames = seqnames,
      start = start,
      end = end,
      score = score,
      name = name,
      fill = name  # Color by track name
    )
  ) +
  scale_fill_manual(values = c("H3K27ac" = "#E63946", "ATAC-seq" = "#457B9D"))

# Using scale_fill_brewer for automatic colors
gghic(cc) +
  geom_track2(
    data = tracks_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name, fill = name)
  ) +
  scale_fill_brewer(palette = "Set2")

# Shared y-axis scaling
gghic(cc) +
  geom_track2(
    data = tracks_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name, fill = name),
    data_range = "maximum"
  )

# Custom track dimensions
gghic(cc) +
  geom_track2(
    data = track_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name),
    width_ratio = 1 / 15,    # Taller tracks
    spacing_ratio = 1         # More spacing
  )

# Fixed y-axis maximum
gghic(cc) +
  geom_track2(
    data = track_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name),
    data_range = 50  # Fix max to 50
  )

# Multiple tracks with different y-axis limits
gghic(cc) +
  geom_track2(
    data = tracks_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name),
    data_range = c(100, 50)  # Different max for each track
  )

# High-resolution rasterization
gghic(cc) +
  geom_track2(
    data = track_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name),
    rasterize = TRUE,
    dpi = 600,
    dev = "ragg"
  )

# Combine with other geoms
gghic(cc) +
  geom_track2(
    data = track_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name),
    fill = "darkblue"
  ) +
  geom_loop("path/to/loops.bedpe")

# Multi-chromosome display
cc_multi <- ChromatinContacts("path/to/cooler.cool",
                              focus = c("chr4", "chr8")) |>
  import()
gghic(cc_multi) +
  geom_track2(
    data = tracks_df,
    aes(seqnames = seqnames, start = start, end = end,
        score = score, name = name, fill = name),
    draw_boundary = TRUE
  )
} # }
```
