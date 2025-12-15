# Add genomic signal tracks to Hi-C plots

Displays continuous genomic signal data (ChIP-seq, ATAC-seq, RNA-seq,
etc.) as tracks below Hi-C heatmaps. Data is read from standard genomic
file formats (BigWig, BED, bedGraph) or provided as GRanges objects.
Each track shows signal intensity as a filled area plot with automatic
y-axis scaling and labeled axes. Multiple tracks can be stacked
vertically with independent or shared scaling.

## Usage

``` r
geom_track(
  mapping = NULL,
  data = NULL,
  stat = StatTrack,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  data_paths = NULL,
  track_grs = NULL,
  width_ratio = 1/20,
  spacing_ratio = 0.5,
  data_range = c("auto", "maximum"),
  fill = "black",
  fontsize = 5,
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

- data_paths:

  Character vector of file paths to genomic signal files. Supported
  formats include BigWig (.bw, .bigWig), BED, and bedGraph files. Files
  are automatically read using rtracklayer. Named vectors are
  recommended to provide track labels (e.g.,
  `c("H3K27ac" = "chip.bw", "ATAC" = "atac.bw")`). Either `data_paths`
  or `track_grs` must be provided. Default is `NULL`.

- track_grs:

  Named list of GRanges objects containing genomic signal data. Each
  GRanges should have a `score` metadata column with signal values.
  Names are used as track labels. Either `data_paths` or `track_grs`
  must be provided. Default is `NULL`.

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

- fill:

  Character string or named character vector specifying fill colors for
  tracks. If named, names should match track names. If unnamed, colors
  are recycled across tracks. Default is `"black"`.

- fontsize:

  Numeric value specifying the font size for track labels and axis
  annotations. Default is `5`.

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

This geom inherits aesthetics from the Hi-C data and requires:

- `seqnames1`, `seqnames2`: Chromosome names

- `start1`, `start2`: Start positions

- `end1`, `end2`: End positions

### Input Data Formats

Tracks can be provided via two methods:

#### File-based input (`data_paths`)

- Reads directly from genomic signal files

- Automatically extracts data for displayed genomic region

- Supports BigWig (.bw), BED, bedGraph formats

- Named vector provides track labels

#### GRanges input (`track_grs`)

- Pre-loaded genomic data as GRanges objects

- Must contain `score` metadata column

- Useful for custom data or pre-processed signals

- Automatically filtered to displayed region

### Track Visualization

Each track displays:

- **Signal area**: Filled polygon showing signal intensity

- **Y-axis**: Left side with minimum (0) and maximum value labels

- **Track label**: Name displayed on the left

- **Baseline**: Zero line at bottom of each track

### Y-Axis Scaling Options

The `data_range` parameter controls how tracks are scaled:

- **Independent scaling** (`"auto"`): Each track scaled to its own
  maximum, best for comparing patterns within each dataset

- **Shared scaling** (`"maximum"`): All tracks use the same scale,
  enabling direct amplitude comparison between tracks

- **Fixed scaling** (numeric): Manually set maximum value(s) for
  consistent scaling across multiple plots or known value ranges

### Color Specification

The `fill` parameter accepts:

- Single color: Applied to all tracks (e.g., `"black"`)

- Named vector: Specific colors for each track (e.g.,
  `c("ChIP" = "blue", "ATAC" = "red")`)

- Unnamed vector: Colors recycled across tracks in order

### Performance Optimization

- **Rasterization** (`rasterize = TRUE`): Converts high-resolution
  tracks to raster images, dramatically reducing file size and rendering
  time

- **DPI setting**: Balance quality vs. file size (150-300 dpi typical)

- **Data filtering**: Only overlapping regions are loaded from files

### Multi-Chromosome Display

When displaying multiple chromosomes:

- Track coordinates automatically adjusted for proper alignment

- Vertical boundaries drawn between chromosomes (if
  `draw_boundary = TRUE`)

- Y-axis labels positioned consistently across chromosomes

## See also

- [`geom_track2()`](https://jasonwong-lab.github.io/gghic/reference/geom_track2.md)
  for a simplified version with direct data frame input

- [`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)
  for creating the base Hi-C plot

- [`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md)
  for gene annotation tracks

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md)
  for the Hi-C heatmap layer

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with BigWig files
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

track1 <- "path/to/H3K27ac.bw"
track2 <- "path/to/ATAC.bw"

# Add tracks using file paths
gghic(cc) +
  geom_track(data_paths = c(track1, track2))

# Named tracks with labels
gghic(cc) +
  geom_track(
    data_paths = c(
      "H3K27ac" = "path/to/H3K27ac.bw",
      "ATAC-seq" = "path/to/ATAC.bw",
      "RNA-seq" = "path/to/RNA.bw"
    )
  )

# Custom colors for each track
gghic(cc) +
  geom_track(
    data_paths = c("H3K27ac" = track1, "ATAC" = track2),
    fill = c("H3K27ac" = "#E63946", "ATAC" = "#457B9D")
  )

# Shared y-axis scaling for comparison
gghic(cc) +
  geom_track(
    data_paths = c("Condition1" = "cond1.bw", "Condition2" = "cond2.bw"),
    data_range = "maximum"
  )

# Fixed y-axis maximum
gghic(cc) +
  geom_track(
    data_paths = c("ChIP-seq" = track1),
    data_range = 100  # Fix max to 100
  )

# Different maximums for each track
gghic(cc) +
  geom_track(
    data_paths = c("Track1" = track1, "Track2" = track2),
    data_range = c(50, 100)  # Track1 max=50, Track2 max=100
  )

# Using GRanges objects
library(rtracklayer)
gr1 <- import("path/to/track1.bw")
gr2 <- import("path/to/track2.bw")

gghic(cc) +
  geom_track(
    track_grs = list("H3K27ac" = gr1, "H3K4me3" = gr2)
  )

# Custom track dimensions
gghic(cc) +
  geom_track(
    data_paths = c("Signal" = track1),
    width_ratio = 1 / 15,   # Taller tracks
    spacing_ratio = 1        # More spacing
  )

# Disable rasterization for vector output
gghic(cc) +
  geom_track(
    data_paths = c("ChIP" = track1),
    rasterize = FALSE
  )

# High-resolution rasterization
gghic(cc) +
  geom_track(
    data_paths = c("Signal" = track1),
    rasterize = TRUE,
    dpi = 600,
    dev = "ragg"
  )

# Multi-chromosome with custom boundaries
cc_multi <- ChromatinContacts("path/to/cooler.cool",
                              focus = c("chr4", "chr8")) |>
  import()
gghic(cc_multi) +
  geom_track(
    data_paths = c("ATAC" = track1),
    draw_boundary = TRUE,
    boundary_colour = "gray50",
    linetype = "dotted"
  )
} # }
```
