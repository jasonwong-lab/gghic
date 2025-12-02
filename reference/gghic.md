# Create Hi-C visualization plot

High-level wrapper function to create publication-ready Hi-C contact map
visualizations with optional genomic features. Automatically handles
data transformation, feature integration, and theme application.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
gghic(
  x,
  scale_method = log10,
  ideogram = FALSE,
  ideogram_width_ratio = 1/30,
  ideogram_fontsize = 10,
  ideogram_colour = "red",
  ideogram_fill = "#FFE3E680",
  annotation = FALSE,
  annotation_style = "basic",
  annotation_width_ratio = 1/50,
  annotation_spacing_ratio = 1/3,
  annotation_fontsize = 10,
  annotation_colour = "#48CFCB",
  annotation_fill = "#48CFCB",
  track = FALSE,
  track_width_ratio = 1/20,
  track_spacing_ratio = 1/2,
  track_fill = "black",
  track_fontsize = 5,
  tad = FALSE,
  tad_colour = "grey",
  loop = FALSE,
  loop_style = "circle",
  loop_colour = "black",
  loop_fill = NA,
  concatemer = FALSE,
  concatemer_width_ratio = 1/100,
  concatemer_spacing_ratio = 1/5,
  expand_xaxis = FALSE,
  expand_left = NULL,
  expand_right = NULL,
  ...
)

# S4 method for class 'data.frame'
gghic(
  x,
  scale_column = "balanced",
  scale_method = log10,
  ideogram = FALSE,
  ideogram_width_ratio = 1/30,
  ideogram_fontsize = 10,
  ideogram_colour = "red",
  ideogram_fill = "#FFE3E680",
  annotation = FALSE,
  annotation_style = "basic",
  annotation_width_ratio = 1/50,
  annotation_spacing_ratio = 1/3,
  annotation_fontsize = 10,
  annotation_colour = "#48CFCB",
  annotation_fill = "#48CFCB",
  track = FALSE,
  track_width_ratio = 1/20,
  track_spacing_ratio = 1/2,
  track_fill = "black",
  track_fontsize = 5,
  tad = FALSE,
  tad_is_0_based = FALSE,
  tad_colour = "grey",
  loop = FALSE,
  loop_style = "circle",
  loop_is_0_based = FALSE,
  loop_colour = "black",
  loop_fill = NA,
  concatemer = FALSE,
  concatemer_width_ratio = 1/100,
  concatemer_spacing_ratio = 1/5,
  expand_xaxis = FALSE,
  expand_left = NULL,
  expand_right = NULL,
  ...
)

# S4 method for class 'GInteractions'
gghic(
  x,
  scale_column = "balanced",
  scale_method = log10,
  ideogram = FALSE,
  ideogram_width_ratio = 1/30,
  ideogram_fontsize = 10,
  ideogram_colour = "red",
  ideogram_fill = "#FFE3E680",
  annotation = FALSE,
  annotation_style = "basic",
  annotation_width_ratio = 1/50,
  annotation_spacing_ratio = 1/3,
  annotation_fontsize = 10,
  annotation_colour = "#48CFCB",
  annotation_fill = "#48CFCB",
  track = FALSE,
  track_width_ratio = 1/20,
  track_spacing_ratio = 1/2,
  track_fill = "black",
  track_fontsize = 5,
  tad = FALSE,
  tad_is_0_based = FALSE,
  tad_colour = "grey",
  loop = FALSE,
  loop_style = "circle",
  loop_is_0_based = FALSE,
  loop_colour = "black",
  loop_fill = NA,
  concatemer = FALSE,
  concatemer_width_ratio = 1/100,
  concatemer_spacing_ratio = 1/5,
  expand_xaxis = FALSE,
  expand_left = NULL,
  expand_right = NULL,
  ...
)
```

## Arguments

- x:

  A `ChromatinContacts` object with imported interaction data, or a
  `GInteractions` object, or a `data.frame`/`tibble` with interaction
  data.

- scale_method:

  Function to apply for data transformation. Common choices: `log10`
  (default), `log2`, `function(x) x` (no transformation).

- ideogram:

  Logical. Add chromosome ideogram track. Default is `FALSE`.

- ideogram_width_ratio:

  Numeric. Height of ideogram relative to heatmap height. Default is
  `1/30`.

- ideogram_fontsize:

  Numeric. Font size for ideogram labels. Default is `10`.

- ideogram_colour:

  Character. Color for highlighted region on ideogram. Default is
  `"red"`.

- ideogram_fill:

  Character. Fill color for highlighted region. Default is `"#FFE3E680"`
  (transparent red).

- annotation:

  Logical. Add gene annotation track. Requires `gtf_path` in `...`.
  Default is `FALSE`.

- annotation_style:

  Character. Style for gene annotation: `"basic"` or `"arrow"`. Default
  is `"basic"`.

- annotation_width_ratio:

  Numeric. Height of annotation track relative to heatmap. Default is
  `1/50`.

- annotation_spacing_ratio:

  Numeric. Spacing between genes. Default is `1/3`.

- annotation_fontsize:

  Numeric. Font size for gene labels. Default is `10`.

- annotation_colour:

  Character. Color for gene features. Default is `"#48CFCB"`.

- annotation_fill:

  Character. Fill color for gene features. Default is `"#48CFCB"`.

- track:

  Logical. Add genomic signal tracks (e.g., ChIP-seq). Requires `tracks`
  in the `ChromatinContacts` object or `data_paths` in `...`. Default is
  `FALSE`.

- track_width_ratio:

  Numeric. Height of track area relative to heatmap. Default is `1/20`.

- track_spacing_ratio:

  Numeric. Spacing between multiple tracks. Default is `1/2`.

- track_fill:

  Character or vector. Colors for track signals. Default is `"black"`.

- track_fontsize:

  Numeric. Font size for track labels. Default is `5`.

- tad:

  Logical. Add TAD (topologically associating domain) boundaries.
  Requires `TADs` in the `ChromatinContacts` object or `tad_path` in
  `...`. Default is `FALSE`.

- tad_colour:

  Character. Color for TAD boundaries. Default is `"grey"`.

- loop:

  Logical. Add chromatin loop arcs. Requires `loops` in the
  `ChromatinContacts` object or `loop_path` in `...`. Default is
  `FALSE`.

- loop_style:

  Character. Style for loops: `"circle"` or `"arc"`. Default is
  `"circle"`.

- loop_colour:

  Character. Color for loop arcs. Default is `"black"`.

- loop_fill:

  Fill color for loop arcs. Default is `NA`.

- concatemer:

  Logical. Add concatemer visualization for multi-way contacts. Requires
  `multi_contacts` in object. Default is `FALSE`.

- concatemer_width_ratio:

  Numeric. Height of concatemer track. Default is `1/100`.

- concatemer_spacing_ratio:

  Numeric. Spacing between concatemers. Default is `1/5`.

- expand_xaxis:

  Logical. Expand x-axis to show full chromosome context. Default is
  `FALSE`.

- expand_left:

  Numeric. Left expansion in base pairs. If `NULL` (default), uses 10×
  resolution.

- expand_right:

  Numeric. Right expansion in base pairs. If `NULL` (default), uses 10×
  resolution.

- ...:

  Additional arguments passed to individual geom functions:

  - For
    [`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md):
    `genome`, `highlight`, `length_ratio`

  - For
    [`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md):
    `gtf_path`, `style`, `maxgap`, `gene_symbols`, `include_ncrna`

  - For
    [`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md):
    `data_paths`, `data_range`, `rasterize`

  - For
    [`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md):
    `tad_path` (path to BED file), `tad_is_0_based` (logical, default
    TRUE), `stroke`

  - For
    [`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md):
    `loop_path` (path to BEDPE file), `loop_is_0_based` (logical,
    default TRUE), `stroke`

  - For
    [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md):
    `draw_boundary`, `rasterize`

- scale_column:

  Character string. Name of the column to use for scaling (e.g.,
  `"balanced"`, `"count"`). Used when input is `GInteractions` or
  `data.frame`.

- tad_is_0_based:

  Logical. Whether TAD coordinates are 0-based. Default is `TRUE`.

- loop_is_0_based:

  Logical. Whether loop coordinates are 0-based. Default is `TRUE`.

## Value

A `ggplot2` object that can be further customized with additional
ggplot2 layers and functions.

## Details

`gghic()` provides a high-level interface for creating Hi-C
visualizations. It automatically:

- Scales interaction data using the specified method

- Applies appropriate coordinate transformations

- Integrates genomic features from the `ChromatinContacts` object

- Sets up axis labels and breaks

- Applies a clean default theme

**For more control**, use individual `geom_*()` functions with
`ggplot2`:

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md) -
  Base heatmap layer

- [`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md) -
  Chromosome ideogram

- [`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md) -
  Gene annotations

- [`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md) -
  Genomic signal tracks

- [`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md) -
  TAD boundaries

- [`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md) -
  Chromatin loops

- [`geom_concatemer()`](https://jasonwong-lab.github.io/gghic/reference/geom_concatemer.md) -
  Multi-way contacts

## See also

- [ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md) -
  Main data class

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md) -
  Base heatmap layer

- [`theme_hic()`](https://jasonwong-lab.github.io/gghic/reference/theme_hic.md) -
  Default theme

- [`scaleData()`](https://jasonwong-lab.github.io/gghic/reference/scaleData.md) -
  Data transformation

## Examples

``` r
if (FALSE) { # \dontrun{
# === Basic Usage ===

# Load and visualize Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool") |>
  import()
gghic(cc)

# Focus on specific region
cc["chr4:0-50000000"] |>
  gghic()

# === With Ideogram ===

# Add chromosome ideogram for context
gghic(cc, ideogram = TRUE, genome = "hg19")

# Customize ideogram appearance
gghic(cc,
  ideogram = TRUE,
  ideogram_width_ratio = 1/25,
  ideogram_colour = "blue",
  ideogram_fill = "#ADD8E680"
)

# === Adding Genomic Features ===

# Add TAD boundaries
tad_file <- "path/to/tads.bed"
gghic(cc, tad = TRUE, tad_path = tad_file, tad_colour = "darkgreen")

# Add chromatin loops
loop_file <- "path/to/loops.bedpe"
gghic(cc, loop = TRUE, loop_path = loop_file, loop_colour = "red")

# Add signal tracks (e.g., ChIP-seq)
track1 <- "path/to/track1.bw"
track2 <- "path/to/track2.bw"
tracks <- GRangesList(
  H3K27ac_rep1 = rtracklayer::import(track1),
  H3K27ac_rep2 = rtracklayer::import(track2)
)
features(cc, "tracks") <- tracks

gghic(cc, track = TRUE, track_fill = c("blue", "red"))

# === Axis Expansion ===

# Expand x-axis to show genomic context
gghic(
  cc,
  expand_xaxis = TRUE,
  expand_left = 1e6,   # Extend 1 Mb to the left
  expand_right = 1e6   # Extend 1 Mb to the right
)
} # }
```
