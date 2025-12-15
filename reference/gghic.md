# Create Hi-C visualization plot

High-level wrapper for publication-ready Hi-C visualizations with
genomic features.

## Usage

``` r
gghic(x, ...)

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

  ChromatinContacts, GInteractions, or data.frame with interactions.

- ...:

  Additional arguments for geom functions.

- scale_method:

  Function for data transformation (default: log10).

- ideogram:

  Logical. Add chromosome ideogram (default: FALSE).

- ideogram_width_ratio:

  Numeric. Ideogram height ratio (default: 1/30).

- ideogram_fontsize:

  Numeric. Ideogram font size (default: 10).

- ideogram_colour:

  Character. Highlight color (default: `"red"`).

- ideogram_fill:

  Character. Highlight fill (default: `"#FFE3E680"`).

- annotation:

  Logical. Add gene annotations (default: FALSE). Requires `gtf_path` in
  `...`.

- annotation_style:

  Character. `"basic"` or `"arrow"` (default: `"basic"`).

- annotation_width_ratio:

  Numeric. Track height (default: 1/50).

- annotation_spacing_ratio:

  Numeric. Gene spacing (default: 1/3).

- annotation_fontsize:

  Numeric. Label size (default: 10).

- annotation_colour:

  Character. Feature color (default: `"#48CFCB"`).

- annotation_fill:

  Character. Feature fill (default: `"#48CFCB"`).

- track:

  Logical. Add genomic tracks (default: FALSE). Requires `tracks` in
  object or `data_paths` in `...`.

- track_width_ratio:

  Numeric. Track height (default: 1/20).

- track_spacing_ratio:

  Numeric. Track spacing (default: 1/2).

- track_fill:

  Character. Track colors (default: `"black"`).

- track_fontsize:

  Numeric. Label size (default: 5).

- tad:

  Logical. Add TAD boundaries (default: FALSE). Requires `TADs` in
  object or `tad_path` in `...`.

- tad_colour:

  Character. TAD color (default: `"grey"`).

- loop:

  Logical. Add loops (default: FALSE). Requires `loops` in object or
  `loop_path` in `...`.

- loop_style:

  Character. `"circle"` or `"arc"` (default: `"circle"`).

- loop_colour:

  Character. Loop color (default: `"black"`).

- loop_fill:

  Fill color (default: NA).

- concatemer:

  Logical. Add multi-way contacts (default: FALSE).

- concatemer_width_ratio:

  Numeric. Track height (default: 1/100).

- concatemer_spacing_ratio:

  Numeric. Spacing (default: 1/5).

- expand_xaxis:

  Logical. Expand x-axis (default: FALSE).

- expand_left:

  Numeric. Left expansion in bp (default: 10×resolution).

- expand_right:

  Numeric. Right expansion in bp (default: 10×resolution).

- scale_column:

  Character string. Name of the column to use for scaling (e.g.,
  `"balanced"`, `"count"`). Used when input is `GInteractions` or
  `data.frame`.

- tad_is_0_based:

  Logical. Whether TAD coordinates are 0-based. Default is `TRUE`.

- loop_is_0_based:

  Logical. Whether loop coordinates are 0-based. Default is `TRUE`.

## Value

ggplot2 object.

## Details

Automatically handles: data scaling, coordinate transformation, feature
integration, axis labels, default theme.

Use individual `geom_*()` functions for more control:
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md),
[`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md),
[`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md),
[`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md),
[`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md),
[`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md),
[`geom_concatemer()`](https://jasonwong-lab.github.io/gghic/reference/geom_concatemer.md).

## See also

[ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md),
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md),
[`theme_hic()`](https://jasonwong-lab.github.io/gghic/reference/theme_hic.md),
[`scaleData()`](https://jasonwong-lab.github.io/gghic/reference/scaleData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
cc <- ChromatinContacts("file.cool") |> import()
gghic(cc)
gghic(cc["chr4:0-50000000"])

# With features
gghic(cc, ideogram = TRUE, genome = "hg19")
gghic(cc, tad = TRUE, tad_path = "tads.bed")
gghic(cc, loop = TRUE, loop_path = "loops.bedpe")

# With tracks
features(cc, "tracks") <- GRangesList(
  H3K27ac = rtracklayer::import("track.bw")
)
gghic(cc, track = TRUE, track_fill = "blue")
} # }
```
