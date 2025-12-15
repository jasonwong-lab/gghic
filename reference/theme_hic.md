# ggplot2 theme optimized for Hi-C contact maps

Creates a clean, publication-ready theme for Hi-C visualization with
optimized defaults: hidden y-axis, fixed coordinate ratio, and custom
color gradient. Provides fine control over axis formatting and limits.

## Usage

``` r
theme_hic(
  hide_y = TRUE,
  coord_ratio = 1,
  scale_fill_gradientn = TRUE,
  breaks = ggplot2::waiver(),
  labels = scales::unit_format(unit = "M", scale = 1e-06),
  xmin = NULL,
  xmax = NULL,
  expand_x = c(0, 0)
)
```

## Arguments

- hide_y:

  Logical. Hide y-axis elements (axis title, text, ticks) (default:
  TRUE). Y-axis typically uninformative for Hi-C plots.

- coord_ratio:

  Numeric. Aspect ratio for
  [`ggplot2::coord_fixed()`](https://ggplot2.tidyverse.org/reference/coord_fixed.html).
  1 maintains equal x/y scaling (default: 1).

- scale_fill_gradientn:

  Logical. Apply custom Hi-C color gradient (white → yellow → orange →
  red → black) (default: TRUE). Set FALSE to use custom scales.

- breaks:

  Numeric vector or
  [`ggplot2::waiver()`](https://ggplot2.tidyverse.org/reference/waiver.html).
  X-axis break positions for
  [`ggplot2::scale_x_continuous()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).
  Default uses automatic breaks.

- labels:

  Function or character vector. X-axis label formatter. Default displays
  genomic coordinates in megabases (e.g., "10M" for 10,000,000 bp).

- xmin:

  Numeric or NULL. Minimum x-axis limit. If NULL, uses data range
  (default: NULL).

- xmax:

  Numeric or NULL. Maximum x-axis limit. If NULL, uses data range
  (default: NULL).

- expand_x:

  Numeric vector (length 2). X-axis expansion in plot units:
  `c(left_expansion, right_expansion)` (default: `c(0, 0)` for no
  expansion).

## Value

List of ggplot2 theme elements and scales that can be added to a ggplot
object.

## Details

### Theme components

- Based on
  [`theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  with minimal grid lines

- Hidden y-axis (configurable)

- Fixed 1:1 coordinate ratio for accurate distance representation

- Custom color gradient optimized for Hi-C data

- Genomic coordinate formatting (Mb units)

### Color gradient

The default gradient transitions through 13 colors from white (no
contacts) to black (maximum contacts), providing excellent dynamic range
for typical Hi-C data.

### Customization

Combine with other ggplot2 theme elements for additional customization:

    gghic(cc) + theme_hic() + theme(text = element_text(size = 14))

## See also

[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md),
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("file.cool", focus = "chr4") |> import()

# Basic usage with default theme
gghic(cc) + theme_hic()

# Custom color scale
gghic(cc) +
  theme_hic(scale_fill_gradientn = FALSE) +
  scale_fill_viridis_c()

# Show y-axis
gghic(cc) + theme_hic(hide_y = FALSE)

# Custom axis limits and breaks
gghic(cc) +
  theme_hic(
    xmin = 1e6,
    xmax = 5e6,
    breaks = seq(1e6, 5e6, by = 1e6)
  )

# Custom axis labels (show kb instead of Mb)
gghic(cc) +
  theme_hic(labels = scales::unit_format(unit = "kb", scale = 1e-3))
} # }
```
