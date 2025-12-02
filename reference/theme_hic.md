# theme_hic

Generate a ggplot2 theme for Hi-C plots.

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

  A logical value indicating whether to hide the y-axis. Default is
  `TRUE`.

- coord_ratio:

  The ratio of the
  [ggplot2::coord_fixed](https://ggplot2.tidyverse.org/reference/coord_fixed.html).
  Default is `1`.

- scale_fill_gradientn:

  A logical value indicating whether to apply the gradient color scale.
  Default is `TRUE`.

- breaks:

  This argument is passed to ggplot2::scale_x_continuous.

- labels:

  This argument is passed to ggplot2::scale_x_continuous.

- xmin:

  The minimum x-axis value. Default is `NULL`.

- xmax:

  The maximum x-axis value. Default is `NULL`.

- expand_x:

  A numeric vector of length 2 to expand the x-axis. Default is
  `c(0, 0)`.

## Value

A ggplot2 theme.

## Details

If either `xmin` or `xmax` is `NULL`, the x-axis will not be expanded.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool", focus = "chr4") |>
  import()

# Apply default theme_hic
gghic(cc) + theme_hic()

# Disable default color scale
gghic(cc) +
  theme_hic(scale_fill_gradientn = FALSE) +
  scale_fill_viridis_c()
} # }
```
