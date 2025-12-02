# Rename aesthetic mappings in a ggplot2 layer

Advanced utility function for remapping aesthetic names in a ggplot2
layer. Used internally by gghic for creating complex multi-layer
visualizations with custom aesthetics.

## Usage

``` r
renameGeomAes(layer, new_aes, mapping = NULL)
```

## Arguments

- layer:

  A ggplot2 layer object.

- new_aes:

  Named character vector. Maps new aesthetic names to original names
  (e.g., `c("fill2" = "fill")`).

- mapping:

  Optional aesthetic mapping. If NULL, created from `new_aes`.

## Value

A modified ggplot2 layer with remapped aesthetics.

## Details

This function creates a new ggproto object that remaps aesthetics,
allowing multiple layers to use different aesthetic names (e.g., `fill`
and `fill2`) with independent scales.

This is primarily an internal function used by
[`geom_hic_under()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic_under.md)
to allow overlaying two heatmaps with independent color scales.

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically used internally
layer <- geom_tile(aes(fill = score))
renamed_layer <- renameGeomAes(layer, c("fill2" = "fill"))
} # }
```
