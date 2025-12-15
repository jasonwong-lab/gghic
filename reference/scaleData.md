# Scale and transform Hi-C interaction data for visualization

Transforms and scales chromatin interaction data to prepare it for
visualization. Applies user-defined scaling functions (e.g., log
transformation) to interaction scores and handles missing values.

## Usage

``` r
scaleData(data, scale_column, scale_method, remove_na = FALSE)
```

## Arguments

- data:

  Input data in one of these formats:

  - ChromatinContacts object with imported interactions

  - GInteractions object with score metadata

  - data.frame or tibble with columns: seqnames1, start1, end1,
    seqnames2, start2, end2, plus score column

- scale_column:

  Character. Name of column containing values to scale. Common options:

  - `"balanced"`: ICE-normalized counts (recommended)

  - `"count"`: raw contact counts

  - Any other numeric metadata column

- scale_method:

  Function to apply for transformation. Common options:

  - `log10`: log10 transformation (default for most Hi-C data)

  - `log2`: log2 transformation

  - `log1p`: log(x + 1) transformation (handles zeros)

  - `identity` or `function(x) x`: no transformation

  - Custom function: any function that takes numeric vector and returns
    numeric vector

- remove_na:

  Logical. Whether to remove rows with NA or infinite values after
  scaling (default: FALSE). Set TRUE to remove missing data that could
  cause visualization issues.

## Value

Tibble (data frame) with standardized columns:

- `seqnames1`, `start1`, `end1`: First anchor coordinates

- `seqnames2`, `start2`, `end2`: Second anchor coordinates

- `score`: Transformed and scaled values

## Details

### Processing steps

1.  Convert input to tibble format

2.  Apply `scale_method` function to `scale_column`

3.  Create new `score` column with transformed values

4.  Optionally remove NA/infinite values

5.  Squish extreme outliers to prevent visualization artifacts

### Recommended scaling

For typical Hi-C data visualization:

- Use `"balanced"` column with `log10` transformation

- Set `remove_na = TRUE` to handle bins with no coverage

### Custom transformations

You can provide any function for scaling:

    # Square root transformation
    scaleData(cc, "count", sqrt)

    # Custom transformation
    scaleData(cc, "balanced", function(x) log2(x + 0.1))

## See also

[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md),
[`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("file.cool") |> import()

# Standard log10 scaling of balanced data
scaled_data <- scaleData(cc, "balanced", log10)

# Raw counts without transformation
scaled_raw <- scaleData(cc, "count", function(x) x)

# Log2 scaling with NA removal
scaled_clean <- scaleData(cc, "balanced", log2, remove_na = TRUE)

# Use with plotting
library(ggplot2)
ggplot() +
  geom_hic(data = scaleData(cc, "balanced", log10),
           aes(seqnames1 = seqnames1, start1 = start1, end1 = end1,
               seqnames2 = seqnames2, start2 = start2, end2 = end2,
               fill = score))
} # }
```
