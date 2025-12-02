# Scale Hi-C interaction data

Transforms and scales chromatin interaction data for visualization.
Applies a scaling function to a specified column and handles missing
values.

## Usage

``` r
scaleData(data, scale_column, scale_method, remove_na = FALSE)
```

## Arguments

- data:

  A `ChromatinContacts`, `GInteractions`, `data.frame`, or `tibble`
  object containing chromatin interaction data.

- scale_column:

  Character string. Name of the column to scale (e.g., `"balanced"`,
  `"count"`).

- scale_method:

  Function to apply for scaling. Common choices include `log10`, `log2`,
  or identity function `function(x) x`.

- remove_na:

  Logical. If `TRUE`, removes rows with `NA` or infinite values in the
  score column. Default is `FALSE`.

## Value

A tibble with columns: `seqnames1`, `start1`, `end1`, `seqnames2`,
`start2`, `end2`, and `score` (the scaled values).

## Details

The function:

1.  Converts input data to a tibble format

2.  Applies the scaling method to the specified column

3.  Creates a `score` column with the transformed values

4.  Optionally removes missing values

5.  Applies out-of-bounds squishing to ensure values are within range

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data
cc <- ChromatinContacts("path/to/cooler.cool") |>
  import()

# Scale using log10 transformation (most common)
scaled_data <- scaleData(cc, "balanced", log10)
head(scaled_data)

# Use raw counts without transformation
scaled_raw <- scaleData(cc, "count", function(x) x)

# Scale with log2 and remove missing values
scaled_clean <- scaleData(cc, "balanced", log2, remove_na = TRUE)

# From GInteractions object
gis <- interactions(cc)
scaled_gis <- scaleData(gis, "balanced", log10)

# From data frame (must have required columns)
df <- as.data.frame(gis)
scaled_df <- scaleData(df, "count", log10)

# Percentile rank transformation
percentile_rank <- function(x) {
  rank(x, na.last = "keep") / sum(!is.na(x))
}
scaled_pct <- scaleData(cc, "balanced", percentile_rank)

# Use with ggplot2 directly
library(ggplot2)
ggplot(scaled_data, aes(x = (start1 + end1) / 2, y = score)) +
  geom_point(alpha = 0.1) +
  labs(title = "Distance decay", x = "Position", y = "Log10(balanced)")
} # }
```
