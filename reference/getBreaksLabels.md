# Generate axis breaks and labels for multi-chromosome plots

Calculates appropriate breaks and labels for x-axis when visualizing
multiple chromosomes simultaneously. Adjusts break density based on
chromosome length.

## Usage

``` r
getBreaksLabels(data)
```

## Arguments

- data:

  A data frame or tibble with columns `seqnames`, `start`, and `end`
  representing genomic coordinates.

## Value

A list with two elements:

- `breaks`: Numeric vector of x-axis break positions

- `labels`: Character vector of formatted labels

## Details

The function:

1.  Groups data by chromosome

2.  Calculates chromosome lengths

3.  Assigns more breaks to longer chromosomes (2-10 breaks per
    chromosome)

4.  Formats labels in megabases (M)

5.  Adds newlines before chromosome transitions for clarity

Useful for multi-chromosome Hi-C visualizations where chromosomes are
displayed side-by-side.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get breaks and labels for multi-chromosome plot
df <- scaleData(cc, "balanced", log10)
breaks_labels <- getBreaksLabels(df)

ggplot(df) +
  geom_hic() +
  scale_x_continuous(
    breaks = breaks_labels$breaks,
    labels = breaks_labels$labels
  )
} # }
```
