# Convert Hypergraph to Long Format for Visualization

Convert Hypergraph to Long Format for Visualization

## Usage

``` r
tidyHypergraph(hg, max_reads = NULL)
```

## Arguments

- hg:

  A hypergraph object from
  [`buildHypergraph()`](https://jasonwong-lab.github.io/gghic/reference/buildHypergraph.md).

- max_reads:

  Integer maximum number of reads to include. Default: NULL (all).

## Value

A tibble with columns: `read_name`, `read_idx`, `bin_id`, `chrom`,
`bin`, `bin_idx`, `n_contacts`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Pore-C pairs data
pairs_df <- readr::read_tsv("path/to/porec_pairs.txt")

# Build hypergraph from Pore-C data
hg <- buildHypergraph(
  pairs = pairs_df,
  bin_size = 100000,
  chrom = "chr22",
  quantile = 0.85,
  min_multiway = 3
)

# Convert to tidy format
tidy_hg <- tidyHypergraph(hg)

# Limit to top 50 reads by contact count
tidy_hg_top50 <- tidyHypergraph(hg, max_reads = 50)

# Use with ggplot2
library(ggplot2)
ggplot(tidy_hg_top50, aes(x = read_idx, y = bin_idx, group = read_idx)) +
  geom_line(aes(color = n_contacts)) +
  geom_point() +
  theme_minimal()
} # }
```
