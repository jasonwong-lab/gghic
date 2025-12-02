# Visualize Multi-way Contacts as Hypergraph

Creates a hypergraph visualization for multi-way chromatin contacts
where genomic bins are displayed as points (nodes) and reads are shown
as vertical lines (hyperedges) connecting multiple bins.

## Usage

``` r
geom_hypergraph(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  line_width = 0.3,
  line_alpha = 0.6,
  point_size = 2,
  point_alpha = 0.8,
  colour_by = "n_contacts",
  palette = "viridis",
  ...
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). Only `x`,
  `y`, and `group` are used.

- data:

  A data frame from `hypergraph_to_tidy()`.

- stat:

  The statistical transformation to use on the data. Default:
  "identity".

- position:

  Position adjustment. Default: "identity".

- na.rm:

  If `FALSE` (default), missing values are removed with a warning.

- show.legend:

  Logical. Should this layer be included in legends?

- inherit.aes:

  If `FALSE`, overrides the default aesthetics.

- line_width:

  Numeric line width for hyperedges. Default: 0.3.

- line_alpha:

  Numeric transparency for lines. Default: 0.6.

- point_size:

  Numeric size of bin points. Default: 2.

- point_alpha:

  Numeric transparency for points. Default: 0.8.

- colour_by:

  Character: "n_contacts" or "none". Default: "n_contacts".

- palette:

  Character: ggplot2 colour palette. Default: "viridis".

- ...:

  Additional arguments passed to layer.

## Details

This geom creates a specialized visualization for multi-way chromatin
contacts from Pore-C/HiPore-C data. The visualization shows:

- **Y-axis**: Genomic bins (sorted by position)

- **X-axis**: Individual reads (ordered by first appearance)

- **Lines**: Connect all bins contacted by a single read

- **Points**: Individual bin positions

- **colour**: Number of contacts per read (optional)

The visualization is similar to a Sankey diagram but specialized for
multi-way genomic interactions.

## Aesthetics

`geom_hypergraph()` understands the following aesthetics: - `x`
(required): Read index or position - `y` (required): Bin position -
`group` (required): Read identifier - `colour` (optional): Number of
contacts per read

## Examples

``` r
if (FALSE) { # \dontrun{
# Build hypergraph
hg <- build_hypergraph(
  pairs = pairs_df,
  bin_size = 100000,
  chrom = "chr22",
  quantile = 0.85,
  min_multiway = 3
)

# Convert to tidy format
df <- hypergraph_to_tidy(hg, max_reads = 100)

# Create visualization
ggplot(df, aes(x = read_idx, y = bin, group = read_idx)) +
  geom_hypergraph(aes(colour = n_contacts)) +
  labs(title = "Multi-way Contacts on chr22") +
  theme_minimal()

# Without colouring
ggplot(df, aes(x = read_idx, y = bin, group = read_idx)) +
  geom_hypergraph(colour_by = "none") +
  theme_minimal()
} # }
```
