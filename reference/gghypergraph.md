# Visualize multi-way chromatin contacts as hypergraph network

Creates an intuitive ggplot2 visualization of multi-way chromatin
contacts from Pore-C or similar long-read technologies. Displays genomic
bins as points and multi-way contacts as connecting lines (hyperedges),
with line width and color encoding contact properties.

## Usage

``` r
gghypergraph(x, ...)

# S4 method for class 'MultiWayContacts'
gghypergraph(
  x,
  point_size = 2,
  line_width = 0.3,
  line_alpha = 0.6,
  color_by = "n_multiways",
  palette = "viridis",
  facet_chrom = TRUE,
  weight_normalization = "none"
)
```

## Arguments

- x:

  MultiWayContacts object with selected hyperedges (after running
  [`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md)).

- ...:

  Additional arguments (not used).

- point_size:

  Numeric. Size of genomic bin points (default: 2).

- line_width:

  Numeric. Base line width for hyperedges. Actual width scaled by
  hyperedge weight (default: 0.3).

- line_alpha:

  Numeric between 0 and 1. Line transparency to reduce overplotting
  (default: 0.6).

- color_by:

  Character. Variable for line coloring:

  - `"n_multiways"`: Color by contact order (2-way, 3-way, etc.)
    (default)

  - Other values: Uniform gray coloring

- palette:

  Character. Color palette name:

  - Viridis options: `"viridis"`, `"magma"`, `"plasma"`, `"inferno"`,
    `"cividis"`

  - ColorBrewer palettes: `"YlOrRd"`, `"Blues"`, etc.

  - Default: `"viridis"`

- facet_chrom:

  Logical. Layout strategy:

  - `TRUE`: Separate facet panel per chromosome (default)

  - `FALSE`: Composite y-axis with all chromosomes stacked

- weight_normalization:

  Character. Deprecated parameter (weights now set in
  [`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md)
  method).

## Value

A ggplot2 object that can be further customized with additional ggplot2
layers and themes.

## Details

### Visualization structure

- **X-axis**: Hyperedge index (each read/contact)

- **Y-axis**: Genomic bin position (genomic coordinate)

- **Points**: Individual genomic bins

- **Lines**: Connect all bins in a single multi-way contact

- **Line width**: Scaled by hyperedge weight

- **Line color**: Contact order (number of fragments per read)

### Interpretation

- Vertical lines: Bins contacted by same read

- Line thickness: Relative importance/weight

- Color gradient: Contact complexity (darker = higher-order)

- Patterns reveal genome organization and contact preferences

### Prerequisites

Must run complete workflow: `import()` →
[`build()`](https://jasonwong-lab.github.io/gghic/reference/build.md) →
[`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md) →
[`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md)
before visualization.

## See also

[`MultiWayContacts()`](https://jasonwong-lab.github.io/gghic/reference/MultiWayContacts-class.md),
[`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md),
[`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md),
[`build()`](https://jasonwong-lab.github.io/gghic/reference/build.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Complete workflow
mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1") |>
  import() |>
  build(bin_size = 1000000L, quantile = 0.85) |>
  tidy(weight_normalization = "log") |>
  select(n_intra = 10, n_inter = 5)

# Basic visualization
gghypergraph(mc)

# Composite layout with custom colors
gghypergraph(mc, facet_chrom = FALSE, palette = "magma")

# Emphasize line patterns
gghypergraph(mc, line_width = 0.5, line_alpha = 0.8, point_size = 3)

# Further customization
gghypergraph(mc) +
  theme_minimal() +
  labs(title = "Multi-way Contacts on Chromosome 1")
} # }
```
