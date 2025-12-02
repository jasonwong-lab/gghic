# Plot Hypergraph Visualization

Creates a hypergraph visualization showing bins as nodes (circles) on
the y-axis and reads as vertical lines connecting multiple bins.

## Usage

``` r
plotHypergraph(
  hg,
  chrom = NULL,
  max_reads = 100,
  point_size = 2,
  line_width = 0.3,
  line_alpha = 0.6,
  color_by = "n_contacts",
  palette = "viridis",
  facet_chrom = TRUE
)
```

## Arguments

- hg:

  A hypergraph object from
  [`buildHypergraph()`](https://jasonwong-lab.github.io/gghic/reference/buildHypergraph.md).

- chrom:

  Character vector of chromosome names to plot. If NULL (default), plots
  all chromosomes in the hypergraph. Use this to subset chromosomes from
  a multi-chromosome hypergraph.

- max_reads:

  Integer maximum number of reads to display. Default: 100.

- point_size:

  Numeric size of bin points. Default: 2.

- line_width:

  Numeric width of read connection lines. Default: 0.3.

- line_alpha:

  Numeric transparency of lines. Default: 0.6.

- color_by:

  Character, how to color reads: "n_contacts" (default) or "none".

- palette:

  Character color palette. Default: "viridis".

- facet_chrom:

  Logical, whether to facet by chromosome for multi-chrom hypergraphs.
  Default: TRUE.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
hg <- buildHypergraph(pairs, bin_size = 1e5, chrom = "chr22")
plotHypergraph(hg, max_reads = 50)

# Multi-chromosome
hg <- buildHypergraph(pairs, bin_size = 1e5, chrom = c("chr21", "chr22"))
plotHypergraph(hg, max_reads = 100, facet_chrom = TRUE)

# Plot only specific chromosomes from multi-chromosome hypergraph
plotHypergraph(hg, chrom = "chr21", max_reads = 50)
plotHypergraph(hg, chrom = c("chr21", "chr22"), facet_chrom = FALSE)
} # }
```
