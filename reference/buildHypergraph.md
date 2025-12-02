# Build Hypergraph from Multi-way Contacts

Analyze multi-way contacts from Pore-C/HiPore-C data by constructing
hypergraph representations where genomic bins are nodes and reads are
hyperedges connecting multiple bins.

## Usage

``` r
buildHypergraph(
  pairs = NULL,
  pairs_file = NULL,
  bin_size,
  chrom = NULL,
  min_contacts = NULL,
  quantile = 0.85,
  min_multiway = 3,
  inter_chrom = FALSE
)
```

## Arguments

- pairs:

  A data frame or tibble with columns: `read_name`, `chrom1`, `pos1`,
  `chrom2`, `pos2`. For large files, use `pairs_file` instead.

- pairs_file:

  Character path to pairs file for chunked reading.

- bin_size:

  Integer bin size in base pairs.

- chrom:

  Character chromosome name(s) to analyze (e.g., "chr22" or c("chr21",
  "chr22")). Use NULL for genome-wide analysis. Default: NULL.

- min_contacts:

  Integer minimum number of contacts per bin pair. Filters out sparse
  pairwise interactions. Default: NULL (no filtering).

- quantile:

  Numeric quantile threshold (0-1) for filtering pairwise contacts.
  Keeps only bin pairs above this quantile. Default: 0.85.

- min_multiway:

  Integer minimum number of bins a read must contact to be included in
  hypergraph. Default: 3.

- inter_chrom:

  Logical, whether to include inter-chromosomal contacts. When FALSE,
  only intra-chromosomal contacts are kept. When TRUE and multiple
  chromosomes are specified in `chrom`, includes both intra- and inter-
  chromosomal contacts between the specified chromosomes only (not with
  other chromosomes in the genome). Default: FALSE (intra-chromosomal
  only).

## Value

A list with class "hypergraph" containing:

- incidence:

  Sparse incidence matrix (bins x reads)

- bins:

  Integer vector of bin IDs (sorted)

- reads:

  Character vector of read names

- contacts_per_read:

  Integer vector of contact counts per read

- bin_size:

  The bin size used

- chrom:

  The chromosome analyzed

## Details

The function performs the following steps:

1.  Bins genomic positions into fixed-size bins

2.  Removes duplicate pairwise contacts within each read

3.  Filters bin pairs by count quantile or minimum threshold

4.  Constructs incidence matrix (bins as rows, reads as columns)

5.  Filters reads by minimum multi-way contact degree

The resulting hypergraph can be visualized with
[`geom_hypergraph()`](https://jasonwong-lab.github.io/gghic/reference/geom_hypergraph.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Build hypergraph from pairs data
hg <- buildHypergraph(
  pairs = pairs_df,
  bin_size = 100000,
  chrom = "chr22",
  quantile = 0.85,
  min_multiway = 3
)

# From large file
hg <- buildHypergraph(
  pairs_file = "contact_matrix.txt",
  bin_size = 100000,
  chrom = "chr22"
)

# Multi-chromosome with inter-chromosomal contacts
hg <- buildHypergraph(
  pairs = pairs_df,
  bin_size = 100000,
  chrom = c("chr21", "chr22"),
  inter_chrom = TRUE
)

# Genome-wide with inter-chromosomal contacts
hg <- buildHypergraph(
  pairs = pairs_df,
  bin_size = 100000,
  chrom = NULL,
  inter_chrom = TRUE
)

# Visualize
plotHypergraph(hg)
} # }
```
