# Select top hyperedges for visualization

Selects top-weighted intra- and inter-chromosomal hyperedges per
chromosome for visualization.

## Usage

``` r
select(x, ...)

# S4 method for class 'MultiWayContacts'
select(
  x,
  n_intra = 5,
  n_inter = 5,
  n_multiways_filter = NULL,
  chroms = NULL,
  append = TRUE,
  ...
)
```

## Arguments

- x:

  MultiWayContacts with tidied hypergraph.

- ...:

  Additional arguments (not used).

- n_intra:

  Integer. Top intra-chromosomal hyperedges per chromosome (default: 5).

- n_inter:

  Integer. Top inter-chromosomal hyperedges per chromosome (default: 5).
  Quota increases if fewer intra hyperedges available.

- n_multiways_filter:

  Integer vector or NULL. Filter by specific contact orders, e.g.,
  `c(3,4)` for 3-way and 4-way only (default: NULL).

- chroms:

  Character or NULL. Filter by chromosome(s) (default: NULL).

- append:

  Logical. Append to existing selection (TRUE) or replace (FALSE,
  default: TRUE).

## Value

MultiWayContacts with `select_hypergraph` slot containing selected
hyperedges plus `type` (intra/inter) and `source_chrom` columns.

## Details

Per chromosome: select top `n_intra` intra-chromosomal and `n_inter`
inter-chromosomal hyperedges by weight. If insufficient intra
hyperedges, quota transferred to inter. Hyperedges deduplicated across
chromosomes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: top 5 intra + 5 inter per chromosome
mc <- MultiWayContacts("sample.pairs.gz") |>
  import() |> build() |> tidy() |> select()

# More hyperedges
mc <- mc |> select(n_intra = 10, n_inter = 10, append = FALSE)

# Filter by contact order
mc <- mc |> select(n_multiways_filter = c(3, 4), append = FALSE)

# Specific chromosomes
mc <- mc |> select(chroms = c("chr1", "chr2"), append = FALSE)
} # }
```
