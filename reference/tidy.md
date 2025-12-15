# Tidy hypergraph to long format

Converts hypergraph incidence matrix to long-format data frame for
analysis and visualization.

## Usage

``` r
tidy(x, ...)

# S4 method for class 'MultiWayContacts'
tidy(x, max_hyperedges = NULL, weight_normalization = "none")
```

## Arguments

- x:

  MultiWayContacts with built hypergraph.

- ...:

  Additional arguments (not used).

- max_hyperedges:

  Integer or NULL. Maximum hyperedges to include (default: NULL for
  all). Selects top by contact order.

- weight_normalization:

  Character. Weight scheme: `"none"` (raw), `"log"`, `"by_order"`, or
  `"minmax"` (default: `"none"`).

## Value

MultiWayContacts with `tidied_hypergraph` slot containing: bin_idx,
hyperedge_idx, bin_id, chrom, bin, n_multiways, weight.

## Details

Converts sparse matrix to long format. When `max_hyperedges` is set,
ranks by contact order and retains top N. Unused bins are removed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Tidy all with raw weights
mc <- MultiWayContacts("sample.pairs.gz") |> import() |> build() |> tidy()

# Top 100 with log weights
mc <- mc |> tidy(max_hyperedges = 100, weight_normalization = "log")
} # }
```
