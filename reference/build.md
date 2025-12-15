# Build hypergraph from multi-way contacts

Constructs hypergraph from pairs data by binning, filtering, and
building sparse incidence matrix.

## Usage

``` r
build(x, ...)

# S4 method for class 'MultiWayContacts'
build(x, bin_size = 1000000L, quantile = 0.85, min_multiway = 2)
```

## Arguments

- x:

  MultiWayContacts object with loaded pairs.

- ...:

  Additional arguments (not used).

- bin_size:

  Integer. Bin size in base pairs (default: 1000000L).

- quantile:

  Numeric (0-1). Quantile threshold for filtering bin pairs by contact
  frequency (default: 0.85).

- min_multiway:

  Integer. Minimum contacts per hyperedge to retain (default: 2).

## Value

MultiWayContacts with populated hypergraph slots: hypergraph (incidence
matrix), weights (normalizations), bin_info, hyperedge_reads, multiways,
and bin_size.

## Details

Steps: bin positions → filter by quantile → build incidence matrix →
deduplicate patterns → filter by min_multiway.

Weight normalizations: raw (frequency), log (log-transformed), by_order
(normalized within order), minmax (scaled to \[0,1\]).

Higher quantile = stricter filtering. Typical: 0.8-0.95.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default parameters
mc <- MultiWayContacts("sample.pairs.gz") |> import() |> build()

# Custom parameters
mc <- mc |> build(bin_size = 500000L, quantile = 0.9, min_multiway = 3)
} # }
```
