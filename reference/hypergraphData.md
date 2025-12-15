# Get hypergraph data from MultiWayContacts object

Retrieves the selected or tidied hypergraph data from a
`MultiWayContacts` object.

## Usage

``` r
hypergraphData(x, which = c("selected", "tidied"))

# S4 method for class 'MultiWayContacts'
hypergraphData(x, which = c("selected", "tidied"))
```

## Arguments

- x:

  A `MultiWayContacts` object.

- which:

  Character. Which hypergraph data to retrieve: `"selected"` (default)
  returns the selected subset for visualization, `"tidied"` returns the
  full tidied long-format data.

## Value

A data.frame containing hypergraph data in long format, or NULL if the
requested data has not been computed yet.

## Details

- `which = "selected"`: Returns data from
  [`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md),
  a subset of hyperedges chosen for visualization. Use this for
  plotting.

- `which = "tidied"`: Returns data from
  [`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md),
  the full hypergraph in long format. Use this for analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
mc <- MultiWayContacts("sample.pairs.gz") |>
  import() |>
  build(bin_size = 1e6, min_multiway = 3) |>
  tidy(max_hyperedges = 100) |>
  select(n = 50)

# Get selected hyperedges for plotting
df <- hypergraphData(mc, "selected")

# Get full tidied data for analysis
df_full <- hypergraphData(mc, "tidied")
} # }
```
