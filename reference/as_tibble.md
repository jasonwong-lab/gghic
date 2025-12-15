# Convert ChromatinContacts to tibble

Converts interaction or feature data from ChromatinContacts to tibble.

## Usage

``` r
as_tibble(x, which = "interactions")

# S4 method for class 'ChromatinContacts'
as_tibble(
  x,
  which = c("interactions", "compartments", "TADs", "loops", "multi_contacts")
)
```

## Arguments

- x:

  ChromatinContacts object.

- which:

  Character. Data to convert: `"interactions"` (default), or feature
  name (`"TADs"`, `"loops"`, `"compartments"`, `"multi_contacts"`).

## Value

Tibble with requested data.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool") |> import()
df <- as_tibble(cc)
tads_df <- as_tibble(cc, which = "TADs")
} # }
```
