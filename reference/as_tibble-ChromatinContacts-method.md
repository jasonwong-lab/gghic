# Convert ChromatinContacts to tibble

Converts interaction data or feature data from a `ChromatinContacts`
object to a tibble (data frame).

## Usage

``` r
# S4 method for class 'ChromatinContacts'
as_tibble(
  x,
  which = c("interactions", "compartments", "TADs", "loops", "multi_contacts")
)
```

## Arguments

- x:

  A `ChromatinContacts` object.

- which:

  Character string. Which data to convert: `"interactions"` (default)
  for Hi-C interactions, or a feature name like `"TADs"`, `"loops"`,
  `"compartments"`, or `"multi_contacts"`.

## Value

A tibble containing the requested data.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool") |> import()

# Convert interactions to tibble
df <- as_tibble(cc)

# Convert specific feature
tads_df <- as_tibble(cc, which = "TADs")
} # }
```
