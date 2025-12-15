# Import chromatin contact data

Imports interaction data from cooler file into ChromatinContacts object.

## Usage

``` r
# S4 method for class 'ChromatinContacts,ANY,ANY'
import(con, format, text, balance_column = "weight", ...)
```

## Arguments

- con:

  ChromatinContacts object with cooler file path.

- format:

  Character. Balance column name for normalization (default:
  `"weight"`). Can be provided as second argument for convenience.

- text:

  Not used (compatibility).

- balance_column:

  Character. Balance column name (overridden by `format` if provided).

- ...:

  Additional arguments (not used).

## Value

ChromatinContacts with `interactions` slot populated by GInteractions
containing bin_id1, bin_id2, count, and balanced columns.

## Details

Loads Hi-C data from cooler file. If focus regions are specified, only
those interactions are loaded. Balance weights (e.g., "weight", "KR",
"VC") are used to normalize raw counts. Warning issued if balance column
not found.

Shorthand: `import(cc, "KR")` is equivalent to
`import(cc, balance_column = "KR")`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic import
cc <- ChromatinContacts("sample.cool") |> import()

# With specific balance column
cc <- import(cc, "KR")

# Focus region (more efficient)
cc <- ChromatinContacts("sample.cool", focus = "chr1:1-5e6") |> import()
} # }
```
