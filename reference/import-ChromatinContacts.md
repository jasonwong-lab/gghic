# Import chromatin contact data

Imports chromatin interaction data from a cooler file into a
ChromatinContacts object. The method loads bin-level interactions and
applies balancing weights if available.

## Usage

``` r
# S4 method for class 'ChromatinContacts,ANY,ANY'
import(con, format, text, balance_column = "weight", ...)
```

## Arguments

- con:

  A `ChromatinContacts` object containing the path to the cooler file
  and optional focus regions.

- format:

  Character string specifying the balance column name to use for
  normalization. This parameter is used as the second positional
  argument for convenience. Default is `"weight"`.

- text:

  Not used. Included for compatibility with the generic.

- balance_column:

  Character string specifying the balance column name to use for
  normalization. If `format` is provided, it overrides this parameter.
  Default is `"weight"`.

- ...:

  Additional arguments (currently not used).

## Value

A `ChromatinContacts` object with the `interactions` slot populated by a
`GInteractions` object containing:

- bin_id1: Bin ID for the first anchor

- bin_id2: Bin ID for the second anchor

- count: Raw interaction count

- balanced: Normalized interaction count (if balance weights are
  available)

## Details

The function loads chromatin interaction data from a cooler file. If a
focus region is specified in the `ChromatinContacts` object, only
interactions within or between those regions are loaded. Otherwise, all
genome-wide interactions are imported.

The balance column (e.g., "weight", "KR", "VC") is used to normalize the
raw counts. If the specified balance column is not found in the cooler
file, a warning is issued and unbalanced counts are returned.

For convenience, the balance column name can be provided as the second
positional argument: `import(cc, "KR")` is equivalent to
`import(cc, balance_column = "KR")`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic import from cooler file
cc <- ChromatinContacts("sample.cool")
cc <- import(cc)

# Using pipe operator
cc <- ChromatinContacts("sample.cool") |> import()

# With specific balance column
cc <- import(cc, balance_column = "KR")

# Convenient syntax for balance column
cc <- import(cc, "VC")

# Import with focus region (more efficient)
cc <- ChromatinContacts(
  "sample.cool",
  focus = "chr1:1000000-5000000"
) |> import()

# Multi-resolution cooler
cc <- ChromatinContacts(
  "sample.mcool",
  resolution = 10000
) |> import()

# Check imported data
interactions(cc)
length(interactions(cc))
} # }
```
