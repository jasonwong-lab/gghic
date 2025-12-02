# Get focus regions of ChromatinContacts object

Retrieves the focus regions (genomic regions of interest) from a
`ChromatinContacts` object.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
focus(x)
```

## Arguments

- x:

  A `ChromatinContacts` object.

## Value

GInteractions object representing the focus regions, or NULL if
genome-wide.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool", focus = "chr1")
focus(cc)
} # }
```
