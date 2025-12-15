# Get focus regions

Retrieves focus regions from ChromatinContacts object.

## Usage

``` r
focus(x)

# S4 method for class 'ChromatinContacts'
focus(x)
```

## Arguments

- x:

  ChromatinContacts object.

## Value

GInteractions (focus regions) or NULL (genome-wide).

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool", focus = "chr1")
focus(cc)
} # }
```
