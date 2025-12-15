# Get interaction data

Retrieves imported Hi-C interactions.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
interactions(x)
```

## Arguments

- x:

  ChromatinContacts object.

## Value

GInteractions (Hi-C data) or NULL (not yet imported).

## Details

Must load with `import()` first. GInteractions contains bin coordinates
and metadata (count, balanced).

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool") |> import()
gis <- interactions(cc)
} # }
```
