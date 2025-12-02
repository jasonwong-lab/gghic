# Get interaction data from ChromatinContacts object

Retrieves the imported Hi-C interaction data as a `GInteractions`
object.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
interactions(x)
```

## Arguments

- x:

  A `ChromatinContacts` object.

## Value

A `GInteractions` object containing the Hi-C interactions, or NULL if
data has not been imported yet.

## Details

The interactions must be loaded first using `import()`. The returned
`GInteractions` object contains bin coordinates and metadata columns
such as `count` and `balanced`.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool") |> import()
gis <- interactions(cc)
} # }
```
