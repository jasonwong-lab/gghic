# Get resolution of ChromatinContacts object

Retrieves bin size in base pairs.

## Usage

``` r
resolution(x)

# S4 method for class 'ChromatinContacts'
resolution(x)
```

## Arguments

- x:

  ChromatinContacts object.

## Value

Integer (resolution in bp) or NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.mcool", resolution = 10000L)
resolution(cc) # Returns 10000
} # }
```
