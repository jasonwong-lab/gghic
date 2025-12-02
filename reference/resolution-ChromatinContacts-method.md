# Get resolution of ChromatinContacts object

Retrieves the resolution (bin size) in base pairs from a
`ChromatinContacts` object.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
resolution(x)
```

## Arguments

- x:

  A `ChromatinContacts` object.

## Value

Integer representing resolution in base pairs, or NULL if not set.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.mcool", resolution = 10000)
resolution(cc) # Returns 10000
} # }
```
