# Get sequence information from ChromatinContacts object

Retrieves the sequence information (chromosome names, lengths) from a
`ChromatinContacts` object.

## Usage

``` r
# S4 method for class 'ChromatinContacts'
seqinfo(x)
```

## Arguments

- x:

  A `ChromatinContacts` object.

## Value

A `Seqinfo` object containing sequence/chromosome information.

## Examples

``` r
if (FALSE) { # \dontrun{
cc <- ChromatinContacts("sample.cool")
seqinfo(cc)
} # }
```
