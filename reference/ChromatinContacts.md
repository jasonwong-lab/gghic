# Create a ChromatinContacts object

Constructor for ChromatinContacts from cooler files (.cool or .mcool).

## Usage

``` r
ChromatinContacts(cooler_path, focus = NULL, resolution = NULL)
```

## Arguments

- cooler_path:

  Character. Path to .cool or .mcool file.

- focus:

  Character, GInteractions, or NULL. Genomic region(s) to focus:

  - Single region: `"chr1"` or `"chr1:1-1000000"`

  - Inter-chromosomal: `"chr1&chr2"` (excludes diagonal)

  - All combinations: `"chr1|chr2"` (includes intra-chromosomal)

  - Multiple: `c("chr1&chr2", "chr3|chr4")`

  - NULL for genome-wide (default)

- resolution:

  Integer. Resolution in base pairs (required for .mcool).

## Value

A ChromatinContacts object.

## See also

[import-ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/import-ChromatinContacts.md),
[`resolution()`](https://jasonwong-lab.github.io/gghic/reference/resolution.md),
[`focus()`](https://jasonwong-lab.github.io/gghic/reference/focus.md),
[`features()`](https://jasonwong-lab.github.io/gghic/reference/features.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Genome-wide .cool file
cc <- ChromatinContacts("file.cool")

# Multi-resolution .mcool
cc <- ChromatinContacts("file.mcool", resolution = 5000L)

# Focus on chromosome
cc <- ChromatinContacts("file.cool", focus = "chr1")

# Inter-chromosomal between chr1 and chr2
cc <- ChromatinContacts("file.cool", focus = "chr1&chr2")
} # }
```
