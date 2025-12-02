# Create a ChromatinContacts object

Constructor function for creating a ChromatinContacts object from a
cooler file (.cool or .mcool format).

## Usage

``` r
ChromatinContacts(cooler_path, focus = NULL, resolution = NULL)
```

## Arguments

- cooler_path:

  Character string. Path to a .cool or .mcool file.

- focus:

  Character string, GInteractions object, or NULL. Specifies the genomic
  region(s) to focus on. Can be specified as:

  - A single region: `"chr1"` or `"chr1:1-1000000"`

  - Multiple regions with `&` (inter-chromosomal pairs only):
    `"chr1&chr2"` or `"chr1:1-1000000&chr2:1-2000000"`

  - Multiple regions with `|` (all combinations including intra):
    `"chr1|chr2"` or `"chr1:1-1000000|chr2:1-2000000"`

  - Multiple specs as a character vector: c("chr1&chr2",
    "chr3:1-2e6\|chr4", "chr5")

  - (Legacy) Single comma-separated string: "chr1&chr2,chr3\|chr4"

  - A GInteractions object

  - NULL for genome-wide (default)

- resolution:

  Integer. Resolution in base pairs. Required for .mcool files, ignored
  for .cool files.

## Value

A ChromatinContacts object.

## Details

The focus parameter allows flexible specification of regions of
interest:

- `&` operator: Creates inter-chromosomal pairs only (excludes diagonal)

- `|` operator: Creates all combinations including intra-chromosomal

- Multiple specifications can be combined with commas

## See also

[import-ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/import-ChromatinContacts.md),
[`Seqinfo::seqinfo()`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html),
[`resolution()`](https://jasonwong-lab.github.io/gghic/reference/generics.md),
[`focus()`](https://jasonwong-lab.github.io/gghic/reference/generics.md),
[`features()`](https://jasonwong-lab.github.io/gghic/reference/generics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create object from .cool file (genome-wide)
cc <- ChromatinContacts("path/to/file.cool")

# Create object from .mcool file with specific resolution
cc <- ChromatinContacts("path/to/file.mcool", resolution = 5000L)

# Focus on a single chromosome
cc <- ChromatinContacts("path/to/file.cool", focus = "chr1")

# Focus on a specific region
cc <- ChromatinContacts("path/to/file.cool", focus = "chr1:1-1000000")

# Focus on inter-chromosomal interactions between chr1 and chr2
cc <- ChromatinContacts("path/to/file.cool", focus = "chr1&chr2")

# Focus on all combinations of chr1 and chr2 (including intra)
cc <- ChromatinContacts("path/to/file.cool", focus = "chr1|chr2")

# Multiple focus specifications
cc <- ChromatinContacts(
  "path/to/file.cool",
  focus = "chr1&chr2,chr3:1-1000000|chr4"
)

# Import the interaction data
cc <- import(cc)
} # }
```
