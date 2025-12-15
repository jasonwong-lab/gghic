# Import multi-way contact data

Imports multi-way chromatin contact data from a .pairs file into a
MultiWayContacts object. The method loads pairwise contacts from
long-read sequencing experiments (e.g., Pore-C).

## Usage

``` r
# S4 method for class 'MultiWayContacts,ANY,ANY'
import(con, format, text, inter_chrom = FALSE, ...)
```

## Arguments

- con:

  A `MultiWayContacts` object containing the path to the pairs file and
  optional focus regions.

- format:

  Logical indicating whether to include inter-chromosomal contacts. This
  parameter is used as the second positional argument for convenience.
  Default is `FALSE`.

- text:

  Not used. Included for compatibility with the generic.

- inter_chrom:

  Logical. If `TRUE`, includes inter-chromosomal contacts. If `FALSE`
  (default), only intra-chromosomal contacts are loaded. If `format` is
  provided, it overrides this parameter.

- ...:

  Additional arguments (currently not used).

## Value

A `MultiWayContacts` object with the `pairs` slot populated by a data
frame containing:

- read_name: Unique identifier for each read

- chrom1, chrom2: Chromosomes for each contact

- pos1, pos2: Genomic positions for each contact

## Details

The function loads pairwise contact data from a .pairs file format
(typically .pairs.gz). If a focus chromosome is specified in the
`MultiWayContacts` object, only contacts involving that chromosome are
loaded.

When `inter_chrom = FALSE`, reads that contain any inter-chromosomal
contacts are completely removed from the dataset to ensure clean
intra-chromosomal analysis.

The pairs file is read using optimized C code for efficient parsing of
large files.

For convenience, the inter_chrom parameter can be provided as the second
positional argument: `import(mc, TRUE)` is equivalent to
`import(mc, inter_chrom = TRUE)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic import (intra-chromosomal only)
mc <- MultiWayContacts("sample.pairs.gz")
mc <- import(mc)

# Using pipe operator
mc <- MultiWayContacts("sample.pairs.gz") |> import()

# Include inter-chromosomal contacts
mc <- import(mc, inter_chrom = TRUE)

# Convenient syntax
mc <- import(mc, TRUE)

# Import with focus on specific chromosome
mc <- MultiWayContacts(
  "sample.pairs.gz",
  focus = "chr1"
) |> import()

# Check imported data
head(mc@pairs)
} # }
```
