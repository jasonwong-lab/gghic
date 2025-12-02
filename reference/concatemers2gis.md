# Convert concatemer reads to pairwise interactions

Converts multi-way contact reads (concatemers) from Pore-C or other
multi-contact sequencing to pairwise interactions in a `GInteractions`
object. Useful for comparing multi-way data with traditional Hi-C.

## Usage

``` r
concatemers2Gis(grs, region = NULL, bin_size = 1e+05, read_group = "read_name")
```

## Arguments

- grs:

  A `GRanges` object containing concatemer reads. Must have a metadata
  column (default `read_name`) identifying which fragments belong to the
  same read.

- region:

  A `GRanges` object specifying the genomic region of interest. Default
  is `NULL`, which uses the entire range of `grs`.

- bin_size:

  Integer. The size of bins (in base pairs) to divide the genomic region
  into. Default is `100000` (100 kb).

- read_group:

  Character string. The name of the metadata column in `grs` that
  defines read groups (typically read names or IDs). Default is
  `"read_name"`.

## Value

A `GInteractions` object containing all pairwise combinations of
fragments within each concatemer read, binned at the specified
resolution.

## Details

For each multi-way contact (concatemer), this function generates all
possible pairwise combinations of the fragments and bins them to the
specified resolution. This allows visualization of Pore-C or other
multi-contact data alongside traditional Hi-C data.

The resulting interaction counts represent the number of concatemer
reads supporting each pairwise interaction, which differs from
traditional Hi-C counting.

**Workflow:**

1.  Groups fragments by read identifier

2.  For each read with multiple fragments, generates all pairwise
    combinations

3.  Bins fragment coordinates to specified resolution

4.  Counts interactions per bin pair

## Examples

``` r
if (FALSE) { # \dontrun{
# Load concatemer data
concatemers <- readRDS("concatemers.rds")

# Convert to pairwise interactions at 100kb resolution
gis <- concatemers2Gis(concatemers, bin_size = 100000)

# Visualize
gghic(gis)

# Focus on specific region
region <- GRanges("chr1:1000000-5000000")
gis_sub <- concatemers2Gis(concatemers, region = region, bin_size = 10000)
} # }
```
