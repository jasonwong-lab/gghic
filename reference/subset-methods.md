# Subset ChromatinContacts objects by various criteria

Flexibly subset ChromatinContacts objects using multiple indexing
methods: numeric indices, logical vectors, character strings (genomic
regions), GRanges, or GInteractions. All attached genomic features
(TADs, loops, tracks) are automatically filtered to match the subset.

## Usage

``` r
# S4 method for class 'ChromatinContacts,numeric'
subsetByOverlaps(x, ranges)

# S4 method for class 'ChromatinContacts,logical'
subsetByOverlaps(x, ranges)

# S4 method for class 'ChromatinContacts,GRanges'
subsetByOverlaps(x, ranges, type = c("within", "any"))

# S4 method for class 'ChromatinContacts,GInteractions'
subsetByOverlaps(x, ranges)

# S4 method for class 'ChromatinContacts,numeric,ANY,ANY'
x[i]

# S4 method for class 'ChromatinContacts,logical,ANY,ANY'
x[i]

# S4 method for class 'ChromatinContacts,GRanges,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'ChromatinContacts,GInteractions,ANY,ANY'
x[i]

# S4 method for class 'ChromatinContacts,character,ANY,ANY'
x[i]
```

## Arguments

- x:

  A `ChromatinContacts` object.

- ranges:

  For `subsetByOverlaps`: GRanges, GInteractions, numeric, or logical
  vector to subset by.

- type:

  Character string specifying overlap type for GRanges subsetting. Can
  be `"within"` or `"any"`. Default is `"any"`.

- i:

  Index, logical vector, character (region string), GRanges, or
  GInteractions for subsetting.

- j, drop:

  Ignored. Included for S4 compatibility.

- ...:

  Additional arguments passed to subsetting methods.

## Value

A subsetted `ChromatinContacts` object. Features (TADs, loops, tracks)
are automatically subsetted to match the new interaction data.

## Details

Subsetting methods:

- **Numeric**: Select interactions by index

- **Logical**: Select interactions matching TRUE values

- **Character**: Specify genomic region(s) as strings (e.g., `"chr1"` or
  `"chr1:1000000-2000000"`)

- **GRanges**: Select interactions overlapping the ranges

- **GInteractions**: Select interactions overlapping the interaction
  pairs

All methods automatically update:

- `interactions` slot

- `focus` slot

- All attached features (TADs, loops, tracks, etc.)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data with features
cc <- ChromatinContacts("sample.cool") |> import()
features(cc, "TADs") <- rtracklayer::import("tads.bed")
features(cc, "loops") <- loops_gi

# Subset by character region string
cc_chr1 <- cc["chr1"]  # Entire chromosome
cc_region <- cc["chr1:1000000-5000000"]  # Specific region
cc_multi <- cc["chr1|chr2"]  # Multiple regions

# Subset by numeric index
cc_first_1000 <- cc[1:1000]  # First 1000 interactions

# Subset by logical vector
gis <- interactions(cc)
high_counts <- gis$count > 100
cc_high <- cc[high_counts]  # Only high-count interactions

# Subset by GRanges
region_gr <- GRanges("chr1:1000000-2000000")
cc_overlap <- cc[region_gr]  # Interactions overlapping region

# Subset by GInteractions
specific_pairs <- GInteractions(
  GRanges("chr1:1000000-1100000"),
  GRanges("chr1:2000000-2100000")
)
cc_pairs <- cc[specific_pairs]

# All features are automatically subset
length(features(cc, "TADs"))  # Original TAD count
length(features(cc_chr1, "TADs"))  # Fewer TADs after subsetting

# By GRanges
roi <- GRanges("chr1:1000000-5000000")
cc_sub <- cc[roi]

# By GInteractions
query_gi <- GInteractions(
  GRanges("chr1:1000000-2000000"),
  GRanges("chr1:3000000-4000000")
)
cc_sub <- cc[query_gi]
} # }
```
