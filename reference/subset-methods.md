# Subset ChromatinContacts objects

Subset `ChromatinContacts` objects using various methods: numeric
indices, logical vectors, character strings (genomic regions), GRanges,
or GInteractions.

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
cc <- ChromatinContacts("sample.cool") |> import()

# By character region
cc_sub <- cc["chr1"]
cc_sub <- cc["chr1:1000000-5000000"]

# By numeric index
cc_sub <- cc[1:1000]

# By logical vector
high_counts <- interactions(cc)$count > 100
cc_sub <- cc[high_counts]

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
