# Set genomic features in ChromatinContacts object

Assigns or updates genomic features (TADs, loops, compartments, tracks)
in a ChromatinContacts object. Features are automatically subsetted to
overlap with loaded interactions, ensuring only relevant data is
retained.

## Usage

``` r
features(x, name, ...) <- value

# S4 method for class 'ChromatinContacts,character,GRanges'
features(x, name, ...) <- value

# S4 method for class 'ChromatinContacts,character,GInteractions'
features(x, name, ...) <- value

# S4 method for class 'ChromatinContacts,character,GRangesList'
features(x, name, ...) <- value

# S4 method for class 'ChromatinContacts,character,list'
features(x, name, ...) <- value
```

## Arguments

- x:

  ChromatinContacts object with imported interaction data.

- name:

  Character. Feature slot name:

  - `"compartments"`: A/B compartment annotations (GRanges)

  - `"TADs"`: Topologically associating domains (GRanges)

  - `"multi_contacts"`: Multi-way contact regions (GRanges)

  - `"loops"`: Chromatin loop calls (GInteractions)

  - `"tracks"`: Genomic signal tracks (GRangesList)

- ...:

  Additional arguments (not currently used).

- value:

  Feature data with appropriate type:

  - GRanges: for compartments, TADs, multi_contacts

  - GInteractions: for loops

  - GRangesList or list: for tracks (will be converted to GRangesList)

## Value

ChromatinContacts object with updated features slot. The modified object
is returned invisibly to support piping.

## Details

### Requirements

Interactions must be imported before assigning features. Use `import()`
first.

### Automatic subsetting

Features are automatically filtered to retain only elements overlapping
the loaded Hi-C interaction data. This prevents visualization errors and
reduces memory usage.

### Track naming

When providing tracks as GRangesList, ensure they are named for proper
labeling in plots. Unnamed tracks receive default names (`track1`,
`track2`, etc.).

## See also

[`features()`](https://jasonwong-lab.github.io/gghic/reference/features.md),
[`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md),
[import-ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/import-ChromatinContacts.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load Hi-C data first
cc <- ChromatinContacts("sample.cool") |> import()

# Add TADs from BED file
features(cc, "TADs") <- rtracklayer::import("TADs.bed")

# Add chromatin loops
features(cc, "loops") <- loops_ginteractions

# Add multiple ChIP-seq tracks
features(cc, "tracks") <- GRangesList(
  H3K27ac = rtracklayer::import("H3K27ac.bw"),
  H3K4me3 = rtracklayer::import("H3K4me3.bw")
)

# Add compartments
features(cc, "compartments") <- compartment_gr

# Chain operations
cc <- ChromatinContacts("sample.cool") |>
  import() |>
  `features<-`("TADs", value = tads) |>
  `features<-`("loops", value = loops)
} # }
```
