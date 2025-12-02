# Set genomic features in ChromatinContacts object

Assigns genomic features (TADs, compartments, multi_contacts, or tracks)
to a `ChromatinContacts` object. The features are automatically confined
to overlap with the imported interaction regions.

Assigns chromatin loop data to a `ChromatinContacts` object. Loops are
automatically confined to overlap with the imported interaction regions.

Assigns genomic signal tracks (e.g., ChIP-seq, ATAC-seq) to a
`ChromatinContacts` object. Tracks are automatically confined to overlap
with the imported interaction regions.

## Usage

``` r
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

  A `ChromatinContacts` object.

- name:

  Character string. Must be `"tracks"`.

- ...:

  Additional arguments. Can include `which` to specify track names.

- value:

  A `GRangesList` object containing track data. Each element should be a
  `GRanges` object representing a single track.

## Value

The modified `ChromatinContacts` object.

The modified `ChromatinContacts` object.

The modified `ChromatinContacts` object.

## Details

Features are automatically subsetted to only include regions that
overlap with the imported interaction data. This ensures consistency
when subsetting the `ChromatinContacts` object.

**Note:** For tracks, use `GRangesList` instead. If a single `GRanges`
is provided for tracks, it will be converted to a `GRangesList` with one
element.

**Important:** Interactions must be imported before features can be
assigned.

Loops are stored as `GInteractions` objects representing paired genomic
regions. They are automatically subsetted to only include loops that
overlap with the imported interaction data.

**Important:** Interactions must be imported before loops can be
assigned.

Tracks should be provided as a `GRangesList` where each element
represents a different signal track. If the list is unnamed, default
names (`track1`, `track2`, etc.) will be assigned.

**Important:** Interactions must be imported before tracks can be
assigned.

## See also

[`features()`](https://jasonwong-lab.github.io/gghic/reference/generics.md)

[`features()`](https://jasonwong-lab.github.io/gghic/reference/generics.md)

[`features()`](https://jasonwong-lab.github.io/gghic/reference/generics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Import interactions first
cc <- ChromatinContacts("sample.cool") |> import()

# Add TADs
tads <- rtracklayer::import("TADs.bed")
features(cc, "TADs") <- tads

# Add compartments
features(cc, "compartments") <- compartments_gr
} # }

if (FALSE) { # \dontrun{
# Import interactions first
cc <- ChromatinContacts("sample.cool") |> import()

# Add loops
loops <- rtracklayer::import("loops.bedpe") |>
  makeGInteractionsFromGRangesPairs()
features(cc, "loops") <- loops
} # }

if (FALSE) { # \dontrun{
# Import interactions first
cc <- ChromatinContacts("sample.cool") |> import()

# Add tracks from BigWig files
track1 <- rtracklayer::import("signal1.bw")
track2 <- rtracklayer::import("signal2.bw")
tracks <- GRangesList(ChIP1 = track1, ChIP2 = track2)
features(cc, "tracks") <- tracks_list
} # }
```
