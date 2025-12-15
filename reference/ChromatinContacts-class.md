# ChromatinContacts S4 class

S4 class representing chromatin contact data from cooler files (.cool or
.mcool). Memory-efficient pointer that loads data only when `import()`
is called.

## Details

Key features:

- Memory-efficient: loads data only when needed

- Flexible focusing: subset regions before or after import

- Feature integration: automatically subsets attached features

- Multi-resolution: supports both .cool and .mcool files

Focus parameter options:

- Single region: `"chr1:1000000-2000000"`

- Whole chromosome: `"chr1"`

- Multiple regions (all combinations): `"chr1 | chr2"`

- Inter-chromosomal only: `"chr1 & chr2"`

- GInteractions object for programmatic control

## Slots

- `cooler_path`:

  Character. Path to cooler file.

- `resolution`:

  Integer or NULL. Resolution in base pairs (required for .mcool files).

- `focus`:

  GInteractions or NULL. Genomic regions of interest (NULL for
  genome-wide).

- `interactions`:

  GInteractions or NULL. Imported Hi-C interaction data.

- `seqinfo`:

  Seqinfo or NULL. Genome assembly sequence information.

- `compartments`:

  GRanges or NULL. A/B compartment annotations.

- `TADs`:

  GRanges or NULL. Topologically associating domains.

- `loops`:

  GInteractions or NULL. Chromatin loop calls.

- `multi_contacts`:

  GRanges or NULL. Multi-way contact regions.

- `tracks`:

  GRangesList. Genomic signal tracks (ChIP-seq, ATAC-seq, etc.).

## See also

[`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md),
[import-ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/import-ChromatinContacts.md),
[`features()`](https://jasonwong-lab.github.io/gghic/reference/features.md),
[`resolution()`](https://jasonwong-lab.github.io/gghic/reference/resolution.md),
[`focus()`](https://jasonwong-lab.github.io/gghic/reference/focus.md),
[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create and import
cc <- ChromatinContacts("data.cool") |> import()

# Focus on specific region
cc <- ChromatinContacts("data.cool", focus = "chr1:1-2e6") |> import()

# Multi-resolution cooler
cc <- ChromatinContacts("data.mcool", resolution = 10000) |> import()

# Add features
features(cc, "TADs") <- tads_granges
features(cc, "loops") <- loops_ginteractions
} # }
```
