# ChromatinContacts S4 class

An S4 class to represent chromatin contact data from Hi-C experiments
stored in cooler format (.cool or .mcool files). This lightweight object
acts as a pointer to Hi-C data without loading the full matrix into
memory until explicitly imported.

## Details

The `ChromatinContacts` class provides an efficient framework for
managing Hi-C data. Key features:

- **Memory efficient**: Only loads data when `import()` is called

- **Flexible focusing**: Subset regions before or after loading

- **Feature integration**: Attach TADs, loops, and tracks that
  auto-subset with the interaction data

- **Multi-resolution support**: Works with both .cool and .mcool files

The `focus` parameter accepts flexible specifications:

- Single region: `"chr1:1000000-2000000"`

- Whole chromosome: `"chr1"`

- Multiple regions with OR: `"chr1 | chr2"` (includes all interactions)

- Inter-regional with AND: `"chr1 & chr2"` (inter-chromosomal only)

- GInteractions object for programmatic definition

## Slots

- `cooler_path`:

  Character string. Path to the .cool or .mcool file.

- `resolution`:

  Integer or NULL. Resolution in base pairs for .mcool files. Required
  for multi-resolution cooler files.

- `focus`:

  GInteractions or NULL. Genomic regions of interest to focus on. Can be
  specified as a GInteractions object, or NULL for genome-wide data.

- `interactions`:

  GInteractions or NULL. The imported Hi-C interaction data. Initially
  NULL until `import()` is called.

- `seqinfo`:

  Seqinfo or NULL. Sequence information for the genome assembly.

- `compartments`:

  GRanges or NULL. A/B compartment annotations.

- `TADs`:

  GRanges or NULL. Topologically associating domains (TADs).

- `loops`:

  GInteractions or NULL. Chromatin loop calls.

- `multi_contacts`:

  GRanges or NULL. Multi-way contact regions (e.g., from Pore-C data).

- `tracks`:

  GRangesList. Additional genomic signal tracks (e.g., ChIP-seq,
  ATAC-seq).

## See also

- [`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md) -
  Constructor function

- [import-ChromatinContacts](https://jasonwong-lab.github.io/gghic/reference/import-ChromatinContacts.md) -
  Load interaction data

- [`features()`](https://jasonwong-lab.github.io/gghic/reference/generics.md) -
  Access or set genomic features

- [`resolution()`](https://jasonwong-lab.github.io/gghic/reference/generics.md) -
  Get resolution

- [`focus()`](https://jasonwong-lab.github.io/gghic/reference/generics.md) -
  Get focus regions

- [`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md) -
  Create visualizations

## Examples

``` r
if (FALSE) { # \dontrun{
# Create object pointing to cooler file
cc <- ChromatinContacts(cooler_path = "data.cool")

# Create with focus on specific region
cc <- ChromatinContacts(
  cooler_path = "data.cool",
  focus = "chr1:1000000-2000000"
)

# Multi-resolution cooler with specific resolution
cc <- ChromatinContacts(
  cooler_path = "data.mcool",
  resolution = 10000
)

# Import the data
cc <- import(cc)

# Add features
features(cc, "TADs") <- tads_granges
features(cc, "loops") <- loops_ginteractions
features(cc, "tracks") <- tracks_grangeslist
} # }
```
