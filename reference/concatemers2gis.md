# Convert multi-way concatemer reads to pairwise Hi-C interactions

Transforms multi-way chromatin contact reads (concatemers from Pore-C,
Tri-C, Multi-C, etc.) into traditional pairwise Hi-C format. Generates
all pairwise combinations of fragments within each read, bins to
specified resolution, and creates a GInteractions object compatible with
standard Hi-C analysis tools.

## Usage

``` r
concatemers2Gis(grs, region = NULL, bin_size = 1e+05, read_group = "read_name")
```

## Arguments

- grs:

  GRanges object containing concatemer fragment coordinates. Must
  include metadata column identifying which fragments belong to the same
  read (specified by `read_group` parameter).

- region:

  GRanges object defining genomic region(s) of interest. Only reads
  overlapping this region are included (default: NULL for genome-wide
  analysis).

- bin_size:

  Integer. Genomic bin size in base pairs for aggregating contacts.
  Typical values: 5000-1000000 (default: 100000 = 100kb).

- read_group:

  Character. Name of metadata column in `grs` that contains read
  identifiers grouping fragments from the same sequencing read (default:
  `"read_name"`).

## Value

GInteractions object containing:

- Binned pairwise interactions from all concatemers

- `count` metadata column: number of reads supporting each interaction

- Compatible with
  [`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
  [`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md),
  and other Hi-C tools

## Details

### Algorithm

1.  Group fragments by read identifier (`read_group`)

2.  For each read with ≥2 fragments, generate all pairwise combinations

3.  Bin both anchors of each pair to specified resolution

4.  Aggregate identical bin pairs, counting supporting reads

5.  Remove self-interactions (same bin pairs)

### Multi-way to pairwise conversion

A read contacting N fragments generates N×(N-1)/2 pairwise interactions.
For example:

- 2-way contact → 1 pair

- 3-way contact → 3 pairs

- 4-way contact → 6 pairs

This inflation should be considered when comparing multi-way and
traditional Hi-C datasets.

### Memory considerations

For large datasets, use `region` parameter to process chromosomes
individually and reduce memory usage.

## See also

[`ChromatinContacts()`](https://jasonwong-lab.github.io/gghic/reference/ChromatinContacts.md),
[`gghic()`](https://jasonwong-lab.github.io/gghic/reference/gghic.md),
[`MultiWayContacts()`](https://jasonwong-lab.github.io/gghic/reference/MultiWayContacts-class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load concatemer data
concatemers <- readRDS("concatemers.rds")

# Convert to pairwise at 100kb resolution
gis <- concatemers2Gis(concatemers, bin_size = 100000)

# Visualize as Hi-C map
gghic(gis)

# Focus on specific region at higher resolution
region <- GRanges("chr1:1000000-5000000")
gis_region <- concatemers2Gis(
  concatemers,
  region = region,
  bin_size = 10000
)

# Process by chromosome to manage memory
chr_gis <- lapply(paste0("chr", 1:22), function(chr) {
  region <- GRanges(seqnames = chr,
                    ranges = IRanges(1, 250000000))
  concatemers2Gis(concatemers, region, bin_size = 100000)
})

# Use with ChromatinContacts
gis_all <- do.call(c, chr_gis)
# Then visualize or analyze
} # }
```
