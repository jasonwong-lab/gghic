# Calculate Resolution and Depth for Pairs Data

Functions to analyze genomic interaction resolution and coverage depth.

Load a pairs file into C memory for fast repeated analysis. The cache
persists until garbage collected or explicitly cleared.

## Usage

``` r
calculateResolutionDepth(pairs, bin_size)

calculateGenomeCoverage(pairs, bin_size, min_contacts = 1000)

findOptimalResolution(
  pairs,
  min_bin = 1000,
  max_bin = 5e+06,
  target_coverage = 0.8,
  min_contacts = 1000
)

readPairsCache(pairs_file)

calcResDepthChunked(pairs_file = NULL, bin_size, cache = NULL)

calcGenomeCovChunked(
  pairs_file = NULL,
  bin_size,
  min_contacts = 1000,
  cache = NULL
)

findOptResChunked(
  pairs_file = NULL,
  min_bin = 1000,
  max_bin = 5e+06,
  target_coverage = 0.8,
  min_contacts = 1000,
  cache = NULL,
  return_cache = FALSE
)
```

## Arguments

- pairs:

  Data frame with columns: `read_name`, `chrom1`, `pos1`, `chrom2`,
  `pos2`. For large files, use `pairs_file` instead.

- bin_size:

  Integer. Bin size in base pairs.

- min_contacts:

  Integer. Minimum contacts per bin (default: 1000).

- min_bin:

  Integer. Minimum bin size for search (default: 1000).

- max_bin:

  Integer. Maximum bin size for search (default: 5000000).

- target_coverage:

  Numeric (0-1). Target coverage fraction (default: 0.8).

- pairs_file:

  Character path to pairs file.

- cache:

  External pointer from `readPairsCache()`. Faster for repeated
  operations.

- return_cache:

  Logical. Return list with bin size and cache pointer (default: FALSE).

## Value

- `calculateResolutionDepth()`: Tibble with chrom, bin, count.

- `calculateGenomeCoverage()`: Coverage fraction (bins ≥ min_contacts).

- `findOptimalResolution()`: Optimal bin size (integer).

- `calcResDepthChunked()`: Same, reads file in chunks.

- `calcGenomeCovChunked()`: Coverage from file.

- `findOptResChunked()`: Optimal resolution for large files.

An external pointer to the cached data.

## Details

Analyzes Hi-C pairs to determine optimal resolution. For large files,
use `_chunked` variants or cache with `readPairsCache()` for speed.

## Examples

``` r
if (FALSE) { # \dontrun{
pairs <- read.table("contact_matrix.txt")
depth_10kb <- calculateResolutionDepth(pairs, bin_size = 10000)
coverage <- calculateGenomeCoverage(pairs, bin_size = 10000)
opt_bin <- findOptimalResolution(pairs, target_coverage = 0.8)

# Large file with cache
cache <- readPairsCache("contact_matrix.txt")
opt_bin <- findOptResChunked(cache = cache, target_coverage = 0.8)
depth <- calcResDepthChunked(cache = cache, bin_size = opt_bin)
} # }

if (FALSE) { # \dontrun{
# Create cache once
cache <- readPairsCache("data.pairs.gz")

# Reuse cache for multiple operations
depth1 <- calcResDepthChunked(cache = cache, bin_size = 10000)
depth2 <- calcResDepthChunked(cache = cache, bin_size = 50000)
coverage <- calcGenomeCovChunked(cache = cache, bin_size = 10000)

# Cache is automatically freed when R session ends
# Or explicitly remove it:
rm(cache)
gc()
} # }
```
