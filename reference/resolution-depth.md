# Calculate Resolution and Depth for Pairs Data

Functions to analyze genomic interaction data resolution and coverage
depth.

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

  A data frame or tibble with columns: `read_name`, `chrom1`, `pos1`,
  `chrom2`, `pos2`. Can also have additional columns like `strand1`,
  `strand2`, etc. For large files, use `pairs_file` parameter instead to
  read in chunks.

- bin_size:

  Integer bin size in base pairs.

- min_contacts:

  Integer minimum number of contacts required for a bin. Default: 1000.

- min_bin:

  Integer minimum bin size for search. Default: 1000.

- max_bin:

  Integer maximum bin size for search. Default: 5000000.

- target_coverage:

  Numeric target coverage fraction (0-1). Default: 0.8.

- pairs_file:

  Character path to pairs file.

- cache:

  External pointer to cached pairs data (from `readPairsCache()`). If
  provided, data is read from cache instead of file for much faster
  repeated operations.

- return_cache:

  Logical, if TRUE, `findOptResChunked()` returns a list with both the
  optimal bin size and the cache pointer. Default: FALSE.

## Value

- `calculateResolutionDepth()`: A tibble with columns `chrom`, `bin`,
  and `count`.

- `calculateGenomeCoverage()`: A numeric value representing the fraction
  of bins with \>= min_contacts.

- `findOptimalResolution()`: An integer representing the optimal bin
  size.

- `calcResDepthChunked()`: Same as above, reads file in chunks.

- `calcGenomeCovChunked()`: Coverage calculation from file.

- `findOptResChunked()`: Optimal resolution for large files.

An external pointer to the cached data.

## Details

These functions analyze Hi-C pairs data to determine optimal genomic
resolution. `calculateResolutionDepth()` bins genomic positions and
counts unique interactions per bin. `calculateGenomeCoverage()` computes
the fraction of bins meeting a contact threshold.
`findOptimalResolution()` uses binary search to find the smallest bin
size achieving target coverage.

For large files that don't fit in memory, use the `_chunked` variants
with `pairs_file` parameter. These read the file in chunks and aggregate
results.

For large datasets, C-accelerated computation is used when available.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load pairs data
pairs <- read.table("contact_matrix.txt",
  header = FALSE,
  col.names = c(
    "read_name", "strand1", "chrom1", "pos1", "frag1",
    "strand2", "chrom2", "pos2", "frag2", "mapq1", "mapq2"
  )
)

# Calculate resolution and depth
depth_10kb <- calculateResolutionDepth(pairs, bin_size = 10000)

# For large file, read in chunks
depth_10kb <- calcResDepthChunked(
  pairs_file = "contact_matrix.txt",
  bin_size = 10000
)

# Calculate coverage
coverage <- calculateGenomeCoverage(pairs, bin_size = 10000)

# Find optimal resolution
opt_bin <- findOptimalResolution(pairs, target_coverage = 0.8)

# For large file - method 1: let function read file each time
opt_bin <- findOptResChunked(
  pairs_file = "contact_matrix.txt",
  target_coverage = 0.8
)

# For large file - method 2: cache once, reuse many times (FASTER!)
cache <- readPairsCache("contact_matrix.txt")

# Find optimal resolution using cache
opt_bin <- findOptResChunked(cache = cache, target_coverage = 0.8)

# Reuse cache for multiple analyses (no file I/O!)
depth_10kb <- calcResDepthChunked(cache = cache, bin_size = 10000)
depth_50kb <- calcResDepthChunked(cache = cache, bin_size = 50000)
coverage <- calcGenomeCovChunked(cache = cache, bin_size = 10000)

# Method 3: Get cache back from findOptResChunked
result <- findOptResChunked(
  pairs_file = "contact_matrix.txt",
  target_coverage = 0.8,
  return_cache = TRUE
)
opt_bin <- result$bin_size
cache <- result$cache

# Now reuse the cache
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
