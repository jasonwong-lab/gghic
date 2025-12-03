# Resolution Depth Analysis for Hi-C/-like Data

## Overview

This vignette demonstrates how to analyze Hi-C pairs data to find
optimal genomic resolutions using gghic’s resolution depth functions.
You’ll learn how to:

1.  **Calculate contact depth** at different bin sizes
2.  **Assess genome coverage** across resolutions
3.  **Find optimal resolution** that balances detail and completeness
4.  **Visualize tradeoffs** between resolution and coverage

### When to Use These Functions

These functions are essential for:

- **Pore-C & long-read data**: Determining the best resolution for long
  DNA-DNA interactions
- **Hi-C QC**: Assessing data quality and sequencing depth
- **Resource planning**: Finding the sweet spot before expensive
  downstream analysis
- **Publication standards**: Communicating data completeness to
  reviewers

## Getting Started

Load required libraries:

``` r
load_pkg <- function(pkgs) {
  for (pkg in pkgs) suppressMessages(require(pkg, character.only = TRUE))
}

load_pkg(c("readr", "dplyr", "rappdirs", "gghic"))
#> Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
#> logical.return = TRUE, : there is no package called 'readr'
```

### Example Data

We’ll create synthetic pairs data for demonstration. In practice, load
from your contact matrix file:

``` r
# Create example pairs dataset
set.seed(42)
n_interactions <- 10000

pairs <- data.frame(
  read_name = paste0("read_", seq_len(n_interactions)),
  chrom1 = sample(c("chr1", "chr2", "chr3"), n_interactions, replace = TRUE),
  pos1 = sample(1e6:10e6, n_interactions, replace = TRUE),
  chrom2 = sample(c("chr1", "chr2", "chr3"), n_interactions, replace = TRUE),
  pos2 = sample(1e6:10e6, n_interactions, replace = TRUE)
)

glimpse(pairs)
#> Rows: 10,000
#> Columns: 5
#> $ read_name <chr> "read_1", "read_2", "read_3", "read_4", "read_5", "read_6", …
#> $ chrom1    <chr> "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr…
#> $ pos1      <int> 7501652, 9834771, 7511178, 9988046, 2376108, 8056778, 991689…
#> $ chrom2    <chr> "chr3", "chr1", "chr1", "chr2", "chr1", "chr2", "chr3", "chr…
#> $ pos2      <int> 4068533, 1789229, 6896113, 6916211, 7349759, 9964631, 150451…
```

**Expected column format**:

- `read_name`: Unique identifier for each interaction
- `chrom1`, `chrom2`: Chromosome names (e.g., “chr1”)
- `pos1`, `pos2`: Genomic positions in base pairs
- Optional: strand, mapping quality, fragment ID (ignored by analysis
  functions)

## Choosing Your Approach

The gghic package provides two complementary approaches for analyzing
pairs data:

| Approach                                                                                                          | Memory          | Speed   | When to Use                            |
|-------------------------------------------------------------------------------------------------------------------|-----------------|---------|----------------------------------------|
| **In-memory** [`calculateResolutionDepth()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md) | O(interactions) | Fastest | Files \< 50% RAM, interactive analysis |
| **Chunked** [`calcResDepthChunked()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md)        | O(unique bins)  | Fast    | Files \> 50% RAM, automated pipelines  |

**In-Memory Processing (Smaller Files)** - Load entire dataset into
memory - Single file read (fastest) - Typical memory: ~4 bytes per
interaction - Best for: Interactive exploration, visualization

**Chunked Streaming (Large Files)** - Read file in streaming fashion -
Minimal memory footprint (~100 MB even for 100 GB files) - C-accelerated
file I/O - Automatic gzip decompression - Best for: Production
pipelines, server environments, very large datasets

**Decision Rule**: If your file size \> available RAM / 2, use chunked
functions.

## Step 1: Calculate Resolution Depth

The
[`calculateResolutionDepth()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md)
function bins genomic positions and counts unique interactions per bin:

``` r
# Calculate depth at 10kb resolution
depth_10kb <- calculateResolutionDepth(pairs, bin_size = 10000)

depth_10kb
#> # A tibble: 2,699 × 3
#>    chrom   bin count
#>    <chr> <int> <int>
#>  1 chr1    845    12
#>  2 chr2    812     6
#>  3 chr3    779     7
#>  4 chr1    846     7
#>  5 chr2    813     6
#>  6 chr3    780     4
#>  7 chr1    847     7
#>  8 chr3    781     8
#>  9 chr2    814     6
#> 10 chr1    848     9
#> # ℹ 2,689 more rows
```

``` r
# Summary statistics
summary(depth_10kb$count)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    1.00    5.00    7.00    7.41    9.00   19.00
```

**What this returns**:

- `chrom`: Chromosome name
- `bin`: Bin number (position / bin_size, rounded up)
- `count`: Number of unique interactions in this bin

**Interpreting the results**:

- High counts = well-sequenced regions
- Low counts = sparse or unmappable regions
- Zero-count bins = gaps (excluded from output)

## Step 2: Assess Genome Coverage

The
[`calculateGenomeCoverage()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md)
function determines what fraction of bins meet a minimum contact
threshold:

``` r
# Calculate coverage at multiple resolutions
resolutions <- c(5e3, 10e3, 25e3, 50e3, 100e3, 500e3)

coverage_df <- data.frame(
  resolution_kb = resolutions / 1e3,
  n_bins = NA,
  coverage = NA
)

for (i in seq_along(resolutions)) {
  depth <- calculateResolutionDepth(pairs, resolutions[i])
  coverage_df$n_bins[i] <- nrow(depth)
  coverage_df$coverage[i] <- calculateGenomeCoverage(
    pairs, resolutions[i],
    min_contacts = 50
  )
}

coverage_df$coverage_pct <- sprintf("%.1f%%", coverage_df$coverage * 100)

coverage_df
#>   resolution_kb n_bins   coverage coverage_pct
#> 1             5   5273 0.00000000         0.0%
#> 2            10   2699 0.00000000         0.0%
#> 3            25   1080 0.00000000         0.0%
#> 4            50    540 0.02777778         2.8%
#> 5           100    270 0.99629630        99.6%
#> 6           500     54 1.00000000       100.0%
```

**Interpreting coverage**:

- **High coverage (\>80%)** at small bins: Well-sequenced, good
  resolution potential
- **Low coverage (\<50%)** at large bins: Sparse data, may need deeper
  sequencing
- **Coverage curve shape**: Shows sequencing depth and data quality

## Step 3: Find Optimal Resolution

Use binary search to find the **smallest bin size** achieving your
target coverage:

``` r
# Find resolution with ~80% coverage
optimal_bin <- findOptimalResolution(
  pairs,
  min_bin = 1e3, # Start at 1 Kb
  max_bin = 1e6, # Search up to 1 Mb
  target_coverage = 0.8,
  min_contacts = 50
)
#> Iteration 1: Testing bin size 500,500 bp: 98.18% coverage
#> Iteration 2: Testing bin size 250,749 bp: 97.30% coverage
#> Iteration 3: Testing bin size 125,874 bp: 97.72% coverage
#> Iteration 4: Testing bin size 63,436 bp: 36.60% coverage
#> Iteration 5: Testing bin size 94,655 bp: 97.92% coverage
#> Iteration 6: Testing bin size 79,045 bp: 88.41% coverage
#> Iteration 7: Testing bin size 71,240 bp: 64.83% coverage
#> Iteration 8: Testing bin size 75,142 bp: 77.13% coverage
#> Iteration 9: Testing bin size 77,093 bp: 83.85% coverage
#> Iteration 10: Testing bin size 76,117 bp: 81.23% coverage
#> Iteration 11: Testing bin size 75,629 bp: 78.61% coverage
#> Iteration 12: Testing bin size 75,873 bp: 81.23% coverage
#> Iteration 13: Testing bin size 75,751 bp: 79.55% coverage
#> Iteration 14: Testing bin size 75,812 bp: 79.27% coverage
#> Iteration 15: Testing bin size 75,842 bp: 80.67% coverage
#> Iteration 16: Testing bin size 75,827 bp: 80.39% coverage
#> Iteration 17: Testing bin size 75,819 bp: 79.83% coverage
#> Iteration 18: Testing bin size 75,823 bp: 80.11% coverage
#> Iteration 19: Testing bin size 75,821 bp: 79.55% coverage
#> Iteration 20: Testing bin size 75,822 bp: 79.55% coverage

cat(sprintf(
  "Optimal bin size: %d bp (%.1f Kb)\n", optimal_bin, optimal_bin / 1000
))
#> Optimal bin size: 75823 bp (75.8 Kb)

# Verify the result
actual_coverage <- calculateGenomeCoverage(
  pairs, optimal_bin,
  min_contacts = 50
)
cat(sprintf("Actual coverage: %.2f%%\n", actual_coverage * 100))
#> Actual coverage: 80.11%
```

**How it works**:

1.  Binary search across the bin size range
2.  Tests each candidate bin size
3.  Returns smallest bin that meets target coverage
4.  Typically completes in 10-15 iterations

**Pro tip**: Start with `target_coverage = 0.8` (80%). Adjust based on
your analysis needs.

## Step 4: Visualize Coverage Distribution

Create a boxplot showing how contact counts vary across resolutions:

``` r
plotResolutionCoverage(
  pairs,
  bin_sizes = c(5e3, 10e3, 25e3, 50e3, 100e3, 500e3),
  min_contacts = 50,
  title = "Contact Distribution by Resolution"
)
```

![Resolution depth analysis
plot](depth_files/figure-html/unnamed-chunk-7-1.png)

Distribution of contacts per bin at different resolutions

**What to look for**:

- **Widening boxes**: Decreasing contacts with larger bins (expected)
- **Outliers (dots)**: Unusually enriched regions (e.g., centromeres,
  pericentromeres)
- **Median position**: Typical coverage at that resolution
- **Whisker length**: Range of contact values (high whiskers = variable
  coverage)

## Step 5: Plot Coverage vs Resolution

Create a line plot showing how coverage changes with resolution:

``` r
plotCoverageCurve(
  pairs,
  bin_sizes = c(5e3, 10e3, 25e3, 50e3, 100e3, 500e3, 1e6),
  min_contacts = 50,
  title = "Coverage vs Resolution"
)
```

![Resolution depth analysis
plot](depth_files/figure-html/unnamed-chunk-8-1.png)

Genome coverage as a function of bin size

**Curve interpretation**:

- **Sharp drop at small bins**: Limited sequencing depth
- **Plateau region**: Where additional bin size doesn’t improve coverage
  much
- **“Knee” of curve**: The optimal resolution sweet spot
- **Steep regions**: Good resolution potential at that bin size

## Processing Large Files (100 GB+)

> **⚠️ Note**: The code examples in this section are provided for
> reference but are **not executed** in this vignette.

For files that don’t fit in R memory, use the `_chunked` variants which
stream the file with C-level I/O:

``` r
# Download example file with caching
cache_dir <- rappdirs::user_cache_dir("gghic")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

pairs_file <- file.path(cache_dir, "test.txt.gz")
download_url <- "https://www.dropbox.com/scl/fi/yc4axg1mf2i9zylg3d0oe/test.txt.gz?rlkey=sdsdhsnixo01koo38d242y4c8&st=t4rn0js0&dl=1"

if (!file.exists(pairs_file)) {
  message("Downloading test data to cache directory...")
  download.file(
    download_url, pairs_file, method = "auto", quiet = TRUE
  )
  message("Downloaded to: ", pairs_file)
} else {
  message("Using cached file: ", pairs_file)
}
```

### Approach 1: Direct File Reading (Simple but Slower)

Read the file each time a function is called:

``` r
# Each function reads the file independently
opt_bin <- findOptResChunked(
  pairs_file = pairs_file,
  min_bin = 1e3,
  max_bin = 5e4,
  target_coverage = 0.8,
  min_contacts = 1000
)

# Reads file again
coverage <- calcGenomeCovChunked(
  pairs_file = pairs_file,
  bin_size = opt_bin,
  min_contacts = 1000
)

# Reads file again
depth <- calcResDepthChunked(
  pairs_file = pairs_file,
  bin_size = opt_bin
)
```

**Problem**: For a 10 GB file, you’re reading 10 GB three times (30 GB
total I/O)!

### Approach 2: Cache Once, Reuse Many Times (FASTEST!)

Read the file once into C memory, then reuse the cached data:

``` r
# 1. READ FILE INTO C MEMORY ONCE ----
# This takes time initially but creates a persistent cache
cache <- readPairsCache(pairs_file)
cat("File loaded into C memory cache\n")

# 2. FIND OPTIMAL RESOLUTION (using cache) ----
# This is now instant - no file I/O!
opt_bin <- findOptResChunked(
  cache = cache,
  min_bin = 1e3,
  max_bin = 5e4,
  target_coverage = 0.8,
  min_contacts = 1000
)
cat("Optimal resolution:", opt_bin / 1e3, "Kb\n")

# 3. CHECK COVERAGE (using cache) ----
# Also instant!
coverage <- calcGenomeCovChunked(
  cache = cache,
  bin_size = opt_bin,
  min_contacts = 1000
)
cat("Coverage:", sprintf("%.1f%%\n", coverage * 100))

# 4. EXTRACT BINNED DEPTHS (using cache) ----
# Still instant!
depth <- calcResDepthChunked(
  cache = cache,
  bin_size = opt_bin
)

# 5. TRY DIFFERENT BIN SIZES (no additional I/O!) ----
depth_5kb <- calcResDepthChunked(cache = cache, bin_size = 5000)
depth_10kb <- calcResDepthChunked(cache = cache, bin_size = 10000)
depth_50kb <- calcResDepthChunked(cache = cache, bin_size = 50000)

# Cache is automatically freed when removed or R session ends
rm(cache)
gc()
```

**Performance**: For a 10 GB file: - Without cache: ~30 seconds × 3
operations = **90 seconds** - With cache: ~30 seconds (once) + ~0.1
seconds × operations = **30 seconds**

**3x faster!** And the speedup increases with more operations.

### Approach 3: Get Cache from findOptResChunked

Let
[`findOptResChunked()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md)
return the cache for subsequent use:

``` r
# Find optimal resolution and get cache back
result <- findOptResChunked(
  pairs_file = pairs_file,
  min_bin = 1e3,
  max_bin = 5e4,
  target_coverage = 0.8,
  min_contacts = 1000,
  return_cache = TRUE # Get the cache back!
)

# Extract results
opt_bin <- result$bin_size
cache <- result$cache

cat("Optimal resolution:", opt_bin / 1e3, "Kb\n")

# Now reuse the cache for other operations
coverage <- calcGenomeCovChunked(
  cache = cache,
  bin_size = opt_bin,
  min_contacts = 1000
)

depth <- calcResDepthChunked(
  cache = cache,
  bin_size = opt_bin
)
```

### When to Use Each Approach

| Approach           | Best For                                                                                                        | I/O Operations |
|--------------------|-----------------------------------------------------------------------------------------------------------------|----------------|
| **Direct File**    | Single analysis, simple scripts                                                                                 | N operations   |
| **Explicit Cache** | Interactive exploration, multiple analyses                                                                      | 1 operation    |
| **Return Cache**   | When starting with [`findOptResChunked()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md) | 1 operation    |

**Rule of thumb**: If you’ll call more than one chunked function, use
the cache!

## Interpreting Your Results

### What the Coverage Curve Tells You

**Sharp drop at small bins** - Data is sparse at high resolution -
Limited sequencing depth or short reads - Consider deeper sequencing or
larger bins

**Plateau region** - Coverage maxes out (can’t improve by increasing bin
size) - Indicates natural sequencing depth ceiling - Optimal resolution
likely at transition point

**High coverage at large bins (\>80%)** - Well-sequenced dataset - Can
use smaller bins safely - Good candidate for publication-quality
analysis

**Low coverage at large bins (\<50%)** - Sparse data - May need
additional sequencing - Consider quality filtering

### Choosing Target Coverage

- **Exploration**: 50-60% (quick, uses fewer bins)
- **Standard analysis**: 75-85% (good balance)
- **Publication**: 85-90%+ (most stringent)

Higher coverage = higher confidence but requires more sequencing.

## Next Steps

Once you’ve found your optimal resolution, explore other gghic functions
for visualization:

- [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md):
  Create interaction scatter plots
- [`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md):
  Highlight chromatin loops
- [`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md):
  Show topologically associating domains
- [`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md):
  Add chromosome context
- [`theme_hic()`](https://jasonwong-lab.github.io/gghic/reference/theme_hic.md):
  Specialized theme for Hi-C plots

See the main gghic vignette for complete examples.

## FAQ

**Q: Why does coverage decrease with larger bins?** A: By
definition—larger bins aggregate more positions but contact counts don’t
increase proportionally. Coverage measures the fraction of bins with
sufficient contacts.

**Q: Should I use chunked functions for small files?** A: No, use
in-memory functions—they’re simpler and faster for files that fit in
RAM.

**Q: Can I use these functions on non-Hi-C data?** A: Yes, any “pairs”
data with (chrom, position) coordinates works: Pore-C, 4C, 5C, etc.

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gghic_0.2.0    rappdirs_0.3.3 dplyr_1.1.4   
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1          
#>   [3] jsonlite_2.0.0              magrittr_2.0.4             
#>   [5] GenomicFeatures_1.62.0      farver_2.1.2               
#>   [7] rmarkdown_2.30              fs_1.6.6                   
#>   [9] BiocIO_1.20.0               ragg_1.5.0                 
#>  [11] vctrs_0.6.5                 memoise_2.0.1              
#>  [13] Rsamtools_2.26.0            RCurl_1.98-1.17            
#>  [15] base64enc_0.1-3             htmltools_0.5.8.1          
#>  [17] S4Arrays_1.10.0             progress_1.2.3             
#>  [19] curl_7.0.0                  Rhdf5lib_1.32.0            
#>  [21] SparseArray_1.10.3          Formula_1.2-5              
#>  [23] rhdf5_2.54.0                sass_0.4.10                
#>  [25] bslib_0.9.0                 htmlwidgets_1.6.4          
#>  [27] desc_1.4.3                  Gviz_1.54.0                
#>  [29] httr2_1.2.1                 cachem_1.1.0               
#>  [31] GenomicAlignments_1.46.0    lifecycle_1.0.4            
#>  [33] pkgconfig_2.0.3             Matrix_1.7-4               
#>  [35] R6_2.6.1                    fastmap_1.2.0              
#>  [37] MatrixGenerics_1.22.0       digest_0.6.39              
#>  [39] colorspace_2.1-2            AnnotationDbi_1.72.0       
#>  [41] S4Vectors_0.48.0            textshaping_1.0.4          
#>  [43] Hmisc_5.2-4                 GenomicRanges_1.62.0       
#>  [45] RSQLite_2.4.5               labeling_0.4.3             
#>  [47] filelock_1.0.3              httr_1.4.7                 
#>  [49] abind_1.4-8                 compiler_4.5.2             
#>  [51] withr_3.0.2                 bit64_4.6.0-1              
#>  [53] htmlTable_2.4.3             S7_0.2.1                   
#>  [55] backports_1.5.0             BiocParallel_1.44.0        
#>  [57] DBI_1.2.3                   biomaRt_2.66.0             
#>  [59] DelayedArray_0.36.0         rjson_0.2.23               
#>  [61] tools_4.5.2                 foreign_0.8-90             
#>  [63] nnet_7.3-20                 glue_1.8.0                 
#>  [65] restfulr_0.0.16             InteractionSet_1.38.0      
#>  [67] rhdf5filters_1.22.0         grid_4.5.2                 
#>  [69] checkmate_2.3.3             cluster_2.1.8.1            
#>  [71] generics_0.1.4              gtable_0.3.6               
#>  [73] BSgenome_1.78.0             tidyr_1.3.1                
#>  [75] ensembldb_2.34.0            data.table_1.17.8          
#>  [77] hms_1.1.4                   utf8_1.2.6                 
#>  [79] XVector_0.50.0              BiocGenerics_0.56.0        
#>  [81] pillar_1.11.1               stringr_1.6.0              
#>  [83] BiocFileCache_3.0.0         lattice_0.22-7             
#>  [85] rtracklayer_1.70.0          bit_4.6.0                  
#>  [87] deldir_2.0-4                biovizBase_1.58.0          
#>  [89] tidyselect_1.2.1            Biostrings_2.78.0          
#>  [91] knitr_1.50                  gridExtra_2.3              
#>  [93] IRanges_2.44.0              Seqinfo_1.0.0              
#>  [95] ProtGenerics_1.42.0         SummarizedExperiment_1.40.0
#>  [97] stats4_4.5.2                xfun_0.54                  
#>  [99] Biobase_2.70.0              matrixStats_1.5.0          
#> [101] stringi_1.8.7               UCSC.utils_1.6.0           
#> [103] lazyeval_0.2.2              yaml_2.3.11                
#> [105] evaluate_1.0.5              codetools_0.2-20           
#> [107] cigarillo_1.0.0             interp_1.1-6               
#> [109] tibble_3.3.0                cli_3.6.5                  
#> [111] rpart_4.1.24                systemfonts_1.3.1          
#> [113] jquerylib_0.1.4             dichromat_2.0-0.1          
#> [115] Rcpp_1.1.0                  GenomeInfoDb_1.46.1        
#> [117] dbplyr_2.5.1                png_0.1-8                  
#> [119] XML_3.99-0.20               parallel_4.5.2             
#> [121] pkgdown_2.2.0               ggplot2_4.0.1              
#> [123] blob_1.2.4                  prettyunits_1.2.0          
#> [125] latticeExtra_0.6-31         jpeg_0.1-11                
#> [127] AnnotationFilter_1.34.0     bitops_1.0-9               
#> [129] txdbmaker_1.6.0             VariantAnnotation_1.56.0   
#> [131] scales_1.4.0                purrr_1.2.0                
#> [133] crayon_1.5.3                rlang_1.1.6                
#> [135] KEGGREST_1.50.0
```
