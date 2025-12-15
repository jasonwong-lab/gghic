# gghic 0.2.1

## Major Changes

- **S4 Method System for Hypergraphs**: Converted hypergraph construction workflow from S3 to S4 methods for better type safety and extensibility.
  - `buildHypergraph()` is now a generic S4 method `build()` with methods for `MultiWayContacts` objects.
  - `tidyHypergraph()` replaced with S4 generic `tidy()` method.
  - `selectHypergraph()` replaced with S4 generic `select()` method.
  - All methods support method dispatch on `MultiWayContacts` class.
  - Maintains backward compatibility through consistent API design.
  - Added `show()` method for `MultiWayContacts` class providing informative object summaries.

## Documentation Improvements

- **Comprehensive Documentation Overhaul**: Significantly enhanced roxygen2 documentation across all user-facing functions.
  - Added detailed `@details` sections explaining algorithms, workflows, and technical considerations.
  - Expanded `@param` descriptions with valid options, defaults, and usage guidance.
  - Included 8-12 practical examples per major function covering common use cases.
  - Added `@seealso` cross-references linking related functions.

- **Enhanced Function Documentation**:
  - `geom_hic()`: Added coordinate system explanation, transformation details, scaling algorithms.
  - `geom_loop()` and `geom_loop2()`: Detailed BEDPE format specification, styling options.
  - `geom_tad()` and `geom_tad2()`: Comprehensive boundary visualization guidance.
  - `geom_track()` and `geom_track2()`: File format support, scaling strategies, comparison table.
  - `geom_annotation()`: Gene model visualization, GTF parsing, style comparison.
  - `geom_concatemer()` and `geom_concatemer2()`: Multi-way contact representation.
  - `geom_hypergraph()`: Network visualization workflow and layout algorithms.
  - `geom_ideogram()`: Cytogenetic band display, genome assembly support.
  - `geom_hic_under()`: Inverted heatmap visualization for comparative displays.
  - `scaleData()`: Transformation functions (log, asinh, sqrt) with mathematical formulas.
  - All S4 classes and methods: Comprehensive slot descriptions and usage examples.

## Bug Fixes

- Fixed unresolved roxygen2 link in `setters-methods.R` (changed `[import()]` to `[import-ChromatinContacts]`).

## Internal Changes

- Improved consistency in function naming and parameter conventions.
- Enhanced code organization for hypergraph workflow methods.

# gghic 0.2.0

## Major Changes

- **ChromatinContacts S4 Class**: Introduced a new S4 class for managing Hi-C/-like data with improved object-oriented design.
  - New `import()` methods for loading data from Cooler files.
  - Setter methods for assigning genomic features (TADs, loops, tracks, compartments).
  - Comprehensive subsetting methods supporting numeric, logical, character, GRanges, and GInteractions indexing.
  - Enhanced `show()` method for better object summaries.

## New Features

- **Multi-way Contact Visualization**: Added `geom_hypergraph()` for visualizing Pore-C/-like multi-way chromatin interactions.
- **Lower Triangle Plots**: New `geom_hic_under()` for plotting flipped Hi-C heatmaps below the diagonal.
- **Resolution Depth Analysis**: Implemented C-based resolution depth calculation with caching support for improved performance.
- **Pairs File Support**: Added `readPairsCache()` for reading and caching pairs files with filtering capabilities.
- **Concatemer Conversion**: New `concatemers2Gis()` function to convert concatemer reads to GInteractions.

## Performance Improvements

- Implemented C code for resolution depth calculations with hash table-based bin tracking.
- Added comprehensive caching system for resolution and depth calculations.
- Memory-efficient pairs file I/O with filtering support.

## Documentation

- Added comprehensive hypergraph vignette for multi-way contact analysis.
- Enhanced main vignette with improved workflow examples and data loading.
- Refactored depth vignette with better caching explanations.
- Updated pkgdown website with better navigation and feature descriptions.

## Enhancements

- Refactored all `geom_*` functions for consistent parameter handling.
- Added `renameGeomAes()` utility for remapping ggplot2 aesthetics (source: <https://github.com/clauswilke/relayer>).
- Better error messages and validation throughout.

## Internal Changes

- Refactored utility functions with consistent naming conventions.
- Added Makevars with compiler flags for C code compilation.

# gghic 0.1.0

- Added core functionalities for visualising genomic interaction data.
- Integrated with ggplot2 for flexible and customisable plotting.
- Included support for heatmaps, chromosome ideograms, gene tracks, and other genomic data tracks.
