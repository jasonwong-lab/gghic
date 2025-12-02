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
