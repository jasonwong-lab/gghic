# Changelog

## gghic 0.2.1

### Major Changes

- **S4 Method System for Hypergraphs**: Converted hypergraph
  construction workflow from S3 to S4 methods for better type safety and
  extensibility.
  - `buildHypergraph()` is now a generic S4 method
    [`build()`](https://jasonwong-lab.github.io/gghic/reference/build.md)
    with methods for `MultiWayContacts` objects.
  - `tidyHypergraph()` replaced with S4 generic
    [`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md)
    method.
  - `selectHypergraph()` replaced with S4 generic
    [`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md)
    method.
  - All methods support method dispatch on `MultiWayContacts` class.
  - Maintains backward compatibility through consistent API design.
  - Added [`show()`](https://rdrr.io/r/methods/show.html) method for
    `MultiWayContacts` class providing informative object summaries.

### Documentation Improvements

- **Comprehensive Documentation Overhaul**: Significantly enhanced
  roxygen2 documentation across all user-facing functions.
  - Added detailed `@details` sections explaining algorithms, workflows,
    and technical considerations.
  - Expanded `@param` descriptions with valid options, defaults, and
    usage guidance.
  - Included 8-12 practical examples per major function covering common
    use cases.
  - Added `@seealso` cross-references linking related functions.
- **Enhanced Function Documentation**:
  - [`geom_hic()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic.md):
    Added coordinate system explanation, transformation details, scaling
    algorithms.
  - [`geom_loop()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop.md)
    and
    [`geom_loop2()`](https://jasonwong-lab.github.io/gghic/reference/geom_loop2.md):
    Detailed BEDPE format specification, styling options.
  - [`geom_tad()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad.md)
    and
    [`geom_tad2()`](https://jasonwong-lab.github.io/gghic/reference/geom_tad2.md):
    Comprehensive boundary visualization guidance.
  - [`geom_track()`](https://jasonwong-lab.github.io/gghic/reference/geom_track.md)
    and
    [`geom_track2()`](https://jasonwong-lab.github.io/gghic/reference/geom_track2.md):
    File format support, scaling strategies, comparison table.
  - [`geom_annotation()`](https://jasonwong-lab.github.io/gghic/reference/geom_annotation.md):
    Gene model visualization, GTF parsing, style comparison.
  - [`geom_concatemer()`](https://jasonwong-lab.github.io/gghic/reference/geom_concatemer.md)
    and
    [`geom_concatemer2()`](https://jasonwong-lab.github.io/gghic/reference/geom_concatemer2.md):
    Multi-way contact representation.
  - [`geom_hypergraph()`](https://jasonwong-lab.github.io/gghic/reference/geom_hypergraph.md):
    Network visualization workflow and layout algorithms.
  - [`geom_ideogram()`](https://jasonwong-lab.github.io/gghic/reference/geom_ideogram.md):
    Cytogenetic band display, genome assembly support.
  - [`geom_hic_under()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic_under.md):
    Inverted heatmap visualization for comparative displays.
  - [`scaleData()`](https://jasonwong-lab.github.io/gghic/reference/scaleData.md):
    Transformation functions (log, asinh, sqrt) with mathematical
    formulas.
  - All S4 classes and methods: Comprehensive slot descriptions and
    usage examples.

### Bug Fixes

- Fixed unresolved roxygen2 link in `setters-methods.R` (changed
  `[import()]` to `[import-ChromatinContacts]`).

### Internal Changes

- Improved consistency in function naming and parameter conventions.
- Enhanced code organization for hypergraph workflow methods.

## gghic 0.2.0

### Major Changes

- **ChromatinContacts S4 Class**: Introduced a new S4 class for managing
  Hi-C/-like data with improved object-oriented design.
  - New `import()` methods for loading data from Cooler files.
  - Setter methods for assigning genomic features (TADs, loops, tracks,
    compartments).
  - Comprehensive subsetting methods supporting numeric, logical,
    character, GRanges, and GInteractions indexing.
  - Enhanced [`show()`](https://rdrr.io/r/methods/show.html) method for
    better object summaries.

### New Features

- **Multi-way Contact Visualization**: Added
  [`geom_hypergraph()`](https://jasonwong-lab.github.io/gghic/reference/geom_hypergraph.md)
  for visualizing Pore-C/-like multi-way chromatin interactions.
- **Lower Triangle Plots**: New
  [`geom_hic_under()`](https://jasonwong-lab.github.io/gghic/reference/geom_hic_under.md)
  for plotting flipped Hi-C heatmaps below the diagonal.
- **Resolution Depth Analysis**: Implemented C-based resolution depth
  calculation with caching support for improved performance.
- **Pairs File Support**: Added
  [`readPairsCache()`](https://jasonwong-lab.github.io/gghic/reference/resolution-depth.md)
  for reading and caching pairs files with filtering capabilities.
- **Concatemer Conversion**: New
  [`concatemers2Gis()`](https://jasonwong-lab.github.io/gghic/reference/concatemers2gis.md)
  function to convert concatemer reads to GInteractions.

### Performance Improvements

- Implemented C code for resolution depth calculations with hash
  table-based bin tracking.
- Added comprehensive caching system for resolution and depth
  calculations.
- Memory-efficient pairs file I/O with filtering support.

### Documentation

- Added comprehensive hypergraph vignette for multi-way contact
  analysis.
- Enhanced main vignette with improved workflow examples and data
  loading.
- Refactored depth vignette with better caching explanations.
- Updated pkgdown website with better navigation and feature
  descriptions.

### Enhancements

- Refactored all `geom_*` functions for consistent parameter handling.
- Added
  [`renameGeomAes()`](https://jasonwong-lab.github.io/gghic/reference/renameGeomAes.md)
  utility for remapping ggplot2 aesthetics (source:
  <https://github.com/clauswilke/relayer>).
- Better error messages and validation throughout.

### Internal Changes

- Refactored utility functions with consistent naming conventions.
- Added Makevars with compiler flags for C code compilation.

## gghic 0.1.0

- Added core functionalities for visualising genomic interaction data.
- Integrated with ggplot2 for flexible and customisable plotting.
- Included support for heatmaps, chromosome ideograms, gene tracks, and
  other genomic data tracks.
