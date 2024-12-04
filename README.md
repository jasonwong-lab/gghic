<!-- markdownlint-configure-file {
  "no-inline-html": {
    "allowed_elements": [
      "p"
    ]
  }
} -->

# gghic

<p align="center"> A versatile R package for exploring and visualizing 3D genome organization </p>

> :warning: `gghic` is in its early stages. Issues and PRs are welcome.

## Author

Minghao Jiang, <mjhk@connect.hku.hk>

## Features

1. [x] Plot triangular heatmaps for genomic interaction data containing multiple chromosomes.
2. [x] Plot chromosome ideograms with highlighted regions above the heatmap.
3. [x] Plot gene/transcript model tracks under the heatmap.
4. [x] Plot tracks of other genomic data, *e.g.*, ChIP-seq peaks, BigWig files, *etc.*, under the heatmap.
5. [x] Plot compartments, TADs, and loops on the heatmap.
6. [x] Allow for generating rasterized heatmaps to handle large datasets.
7. [ ] Support duckplyr for faster data manipulation.
8. [ ] ...

## Installation

To install the `gghic` package from GitHub, use the following command:

```r
devtools::install_github("jasonwong-lab/gghic", build_vignettes = TRUE)
```

## Usage

You can find the vignettes for `gghic` by running the following command:

```r
browseVignettes("gghic")
```

## License

**`gghic`** is licensed under the [GNU General Public License v3](LICENSE.md).
