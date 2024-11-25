<!-- markdownlint-configure-file {
  "no-inline-html": {
    "allowed_elements": [
      "font"
    ]
  }
} -->

# gghic

<font size="2"> *An R package for easy visualization of genomic interaction data* </font>

> :warning: This R package is in its early stages. Issues and PRs are welcome.

## Author

Minghao Jiang, <mjhk@connect.hku.hk>

## Features

1. [x] Visualise genomic interaction data of multiple chromosomes.
2. [x] Plot chromosome ideograms with highlighted regions above the heatmap.
3. [x] Plot gene tracks under the heatmap easily with a GTF file.
4. [x] Plot tracks of other genomic data, *e.g.*, ChIP-seq peaks, BigWig files, *etc.*, under the heatmap.
5. [x] Plot compartments, TADs, and loops on the heatmap.
6. [ ] Support duckplyr for data manipulation.
7. [ ] ...

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
