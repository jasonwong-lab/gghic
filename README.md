<!-- markdownlint-configure-file {
  "no-inline-html": {
    "allowed_elements": [
      "p", "a", "img"
    ]
  }
} -->

# gghic <img src="https://raw.githubusercontent.com/jasonwong-lab/gghic/refs/heads/master/man/figures/logo.png" alt="gghic's logo" height="180" align="right" />

<p align="center"> A versatile R package for exploring and visualizing 3D genome organization </p>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![GitHub issues](https://img.shields.io/github/issues/jasonwong-lab/gghic)](https://github.com/jasonwong-lab/gghic/issues) [![GitHub stars](https://img.shields.io/github/stars/jasonwong-lab/gghic)](https://github.com/jasonwong-lab/gghic/stargazers)

> :warning: `gghic` is in its early stages. Issues and PRs are welcome.

## Introduction

**gghic** is a versatile R package for creating flexible, publication-ready visualizations of 3D genome organization data. With gghic, you can easily explore and present Hi-C/-like contact maps, topologically associating domains (TADs), chromatin loops, gene annotations, and other genomic features in a unified and customizable framework.

Whether you are analyzing large-scale Hi-C experiments, visualizing multi-chromosome interactions, or integrating additional genomic tracks (such as ChIP-seq or BigWig data), gghic provides a tidyverse-friendly and extensible toolkit for your research.

## Features

1. [x] Introduce new S4 classes `ChromatinContacts` and `MultiWayContacts` for managing pair-wise and multi-way contact data.
2. [x] Plot triangular heatmaps for genomic interaction data containing multiple chromosomes.
3. [x] Plot chromosome ideograms with highlighted regions above the heatmap.
4. [x] Plot gene/transcript model tracks under the heatmap.
5. [x] Plot tracks of other genomic data, *e.g.*, ChIP-seq peaks, BigWig files, concatemers, *etc.*, under the heatmap.
6. [x] Plot compartments, TADs, and loops on the heatmap.
7. [x] Allow for generating rasterized heatmaps and tracks to handle large datasets.
8. [x] Plot concatemers indicating multi-way contacts.
9. [x] Use concatemers only to plot pairwise interactions.
10. [x] Plot a second heatmap below the main heatmap (lower triangle).
11. [x] Check and suggest optimal resolutions based on sequencing depth.
12. [x] Build hypergraph representations from multi-way contact data.
13. [ ] Support duckplyr for faster data manipulation.
14. [ ] ...

## Installation

To install the `gghic` package from GitHub, use the following command:

```r
devtools::install_github("jasonwong-lab/gghic")
```

## Usage

- 📖 **Documentation:**
  Visit the <a href="https://jasonwong-lab.github.io/gghic/" target="_blank">gghic website</a> for comprehensive documentation, tutorials, and examples.

## Citation

```r
citation("gghic")
```

## Authors

Minghao Jiang, <mjhk@connect.hku.hk> [![gghic contributors](https://contrib.rocks/image?repo=jasonwong-lab/gghic)](https://github.com/jasonwong-lab/gghic/graphs/contributors)

## License

**`gghic`** is licensed under the [GNU General Public License v3](LICENSE.md).

---

[![Star History Chart](https://api.star-history.com/svg?repos=jasonwong-lab/gghic&type=Date)](https://star-history.com/#jasonwong-lab/gghic&Date)
