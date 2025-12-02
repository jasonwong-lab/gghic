# gghic

**gghic** is a powerful and flexible R package for visualizing 3D genome organization from Hi-C and multi-way contact data. Built on `ggplot2`, it provides publication-ready visualizations with minimal code while maintaining full customization capabilities.

## Features

- 🎨 **Intuitive API** - Simple high-level functions with sensible defaults
- 🧬 **Comprehensive** - Visualize Hi-C heatmaps, TADs, loops, genes, and signal tracks
- 🚀 **Memory Efficient** - Focused data loading for large datasets
- 📊 **Publication Ready** - Generate high-quality figures with minimal effort
- 🔧 **Customizable** - Full ggplot2 compatibility for advanced styling
- 🧪 **Multi-way Contacts** - Native support for Pore-C and other multi-contact data

## Installation

### From GitHub

Install the development version using `devtools`:

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("devtools", quietly = TRUE)) install.packages("devtools")

# Install gghic
devtools::install_github("jasonwong-lab/gghic")
```

## Documentation

- [Get started with gghic for Hi-C/-like data visualization](articles/gghic.html)
- [Hi-C/-like data quality assessment](articles/depth.html)
- [Pore-C/-like multi-way contact analysis](articles/hypergraph.html)

## Getting Help

- 🐛 **Bug Reports:** [GitHub Issues](https://github.com/jasonwong-lab/gghic/issues)
- 💡 **Feature Requests:** [GitHub Discussions](https://github.com/jasonwong-lab/gghic/discussions)
- 📧 **Contact:**
  - Open an issue for questions and support
  - Email <mjhk@connect.hku.hk> for direct inquiries
