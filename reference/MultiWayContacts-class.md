# MultiWayContacts S4 class

S4 class for multi-way chromatin contacts from long-read sequencing
(Pore-C, Micro-C XL). Represents interactions as hypergraphs where reads
connect multiple genomic loci.

Constructor for MultiWayContacts from .pairs files.

## Usage

``` r
MultiWayContacts(pairs_path = NULL, focus = NULL)
```

## Arguments

- pairs_path:

  Character. Path to .pairs or .pairs.gz file.

- focus:

  Character or NULL. Chromosome(s) to analyze (NULL for genome-wide).

## Details

Workflow: Create → Import → Build → Tidy → Select → Visualize

Key features:

- Hypergraph representation: reads as hyperedges connecting bins

- Sparse matrix storage for memory efficiency

- Multiple weight normalizations

- Flexible filtering by order, chromosome, and weight

## Slots

- `pairs_path`:

  Character. Path to .pairs file.

- `pairs`:

  Data.frame or NULL. Loaded pairs data.

- `bin_size`:

  Integer or NULL. Genomic bin size in base pairs.

- `focus`:

  Character or NULL. Chromosome(s) to analyze.

- `hypergraph`:

  dgCMatrix or NULL. Sparse incidence matrix (bins × hyperedges).

- `weights`:

  Data.frame. Hyperedge weights with normalizations (raw, log, by_order,
  minmax).

- `bin_info`:

  Data.frame. Genomic bin information.

- `hyperedge_reads`:

  SimpleList. Read names per hyperedge.

- `multiways`:

  Numeric. Contact count per hyperedge.

- `tidied_hypergraph`:

  Data.frame or NULL. Long-format hypergraph.

- `select_hypergraph`:

  Data.frame or NULL. Selected hyperedges for plotting.

## See also

`MultiWayContacts()`,
[import-MultiWayContacts](https://jasonwong-lab.github.io/gghic/reference/import-MultiWayContacts.md),
[`build()`](https://jasonwong-lab.github.io/gghic/reference/build.md),
[`tidy()`](https://jasonwong-lab.github.io/gghic/reference/tidy.md),
[`select()`](https://jasonwong-lab.github.io/gghic/reference/select.md),
[`gghypergraph()`](https://jasonwong-lab.github.io/gghic/reference/gghypergraph.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
mc <- MultiWayContacts("sample.pairs.gz")

# Focus on chromosome
mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1")

# Complete workflow
mc <- MultiWayContacts("sample.pairs.gz", focus = "chr1") |>
  import() |>
  build(bin_size = 1000000L) |>
  tidy() |>
  select(n_intra = 10, n_inter = 5) |>
  gghypergraph()
} # }
```
