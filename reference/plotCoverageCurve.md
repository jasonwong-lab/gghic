# Plot Coverage Fraction by Resolution

Create a line plot showing how the fraction of bins meeting contact
threshold changes with resolution.

## Usage

``` r
plotCoverageCurve(
  pairs,
  bin_sizes = c(1000, 2000, 5000, 10000, 20000, 25000, 50000, 1e+05, 150000, 2e+05,
    250000, 5e+05, 1e+06),
  min_contacts = 1000,
  title = "Genome Coverage by Resolution"
)
```

## Arguments

- pairs:

  A data frame or tibble with pairs data.

- bin_sizes:

  Integer vector of bin sizes to test. See
  [plotResolutionCoverage](https://jasonwong-lab.github.io/gghic/reference/plotResolutionCoverage.md).

- min_contacts:

  Integer minimum contacts threshold. Default: 1000.

- title:

  Character plot title.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plotCoverageCurve(pairs)
print(p)
} # }
```
