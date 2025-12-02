# Get all genomic features from ChromatinContacts object

Retrieves all genomic features (TADs, loops, compartments,
multi-contacts) as a SimpleList.

## Usage

``` r
# S4 method for class 'ChromatinContacts,missing'
features(x)
```

## Arguments

- x:

  A `ChromatinContacts` object.

## Value

A `SimpleList` containing all features: compartments, TADs, loops, and
multi_contacts.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all features
all_features <- features(cc)
all_features$TADs
all_features$loops
} # }
```
