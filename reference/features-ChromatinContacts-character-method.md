# Get specific genomic features from ChromatinContacts object

Retrieves a specific type of genomic feature (TADs, loops, compartments,
or multi-contacts) from a `ChromatinContacts` object.

## Usage

``` r
# S4 method for class 'ChromatinContacts,character'
features(x, name)
```

## Arguments

- x:

  A `ChromatinContacts` object.

- name:

  Character string. Name of the feature to retrieve. Must be one of:
  `"compartments"`, `"TADs"`, `"loops"`, or `"multi_contacts"`.

## Value

The requested feature as GRanges (for TADs, compartments,
multi_contacts) or GInteractions (for loops).

## Examples

``` r
if (FALSE) { # \dontrun{
# Get TADs
tads <- features(cc, "TADs")

# Get loops
loops <- features(cc, "loops")
} # }
```
