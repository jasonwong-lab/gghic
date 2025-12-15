# Get genomic features

Retrieves specific or all genomic features.

## Usage

``` r
features(x, name)

# S4 method for class 'ChromatinContacts,character'
features(x, name)

# S4 method for class 'ChromatinContacts,missing'
features(x)
```

## Arguments

- x:

  ChromatinContacts object.

- name:

  Character. Feature name: `"compartments"`, `"TADs"`, `"loops"`, or
  `"multi_contacts"`. If missing, returns all as SimpleList.

## Value

GRanges (TADs, compartments, multi_contacts), GInteractions (loops), or
SimpleList (all features).

## Examples

``` r
if (FALSE) { # \dontrun{
tads <- features(cc, "TADs")
loops <- features(cc, "loops")
all_features <- features(cc)
} # }
```
