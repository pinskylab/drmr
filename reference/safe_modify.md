# Modifying a named list

Safely modifies a named `list`.

## Usage

``` r
safe_modify(original, replacements)
```

## Arguments

- original:

  a named `list` with "original" parameters.

- replacements:

  a named `list` containing elements to be modified in the `original`
  `list`.

## Value

The updated `original` list.

## Details

This function returns an error if any of the names found in the
`replacements` object are not in the `list` to be modified.

## Author

lcgodoy
