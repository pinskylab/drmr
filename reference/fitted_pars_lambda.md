# Retrieve parameters needed for predictions

This function identifies the parameters necessary for carrying out
predictions based on the data used to fit a DRM (or SDM).

## Usage

``` r
fitted_pars_lambda(data_list)
```

## Arguments

- data_list:

  a `list` used as input for model fitting. Typically, the output from
  the
  [make_data](https://pinskylab.github.io/drmr/reference/make_data.md)
  function.

## Value

a `character` vector of labels indicating the parameters necessary for
the predictions.

## Author

lcgodoy
