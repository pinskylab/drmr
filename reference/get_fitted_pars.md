# Retrieve parameters needed for predictions

This function identifies the parameters necessary for carrying out
predictions based on the data used to fit a DRM (or SDM).

## Usage

``` r
get_fitted_pars(data_list, model = "drm")
```

## Arguments

- data_list:

  a `list` used as input for model fitting. Typically, the output from
  the
  [make_data](https://pinskylab.github.io/drmr/reference/make_data.md)
  function.

- model:

  a `character` indicating which model predictions are sought for. This
  input admits two possible entries: "drm" (default) or "sdm".

## Value

a `character` vector of labels indicating the parameters necessary for
the predictions.

## Author

lcgodoy
