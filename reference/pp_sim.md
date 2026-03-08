# Generate samples from the prior predictive distribution of model parameters

Generates samples from the prior distributions of the model parameters.
This is primarily used for prior predictive checks or to generate
initial values for MCMC.

## Usage

``` r
pp_sim(dat, model = "drm", selectivity = NULL, ...)
```

## Arguments

- dat:

  A `list` containing the data and prior parameters, typically generated
  by
  [`make_data()`](https://pinskylab.github.io/drmr/reference/make_data.md)
  or
  [`make_data_sdm()`](https://pinskylab.github.io/drmr/reference/make_data_sdm.md).

- model:

  A `character` string specifying the model type. Must be either "drm"
  or "sdm". Defaults to "drm".

- selectivity:

  a numeric `vector` with the same length as the number of age-groups.

- ...:

  parameters passed on to `pop_dyn`.

## Value

A `list` containing samples drawn from the prior distributions of the
model parameters.

## See also

[`prior_inits()`](https://pinskylab.github.io/drmr/reference/prior_inits.md)

## Author

lcgodoy
