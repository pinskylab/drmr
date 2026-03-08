# Generate a random sample from a model's predictive distribution given a set of parameters

Generates samples from the predictive distribution of a model given a
set of parameters. This is primarily used for prior predictive checks
and simulation studies.

## Usage

``` r
model_sim(dat, model, selectivity, pars, ...)
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

- pars:

  a named `list` of model parameters.

- ...:

  parameters passed on to `pop_dyn`.

## Value

A `list` containing samples drawn from the prior distributions of the
model parameters.

## Details

See [this
link](https://forum.posit.co/t/how-to-solve-no-visible-binding-for-global-variable-note/28887/3)
for the [rlang::.data](https://rlang.r-lib.org/reference/dot-data.html)
import.

## See also

[`prior_inits()`](https://pinskylab.github.io/drmr/reference/prior_inits.md)

## Author

lcgodoy
