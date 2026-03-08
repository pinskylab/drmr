# Generate initial values for MCMC from the prior

Generates initial values for Markov Chain Monte Carlo (MCMC) by sampling
from the prior distributions of the model parameters. This helps in
starting the MCMC chains from different points in the parameter space.

## Usage

``` r
prior_inits(dat, chains, model = "drm")
```

## Arguments

- dat:

  A list containing the data and prior parameters, typically generated
  by
  [`make_data()`](https://pinskylab.github.io/drmr/reference/make_data.md)
  or
  [`make_data_sdm()`](https://pinskylab.github.io/drmr/reference/make_data_sdm.md).

- chains:

  An integer specifying the number of MCMC chains to initialize.
  Defaults to 4.

- model:

  A `character` string specifying the model type. Must be either "drm"
  or "sdm". Defaults to "drm".

## Value

A `list` of lists, where each inner list contains initial values for one
MCMC chain. The structure of each inner list is the same as the output
of
[`prior_sample()`](https://pinskylab.github.io/drmr/reference/prior_sample.md).

## See also

[`prior_sample()`](https://pinskylab.github.io/drmr/reference/prior_sample.md)

## Author

lcgodoy
