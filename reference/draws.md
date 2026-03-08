# Draws method for `adrm` and `sdm` objects.

Mirrors the behavior of `cmdstanr`'s `$draws` method to retrieve the
samples from the posterior distribution of the models' parameters.

## Usage

``` r
draws(x, ...)

# S3 method for class 'drmrmodels'
draws(
  x,
  variables = NULL,
  inc_warmup = FALSE,
  format = getOption("cmdstanr_draws_format", "draws_array"),
  ...
)
```

## Arguments

- x:

  An object of class `adrm` or `sdm`.

- ...:

  currently ignored.

- variables:

  a string vector with the name of the variables for which we want to
  obtain posteriors samples (defaults to the same variables as the
  `summary` methods).

- inc_warmup:

  a boolean indicating whether the warmup samples should be retrieved as
  well. Defaults to `FALSE`.

- format:

  A string. See cmdstanr [\$draws
  documentation](https://mc-stan.org/cmdstanr/reference/fit-method-draws.html).

## Value

an object of class `"draws"` containing the posterior samples from
specified parameters.

## Author

lcgodoy
