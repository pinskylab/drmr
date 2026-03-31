# Draws method for `adrm` and `sdm` objects.

Mirrors the behavior of `cmdstanr`'s `$draws` method to retrieve the
samples from the posterior distribution of the models' parameters.

## Usage

``` r
mcmc_diag(x, ...)

# S3 method for class 'drmrmodels'
mcmc_diag(x, ...)
```

## Arguments

- x:

  An object of class `adrm` or `sdm`.

- ...:

  currently ignored.

## Value

an object of class `"drmrdiag"`

## Author

lcgodoy
