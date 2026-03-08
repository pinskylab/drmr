# Computing the log-likelihood function fot `adrm` and `sdm` objets.

Computing the log-likelihood function fot `adrm` and `sdm` objets.

## Usage

``` r
log_lik(x, cores = 1)

# S3 method for class 'adrm'
log_lik(x, cores = 1)

# S3 method for class 'sdm'
log_lik(x, cores = 1)
```

## Arguments

- x:

  an object of class `adrm` or `sdm`. These objects are generated as the
  output of the
  [`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md)
  and
  [`fit_sdm()`](https://pinskylab.github.io/drmr/reference/fit_sdm.md)
  functions, respectively,

- cores:

  number of threads used to calculate the log-likelihood function. If
  four chains were used in the `fit_drm` (of `fit_sdm`), then four (or
  less) threads are recommended.

## Value

an object of class `"CmdStanGQ"` containing the log-likelihood function
evaluated at each data point given each sample from the posterior.

## Author

lcgodoy
