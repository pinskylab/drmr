# (out-of-sample) Expected Log-posterior Density (ELPD) based on `adrm` and `sdm` objects.

Considering a new dataset (across the same sites), computes the
out-of-sample ELPD based on the DRM passed as `drm`.

## Usage

``` r
elpd(x, ...)

# S3 method for class 'adrm'
elpd(x, new_data, past_data, f_test, seed = 1, cores = 1, ...)

# S3 method for class 'sdm'
elpd(x, new_data, seed = 1, cores = 1, ...)
```

## Arguments

- x:

  A `list` object containing the output from the
  [`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md)
  (or
  [`fit_sdm()`](https://pinskylab.github.io/drmr/reference/fit_sdm.md))
  function.

- ...:

  additional parameters to be passed to `elpd`

- new_data:

  a `data.frame` with the dataset at which we wish to obtain
  predictions. Note that, this `data.frame` must contain the response
  variable used when fitting the DRM as well.

- past_data:

  a `data.frame` with the dataset last year used in model fitting. Only
  needed when `f_test` is not missing or when estimating survival.

- f_test:

  a `matrix` informing the instantaneous fishing mortality rates at each
  age (columns) and timepoint (rows).

- seed:

  a seed used for the forecasts. Forecasts are obtained through Monte
  Carlo samples from the posterior predictive distribution. Therefore, a
  `seed` is needed to ensure the results' reproducibility.

- cores:

  number of threads used for the forecast. If four chains were used in
  the `drm`, then four (or less) threads are recommended.

## Value

an object of class `"CmdStanGQ"` containing the ELPD function evaluated
at each data point given each sample from the posterior.

## Details

The current version of the code assumes the data where forecasts are
needed are ordered by "site" and "site" and, in addition, the sites MUST
be the same as the ones used to obtain the parameters' estimates from
the the `drm` object.

## Author

lcgodoy
