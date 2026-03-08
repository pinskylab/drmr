# Age-specific expected densities based on DRM.

Age-specific expected densities based on DRM.

## Usage

``` r
ages_edens(drm, cores = 1)

lambda_drm(...)
```

## Arguments

- drm:

  A `list` object containing the output from the
  [fit_drm](https://pinskylab.github.io/drmr/reference/fit_drm.md)
  function.

- cores:

  number of threads used for the forecast. If four chains were used in
  the `drm`, then four (or less) threads are recommended.

- ...:

  params to be passed to `age_edens`

## Value

an object of class `"CmdStanGQ"` containing samples for the posterior
predictive distribution for forecasting.

## Author

lcgodoy
