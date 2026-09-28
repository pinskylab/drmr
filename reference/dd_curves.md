# Stock-Recruitment Curves for DD Models

Computes expected recruitment across a grid of Spawning Stock Biomass
(density) values, optionally conditioned on environmental regimes.

## Usage

``` r
dd_curves(object, ...)

# S3 method for class 'adrm'
dd_curves(
  object,
  density_grid = NULL,
  newdata = NULL,
  n_pts = 100,
  summary = TRUE,
  prob = 0.95,
  ...
)
```

## Arguments

- object:

  An object of class `adrm`, typically the output of
  [`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md).

- ...:

  Additional arguments.

- density_grid:

  A numeric vector of density values to evaluate. If `NULL`, a default
  sequence is generated from 0 to the maximum observed response.

- newdata:

  An optional `data.frame` of environmental covariates. Each row
  represents a distinct environmental regime.

- n_pts:

  Integer. Number of points for the default `density_grid` when it is
  not provided. Defaults to 100.

- summary:

  Logical. If `TRUE`, returns posterior quantiles. If `FALSE`, returns
  the raw posterior draws.

- prob:

  Numeric. Credible interval probability mass. Defaults to 0.95.

## Value

A `data.frame` of class `dd_curve` (if `summary = TRUE`) containing the
Stock-Recruitment evaluations.

## Author

lcgodoy
