# Marginal Relationships with Covariates

Evaluates and summarizes the marginal relationships between explanatory
variables and recruitment, survival, or absence probability from a
fitted DRM model.

## Usage

``` r
marg(object, ...)

# S3 method for class 'adrm'
marg(
  object,
  process = c("rec", "surv", "pabs"),
  variable,
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

  Additional arguments passed to methods.

- process:

  A character string indicating the process to evaluate: `"rec"`
  (recruitment), `"surv"` (survival), or `"pabs"` (probability of
  absence).

- variable:

  A character vector with the name(s) of the focal variable(s) to
  examine.

- newdata:

  An optional `data.frame` containing the values for the focal
  variable(s). If `NULL`, a grid is generated automatically based on the
  observed range in the model matrix.

- n_pts:

  An integer specifying the number of points to generate for the
  sequence of each focal variable when `newdata` is `NULL`. Default is
  100.

- summary:

  Logical. If `TRUE` (the default), returns the quantiles of the
  posterior predictions. If `FALSE`, returns the raw posterior draws.

- prob:

  A numeric scalar in \\(0, 1)\\ specifying the probability mass of the
  equal-tailed credible interval. Defaults to `0.9`, which produces a
  90\\ together with the median.

## Value

A `data.frame` with the posterior summaries (or draws) for the specified
process. If `summary = TRUE`, it also receives the class `marg_adrm` to
enable automated plotting.

## Details

The `marg` function computes the predicted relationships across a
sequence of values for a focal variable (or variables), holding all
other non-focal variables in the model matrix at zero.

When `summary = TRUE`, the function calculates an equal-tailed credible
interval and the median using
[`quantile2`](https://mc-stan.org/posterior/reference/quantile2.html),
which is highly optimized for posterior draws.

## Author

lcgodoy
