# Plot Diagnostics for DRM Models

Provides trace and density plots for the posterior draws of a DRM (or
SDM) model.

## Usage

``` r
# S3 method for class 'drmrmodels'
plot(x, variables = NULL, type = c("trace", "density"), ask = NULL, ...)
```

## Arguments

- x:

  An object of class `drmrmodels`. Typically, the output of a `fit_drm`
  or `fit_sdm` call.

- variables:

  An optional character vector of parameter names to plot. If `NULL`,
  plots all.

- type:

  A character vector specifying which plots to draw: `"trace"`,
  `"density"`, or both.

- ask:

  Logical; if `TRUE`, the user is prompted before displaying the next
  plot. If `NULL` (the default), it is calculated dynamically based on
  the current graphical layout.

- ...:

  Additional graphical parameters passed to `plot`.
