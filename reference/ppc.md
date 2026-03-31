# Posterior Predictive Checks for DRM Models

Generates posterior predictive check plots comparing observed data to
simulated data from the posterior predictive distribution.

## Usage

``` r
ppc(x, ...)

# S3 method for class 'drmrmodels'
ppc(
  x,
  type = c("density", "ecdf"),
  npost = 50,
  transform = c("none", "log1p"),
  ...
)
```

## Arguments

- x:

  An object of class `drmrmodels`.

- ...:

  Additional graphical parameters passed to `plot`.

- type:

  Character string indicating the type of plot: `"density"` (default) or
  `"ecdf"`.

- npost:

  Integer specifying the number of posterior draws to plot. Default
  is 50. Usually, we do not use the total number of draws here because
  the plot tends to get too heavy.

- transform:

  Character string specifying the scale for plotting. Either `"none"`
  (default) or `"log1p"`.
