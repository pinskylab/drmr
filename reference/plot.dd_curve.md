# Plot Stock-Recruitment Curves for ADRM Objects

Automatically plots the stock-recruitment curves computed by
[`dd_curves()`](https://pinskylab.github.io/drmr/reference/dd_curves.md).
Uses ggplot2 if available, otherwise base R graphics.

## Usage

``` r
# S3 method for class 'dd_curve'
plot(x, ...)
```

## Arguments

- x:

  An object of class `dd_curve`, usually the output of
  [`dd_curves()`](https://pinskylab.github.io/drmr/reference/dd_curves.md).

- ...:

  Additional arguments passed to the underlying plotting functions.
