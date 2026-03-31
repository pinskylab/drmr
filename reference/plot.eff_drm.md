# Plot Effects_Drminal Relationships for ADRM Objects

Automatically plots the effects_drminal relationship computed by
[`effects_drm()`](https://pinskylab.github.io/drmr/reference/effects_drm.md).
If `ggplot2` is installed, it returns a ggplot object; otherwise, it
falls back to base R graphics.

## Usage

``` r
# S3 method for class 'eff_drm'
plot(x, rug_data = NULL, ...)
```

## Arguments

- x:

  An object of class `eff_drm`, usually the output of
  [`effects_drm()`](https://pinskylab.github.io/drmr/reference/effects_drm.md).

- rug_data:

  An optional `data.frame` containing the original data to add a rug
  plot to the x-axis.

- ...:

  Additional arguments passed to the underlying plotting functions.

## Details

This default plotting method is restricted to outputs where exactly one
focal variable was evaluated (`length(variable) == 1`). Visualizing more
complicated cases (e.g., 2D interactions) should be handled manually by
the user.
