# Plot Marginal Relationships for ADRM Objects

Automatically plots the marginal relationship computed by
[`marg()`](https://pinskylab.github.io/drmr/reference/marg.md). If
`ggplot2` is installed, it returns a ggplot object; otherwise, it falls
back to base R graphics.

## Usage

``` r
# S3 method for class 'marg_adrm'
plot(x, rug_data = NULL, ...)
```

## Arguments

- x:

  An object of class `marg_adrm`, usually the output of
  [`marg()`](https://pinskylab.github.io/drmr/reference/marg.md).

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
