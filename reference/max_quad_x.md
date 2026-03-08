# Value of a covariate that maximizes the response variable in a quadratic model.

Consider a linear predictor having linear and square terms associated
with a variable \\x\\. Assume this variable was centered before being
included in the linear predictor. This function returns the value of
\\x\\ (on its original scale) such that the linear predictor is
maximized (or minimized).

## Usage

``` r
max_quad_x(beta1, beta2, offset = 0)
```

## Arguments

- beta1:

  A `numeric` regression coefficient associated with the linear term.

- beta2:

  A `numeric` regression coefficient associated with the quadratic term.

- offset:

  a `numeric` representing the "center" of \\x\\.

## Value

A `numeric` value representing the uncentered \\x\\ that maximizes (or
minimizes) the linear predictor.

## Author

Lucas Godoy
