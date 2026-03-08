# Regression coefficient for non-centered variable

Consider a linear predictor having linear and square terms associated
with a variable \\x\\. Assume this variable was centered before being
included in the linear predictor. This function recovers the regression
coefficient associated with the linear term as if the variable was not
centered.

## Usage

``` r
fix_linbeta(beta1, beta2, offset)
```

## Arguments

- beta1:

  A `numeric` regression coefficient associated with the linear term.

- beta2:

  A `numeric` regression coefficient associated with the quadratic term.

- offset:

  a `numeric` representing the "center" of \\x\\.

## Value

A `numeric` value representing the regression coefficient of the linear
term for the model where \\x\\ is not centered.

## Author

Lucas Godoy
