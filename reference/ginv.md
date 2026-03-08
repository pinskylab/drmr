# Generalized Inverse

Generalized Inverse

## Usage

``` r
ginv(X, tol = sqrt(.Machine$double.eps))
```

## Arguments

- X:

  A matrix we wish to invert.

- tol:

  A relative tolerance to detect zero singular values.

## Value

The generalized inverse of `X`.

## Details

this function is taken from the package `MASS`
