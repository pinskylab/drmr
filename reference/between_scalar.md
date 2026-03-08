# Check if elements in x are between corresponding elements in lb and ub

Check if elements in x are between corresponding elements in lb and ub

## Usage

``` r
between_scalar(x, lb, ub)
```

## Arguments

- x:

  A numeric vector.

- lb:

  A numeric vector of lower bounds.

- ub:

  A numeric vector of upper bounds.

## Value

A logical vector of the same length as x, indicating whether each
element of x is between the corresponding elements of lb and ub.

## Examples

``` r
between(1:5, 1, 5)
#> [1] TRUE TRUE TRUE TRUE TRUE
between(1:5, 2, 4)
#> [1] FALSE  TRUE  TRUE  TRUE FALSE
```
