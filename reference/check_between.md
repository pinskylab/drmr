# Check if x is between lb and ub

Check if x is between lb and ub

## Usage

``` r
check_between(x, lb, ub)
```

## Arguments

- x:

  A numeric vector.

- lb:

  A numeric vector of lower bounds.

- ub:

  A numeric vector of upper bounds.

## Value

Returns `NULL` invisibly. This function stops execution if any of the
following conditions are met:

- lb is greater than or equal to ub.

- The lengths of lb and ub are not equal.

- The lengths of lb, ub, and x are not equal.
