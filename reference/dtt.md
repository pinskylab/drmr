# Density of a truncated Student's t distribution

Calculates the density of a truncated Student's t distribution.

## Usage

``` r
dtt(x, mean = 0, sd = 1, df = 3, range = c(0, Inf), log = FALSE)
```

## Arguments

- x:

  A `numeric vector` at which to evaluate the density.

- mean:

  A `numeric` scalar representing the mean of the underlying t
  distribution. Defaults to 0.

- sd:

  A `numeric` scalar representing the standard deviation of the
  underlying t distribution. Defaults to 1.

- df:

  A `numeric` scalar representing the degrees of freedom of the
  underlying t distribution. Defaults to 3.

- range:

  A `numeric` vector of length 2 specifying the lower and upper
  truncation bounds. Defaults to `c(0, Inf)`, indicating truncation from
  below at 0.

- log:

  `Logical`; if `TRUE`, the natural logarithm of the density is
  returned. Defaults to `FALSE`.

## Value

A numeric vector of the same length as `x`, containing the (log) density
values of the truncated t distribution.

## See also

[`stats::dt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::pt()`](https://rdrr.io/r/stats/TDist.html)

## Author

lcgodoy
