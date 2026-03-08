# Density of a truncated Normal distribution

Calculates the density of a truncated Normal distribution.

## Usage

``` r
dtn(x, mean = 0, sd = 1, range = c(0, Inf), log = FALSE)
```

## Arguments

- x:

  A `numeric vector` at which to evaluate the density.

- mean:

  A `numeric` scalar representing the mean of the underlying normal
  distribution. Defaults to 0.

- sd:

  A `numeric` scalar representing the standard deviation of the
  underlying normal distribution. Defaults to 1.

- range:

  A `numeric` vector of length 2 specifying the lower and upper
  truncation bounds. Defaults to `c(0, Inf)`, indicating truncation from
  below at 0.

- log:

  `Logical`; if `TRUE`, the natural logarithm of the density is
  returned. Defaults to `FALSE`.

## Value

A numeric vector of the same length as `x`, containing the (log) density
values of the truncated normal distribution.

## See also

[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html),
[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html)

## Author

lcgodoy
