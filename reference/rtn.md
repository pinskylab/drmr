# Random number generation from a truncated Normal distribution

Generates random numbers from a truncated Normal distribution.

## Usage

``` r
rtn(n, mean = 0, sd = 1, range = c(0, Inf))
```

## Arguments

- n:

  An `integer` specifying the number of random samples to generate.

- mean:

  A `numeric` scalar representing the mean of the underlying normal
  distribution. Defaults to 0.

- sd:

  A `numeric` scalar representing the standard deviation of the
  underlying normal distribution. Defaults to 1.

- range:

  A `numeric vector` of length 2 specifying the lower and upper
  truncation bounds. Defaults to `c(0, Inf)`, indicating truncation from
  below at 0.

## Value

A `numeric vector` of length `n` containing random numbers drawn from
the specified truncated normal distribution.

## See also

[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html),
[`stats::qnorm()`](https://rdrr.io/r/stats/Normal.html)

## Author

lcgodoy
