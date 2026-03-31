# Random number generation from a truncated Student's t distribution

Generates random numbers from a truncated Student's t distribution.

## Usage

``` r
rtt(n, mean = 0, sd = 1, df = 3, range = c(0, Inf))
```

## Arguments

- n:

  An `integer` specifying the number of random numbers to generate.

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

  A `numeric vector` of length 2 specifying the lower and upper
  truncation bounds. Defaults to `c(0, Inf)`, indicating truncation from
  below at 0.

## Value

A `numeric vector` of length `n` containing random numbers drawn from
the specified truncated t distribution.

## Details

For details on the method used, see:
<https://stats.stackexchange.com/questions/567944/how-can-i-sample-from-a-shifted-and-scaled-student-t-distribution-with-a-specifi>

## See also

[`stats::rt()`](https://rdrr.io/r/stats/TDist.html),
[`stats::qt()`](https://rdrr.io/r/stats/TDist.html)

## Author

lcgodoy
