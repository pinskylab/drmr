# Simulate response variable

Given the expected value `x` and a list of parameters (`pars`), this
function simulates one realization of the probability density function
which is also specified by the list `pars`.

## Usage

``` r
sim_dens(x, pars)
```

## Arguments

- x:

  an expected value for a simulated value.

- pars:

  a `list` of parameters.

## Value

a `numeric` value representing a realization of the probability
distribution with mean `x` and parameters `pars`.

## Author

lcgodoy
