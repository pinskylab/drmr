# Simulate AR(1)

Simulate AR(1)

## Usage

``` r
sim_ar(pars, n_time)
```

## Arguments

- pars:

  a named `list` with two elements: `alpha` (representing the temporal
  autocorrelation parameter) and `tau` representing the conditional
  standard deviation of the AR(1) process.

- n_time:

  number of time points

## Value

a `vector` of length `n_time` representing a realization of a zero-mean
AR(1) process.

## Author

lcgodoy
