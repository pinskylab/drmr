# Simulate log-recruitment

Simulate log-recruitment

## Usage

``` r
sim_log_rec(n_patches, n_time, x_rec, pars, ar_time = TRUE)
```

## Arguments

- n_patches:

  number of patches

- n_time:

  number of timepoints

- x_rec:

  matrix of environmental factors affecting recruitment

- pars:

  a named `list` of parameters used to simulate log-recruitment. It must
  contain a vector named `"beta_r"` with length equal to the number of
  columnts in `x_rec`. In addition, if `ar_time = TRUE`, the list must
  also contain a named `vector` called "ar". This named `vector` must
  contain an element called `alpha` (the autocorrelation parameter) and
  another called `tau` (the conditional SD).

- ar_time:

  a `boolean` indicating whether an AR(1) term should be included to the
  log-recruitment.

## Value

a `matrix` with `n_patches` columns and `n_time` rows representing the
log-recruitment at each patch/site and time.

## Author

lcgodoy
