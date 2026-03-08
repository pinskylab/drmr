# Generate the "survival" terms

Generate the "survival" terms

## Usage

``` r
make_surv(n_patches, n_time, x_sv, pars)
```

## Arguments

- n_patches:

  number of patches

- n_time:

  number of timepoints

- x_sv:

  matrix of environmental factors affecting survival

- pars:

  a named `list` of parameters used to simulate log-recruitment. It must
  contain a vector named `"beta_s"` with length equal to the number of
  columnts in `x_sv`.

## Value

a `matrix` with `n_patches` columns and `n_time` rows representing the
log-survival at each patch/site and time.

## Author

lcgodoy
