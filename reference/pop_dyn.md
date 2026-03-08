# Simulate population dynamics

Given a set of parameters, simulate population dynamics.

## Usage

``` r
pop_dyn(
  n_patches,
  n_time,
  n_ages,
  x_rec,
  ar_time,
  f_a_t,
  x_sv = NULL,
  m = 0.25,
  pars,
  init,
  init_type,
  movement = FALSE,
  adj_mat = NULL,
  mov_age = NULL
)
```

## Arguments

- n_patches:

  number of patches

- n_time:

  number of timepoints

- n_ages:

  a `integer` indicating the number of (assumed) age-classes.

- x_rec:

  matrix of environmental factors affecting recruitment

- ar_time:

  a `boolean` indicating whether an AR(1) term should be included to the
  log-recruitment.

- f_a_t:

  fishing mortality

- x_sv:

  matrix of environmental factors affecting survival

- m:

  a `numeric` value indicating the natural mortality instantaneous rate.

- pars:

  a named `list` of model parameters.

- init:

  initialization vector of length `n_ages - 1`

- init_type:

  type of initialization (integer between 0 and 5)

- movement:

  a `boolean` indicating whether movement should be applied or not. If
  TRUE, than `pars` must have an element called `zeta` indicating the
  probability of staying at a given patch between two timepoints.

- adj_mat:

  a `n_patches` by `n_patches` row-standardized adjacency `matrix`.

- mov_age:

  a `vector` of ages at which movement starts.

## Value

an array of expected densities per age-group, patch, and timepoint.

## Author

lcgodoy
