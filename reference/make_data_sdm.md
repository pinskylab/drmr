# Make data for SDM stan models

This function creates the `list` used as the input for the `stan` model.

## Usage

``` r
make_data_sdm(
  y,
  time,
  site,
  z,
  x,
  adj_mat = matrix(0, ncol = 1, nrow = 1),
  .toggles,
  .priors,
  family = "gamma",
  reorder = TRUE,
  phi_hat = FALSE
)
```

## Arguments

- y:

  a `numeric vector` of species' densities.

- time:

  a `vector` indicating the time point associated to each element of
  `y`.

- site:

  a `vector` indicating the sites associated to each element of `y`.

- z:

  a design `matrix` of variables associated to the probability of
  absence at each site/time.

- x:

  a design `matrix` of variables associated to the non-zero densities.

- adj_mat:

  an adjacency `matrix` of dimensions `sites` \\\times\\ `sites`. Its
  elements are 1 if two sites are neighbors and zero otherwise.

- .toggles:

  a `list` of toggles for model components. The components are:

  - `cloglog`: 1 to use the complementary log-log and 0 for the logit
    link function for the absence probabilities.

  - `ar_re`: 1 to incorporate an AR(1) process to density and 0
    otherwise.

  - `iid_re`: 1 to incorporate a site specific IID random effect to
    density and 0 otherwise.

  - `sp_re`: 1 to incorporate a site specific ICAR random effect to
    density and 0 otherwise.

- .priors:

  a `list` of priors hyperparameters.

- family:

  a `character` specifying the family of the probability distribution
  assumed for density. The options are:

  - `"gamma"` (default): gamma parametrized in terms of its mean;

  - `"lognormal"`: log-normal parametrized in terms of its mean;

  - `"loglogistic"`: log-logistic parametrized in terms of its median
    (usual parametrization);

  - `"lognormal_legacy"`: log-normal with its usual parametrization;

- reorder:

  a `boolean` telling whether the data needs to be reordered. The
  default is TRUE and means the data points will be ordered by site and
  time, respectively.

- phi_hat:

  a `boolean` indicating whether the prior on `phi` should be determined
  through the data (using
  [`get_phi_hat()`](https://pinskylab.github.io/drmr/reference/get_phi_hat.md)).

## Value

A `list` containing the data and settings to be used as input for the
Stan model.

## Author

lcgodoy
