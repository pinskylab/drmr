# Make data for DRM stan models

This function creates the `list` used as the input for the `stan` model.

## Usage

``` r
make_data(
  y,
  time,
  site,
  init_data = numeric(0),
  f_mort,
  m = -log(0.7),
  x_t,
  x_m,
  x_r,
  n_ages = 2,
  age_selectivity,
  ages_movement,
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

- init_data:

  an optional vector (of length n_ages - 1) to initialize the population
  dynamics.

- f_mort:

  an optional `matrix` informing the instantaneous fishing mortality
  rates at each age (columns) and timepoint (rows).

- m:

  a `numeric` value corresponding to the instantaneous natural mortality
  rate. The default value for this is `-log(.7)`, as it implies a
  survival rate of 0.70 between age classes.

- x_t:

  a design `matrix` of variables associated to the probability of
  absence at each site/time.

- x_m:

  a design `matrix` of variables associated to survival.

- x_r:

  a design `matrix` of variables associated to recruitment.

- n_ages:

  an `integer` indicating the number of ages for the underlying
  population dynamic model.

- age_selectivity:

  a `numeric vector` with `n_ages` elements, where each element
  indicates the selectivity of a respective age. All the elements of
  this vector must lie between 0 and 1.

- ages_movement:

  An `integer` or a `numeric vector` specifying the ages at which
  individuals of the focal species are assumed to move. If
  `ages_movement` is an integer, individuals younger than this age are
  considered static (non-moving). If `ages_movement` is a numeric vector
  of length `n_ages`, it indicates movement capability for each age
  group. A value of `0` indicates the corresponding age group is static,
  while `1` indicates movement is allowed. For example,
  `c(0, 0, 1, 1, 0)` specifies that age groups 1, 2, and 5 are static,
  while 3 and 4 are mobile.

- adj_mat:

  an adjacency `matrix` of dimensions `sites` \\\times\\ `sites`. Its
  elements are 1 if two sites are neighbors and zero otherwise.

- .toggles:

  a `list` of toggles for model components. The components are:

  - `rho_mu`: 1 to explicitly relate rho to mu and 0 otherwise.

  - `cloglog`: 1 to use the complementary log-log and 0 for the logit
    link function for the absence probabilities.

  - `movement`: 1 to allow for (adjacent) movement; 0 for static.

  - `est_surv`: 1 to estimate mortality and 0 otherwise.

  - `est_init`: 1 to estimate initial values for lambda and 0 otherwise.

  - `minit`: 1 to use mortality to estimate initial age classes and 0
    otherwise.

  - `ar_re`: a `character`. It assumes one of the following values:
    "none" no AR, "rec" AR(1) for recruitment, "surv" AR(1) for survival
    (only works when `est_surv` is on), "dens" AR(1) for density.

  - `iid_re`: a `character`. It assumes one of the following values:
    "none" no IID random effect, "rec" IID random effect for
    recruitment, "surv" IID random effect for survival (only works when
    `est_surv` is on), "dens" IID random effect for density.

  - `sp_re`: a `character`. It assumes one of the following values:
    "none" no ICAR random effect, "rec" ICAR random effect for
    recruitment, "surv" ICAR random effect for survival (only works when
    `est_surv` is on), "dens" ICAR random effect for density.

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
