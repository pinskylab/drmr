# Estimate phi

Estimate phi

## Usage

``` r
get_phi_hat(y, family)
```

## Arguments

- y:

  a `numeric vector` of species' densities.

- family:

  a `character` specifying the family of the probability distribution
  assumed for density. The options are:

  - `"gamma"` (default): gamma parametrized in terms of its mean;

  - `"lognormal"`: log-normal parametrized in terms of its mean;

  - `"loglogistic"`: log-logistic parametrized in terms of its median
    (usual parametrization);

  - `"lognormal_legacy"`: log-normal with its usual parametrization;

## Value

A `numeric` scalar representing an estimate for phi.

## Author

lcgodoy
