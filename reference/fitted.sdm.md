# SDM fitted values

Calculates the log-likelihood associated to a `sdm` obkect.

## Usage

``` r
# S3 method for class 'sdm'
fitted(object, type = "predictive", cores = 1, ...)
```

## Arguments

- object:

  A `sdm` object containing the output from the
  [fit_sdm](https://pinskylab.github.io/drmr/reference/fit_sdm.md)
  function.

- type:

  type of predictions to be computed. Admitted values are

  - `"predictive"` (default): posterior predictive distribution;

  - `"expected"`: theoretical mean of the posterior predictive
    distribution;

  - `"latent"`: latent density (i.e., disconsidering the observation
    error);

- cores:

  number of threads used for the forecast. If four chains were used in
  the `sdm`, then four (or less) threads are recommended.

- ...:

  additional parameters to be passed to `$generated_quantities`

## Value

an object of class `"CmdStanGQ"` containing samples for the posterior
predictive distribution for the observed data.

## Author

lcgodoy
