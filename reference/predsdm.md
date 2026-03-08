# Predictions based on SDM.

Considering a new dataset (across the same sites), computes predictions
based on the SDM passed as `sdm`.

## Usage

``` r
# S3 method for class 'sdm'
predict(object, new_data, type = "predictive", seed = 1, cores = 1, ...)

predict_sdm(...)
```

## Arguments

- object:

  An `sdm` object containing the output of a
  [fit_sdm](https://pinskylab.github.io/drmr/reference/fit_sdm.md) call.

- new_data:

  a `data.frame` with the dataset at which we wish to obtain
  predictions.

- type:

  type of predictions to be computed. Admitted values are

  - `"predictive"` (default): posterior predictive distribution;

  - `"expected"`: theoretical mean of the posterior predictive
    distribution;

  - `"latent"`: latent density (i.e., disregarding the observation
    error);

- seed:

  a seed used for the predictions. predictions are obtained through
  Monte Carlo samples from the posterior predictive distribution.
  Therefore, a `seed` is needed to ensure the results' reproducibility.

- cores:

  number of threads used to compute the predictions. If four chains were
  used in the `sdm`, then four (or less) threads are recommended.

- ...:

  params to be passed to `predict.sdm`

## Value

An object of class `CmdStanGQ` (from the `instantiate` package)
containing samples for the posterior predictive distribution for
predictions.

## Details

The current version of the code assumes the data where predictions are
needed is ordered by "site" and "site" and, in addition, its sites MUST
be the same as the ones used to obtain the parameters' estimates from
the the `sdm` object.

## See also

[`fit_sdm()`](https://pinskylab.github.io/drmr/reference/fit_sdm.md)

## Author

lcgodoy
