# Predictions based on DRM.

Considering a new dataset (across the same sites), computes predictions
based on the DRM passed as `drm`.

## Usage

``` r
# S3 method for class 'adrm'
predict(
  object,
  new_data,
  past_data,
  f_test,
  type = "predictive",
  seed = 1,
  cores = 1,
  ...
)

predict_drm(...)
```

## Arguments

- object:

  An `adrm` object containing the output from the
  [`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md)
  function.

- new_data:

  a `data.frame` with the dataset at which we wish to obtain
  predictions.

- past_data:

  a `data.frame` with the dataset from the last year used in model
  fitting. Only needed when `f_test` is not missing or when estimating
  survival.

- f_test:

  a `matrix` informing the instantaneous fishing mortality rates at each
  age (columns) and timepoint (rows).

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
  used in the `drm`, then four (or less) threads are recommended.

- ...:

  params to be passed to `predict.adrm`

## Value

an object of class `"CmdStanGQ"` containing samples for the posterior
predictive distribution for predictions.

## Details

The current version of the code assumes the data where predictions are
needed is ordered by "site" and "site" and, in addition, its sites MUST
be the same as the ones used to obtain the parameters' estimates from
the the `drm` object.

## See also

[`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md)

## Author

lcgodoy
