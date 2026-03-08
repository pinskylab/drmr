# Fit a GLM based SDM

Fit the SDM Stan model

## Usage

``` r
fit_sdm(
  .data,
  y_col,
  time_col,
  site_col,
  family = "gamma",
  formula_zero = ~1,
  formula_dens = ~1,
  seed,
  init = "cmdstan_default",
  algorithm = "nuts",
  algo_args = list(),
  ...
)
```

## Arguments

- .data:

  A `data frame` containing the data for the model.

- y_col:

  A `character` specifying the name of the column in `.data` that
  contains the response variable.

- time_col:

  A `character` specifying the name of the column in `.data` that
  contains the time variable.

- site_col:

  A `character` specifying the name of the column in `.data` that
  contains the site variable.

- family:

  a `character` specifying the family of the probability distribution
  assumed for density. The options are:

  - `"gamma"` (default): gamma parametrized in terms of its mean;

  - `"lognormal"`: log-normal parametrized in terms of its mean;

  - `"loglogistic"`: log-logistic parametrized in terms of its median
    (usual parametrization);

  - `"lognormal_legacy"`: log-normal with its usual parametrization;

- formula_zero:

  A `formula` specifying the model for the zero inflation component.
  Defaults to `~ 1` (intercept only).

- formula_dens:

  A `formula` specifying the model for the non-zero density component.
  Defaults to `~ 1` (intercept only).

- seed:

  An `integer` specifying the random number seed.

- init:

  A scalar specifying the initialization method. The default
  ("cmdstan_default") lets `cmdstan` initialize parameters. Other
  options include: a scalar greater than zero, say `x`, which
  initializes all parameters uniformly between `-x` and `x`; `0`, which
  initializes all parameters at `0`; or "prior", which initializes
  parameters by sampling from their priors.

- algorithm:

  a `character` specifying the algorithm used for inference. Default is
  `nuts` (the default MCMC in Stan). The remaining options are different
  flavors of variational bayes algorithms: "vb" (for ADVI), "pathfinder"
  (for Pathfinder), "laplace" (normal approximation centered at the mode
  of the posterior) or "optimize" for (penalized) MLEs.

- algo_args:

  a `list` with arguments for the sampling algorithms. For instance,
  `tol_rel_obj` for variational inference.

- ...:

  Passed on to the
  [`make_data()`](https://pinskylab.github.io/drmr/reference/make_data.md)
  function used to build the input `list` for our `cmdstanr` model.

## Value

An object of class `sdm` which is a `list` containing the MCMC draws,
the model data, the linear predictors formulas, and the (response, time,
site) column names.

- `stanfit`: The MCMC draws from the fitted model.

- `data`: The data used to fit the model (as a list).

- `formulas`: The formulas used to create design matrices.

- `cols`: Important column names.

## See also

[`make_data_sdm()`](https://pinskylab.github.io/drmr/reference/make_data_sdm.md)

Other models:
[`fit_drm()`](https://pinskylab.github.io/drmr/reference/fit_drm.md)

## Author

lcgodoy

## Examples

``` r
if (instantiate::stan_cmdstan_exists()) {
  data(sum_fl)
  fit_sdm(.data = sum_fl,
          y_col = "y",
          time_col = "year",
          site_col = "patch",
          seed = 2025)$stanfit$summary()
}
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 1 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 1 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 1 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 1 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 1 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 1 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 1 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 1 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 1 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 1 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 1 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 1 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 1 finished in 2.1 seconds.
#> Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 2 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 2 '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 193, column 2 to line 196, column 67)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 2 '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 193, column 2 to line 196, column 67)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 2 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 2 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 2 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 2 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 2 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 2 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 2 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 2 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 2 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 2 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 2 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 2 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 2 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 2 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 2 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 2 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 2 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 2 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 2 finished in 1.9 seconds.
#> Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 3 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 3 '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 193, column 2 to line 196, column 67)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 3 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 3 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 3 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 3 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 3 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 3 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 3 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 3 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 3 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 3 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 3 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 3 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 3 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 3 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 3 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 3 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 3 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 3 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 3 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 3 finished in 2.1 seconds.
#> Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 4 Iteration:  100 / 2000 [  5%]  (Warmup) 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 190, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 4 '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b430afe1a1.stan', line 193, column 2 to line 196, column 67)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Iteration:  200 / 2000 [ 10%]  (Warmup) 
#> Chain 4 Iteration:  300 / 2000 [ 15%]  (Warmup) 
#> Chain 4 Iteration:  400 / 2000 [ 20%]  (Warmup) 
#> Chain 4 Iteration:  500 / 2000 [ 25%]  (Warmup) 
#> Chain 4 Iteration:  600 / 2000 [ 30%]  (Warmup) 
#> Chain 4 Iteration:  700 / 2000 [ 35%]  (Warmup) 
#> Chain 4 Iteration:  800 / 2000 [ 40%]  (Warmup) 
#> Chain 4 Iteration:  900 / 2000 [ 45%]  (Warmup) 
#> Chain 4 Iteration: 1000 / 2000 [ 50%]  (Warmup) 
#> Chain 4 Iteration: 1001 / 2000 [ 50%]  (Sampling) 
#> Chain 4 Iteration: 1100 / 2000 [ 55%]  (Sampling) 
#> Chain 4 Iteration: 1200 / 2000 [ 60%]  (Sampling) 
#> Chain 4 Iteration: 1300 / 2000 [ 65%]  (Sampling) 
#> Chain 4 Iteration: 1400 / 2000 [ 70%]  (Sampling) 
#> Chain 4 Iteration: 1500 / 2000 [ 75%]  (Sampling) 
#> Chain 4 Iteration: 1600 / 2000 [ 80%]  (Sampling) 
#> Chain 4 Iteration: 1700 / 2000 [ 85%]  (Sampling) 
#> Chain 4 Iteration: 1800 / 2000 [ 90%]  (Sampling) 
#> Chain 4 Iteration: 1900 / 2000 [ 95%]  (Sampling) 
#> Chain 4 Iteration: 2000 / 2000 [100%]  (Sampling) 
#> Chain 4 finished in 1.9 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 2.0 seconds.
#> Total execution time: 8.2 seconds.
#> 
#> Warning: 39 of 4000 (1.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> # A tibble: 706 × 10
#>    variable       mean    median     sd    mad        q5      q95  rhat ess_bulk
#>    <chr>         <dbl>     <dbl>  <dbl>  <dbl>     <dbl>    <dbl> <dbl>    <dbl>
#>  1 lp__      -1366.    -1366.    1.47   1.26   -1369.    -1.36e+3  1.01     711.
#>  2 lxi[1]       -2.59     -2.25  1.37   1.19      -5.31  -9.37e-1  1.03     121.
#>  3 phi[1]        0.721     0.718 0.0557 0.0534     0.634  8.17e-1  1.00    1313.
#>  4 beta_r[1]     3.90      3.90  0.0772 0.0764     3.77   4.02e+0  1.00    1543.
#>  5 beta_t[1]    -0.128    -0.243 0.483  0.450     -0.694  8.21e-1  1.03     109.
#>  6 xi[1]        -0.138    -0.106 0.123  0.116     -0.392 -4.96e-3  1.03     121.
#>  7 rho[1]        0.340     0.339 0.0244 0.0240     0.300  3.80e-1  1.00    3993.
#>  8 rho[2]        0.340     0.339 0.0244 0.0240     0.300  3.80e-1  1.00    3993.
#>  9 rho[3]        0.340     0.339 0.0244 0.0240     0.300  3.80e-1  1.00    3993.
#> 10 rho[4]        0.340     0.339 0.0244 0.0240     0.300  3.80e-1  1.00    3993.
#> # ℹ 696 more rows
#> # ℹ 1 more variable: ess_tail <dbl>
```
