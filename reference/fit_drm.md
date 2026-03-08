# Fit the dynamic range model.

Fit the DRM Stan model

## Usage

``` r
fit_drm(
  .data,
  y_col,
  time_col,
  site_col,
  family = "gamma",
  formula_zero = ~1,
  formula_rec = ~1,
  formula_surv = NULL,
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

- formula_rec:

  A `formula` specifying the model for the recruitment component.
  Defaults to `~ 1` (intercept only).

- formula_surv:

  A `formula` specifying the model for the survival component. If `NULL`
  (the default), no survival component is included.

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

An object of class `adrm` which is a `list` containing the MCMC draws,
the model data, the linear predictors formulas, and the (response, time,
site) column names. Specifically:

- `stanfit`: The MCMC draws from the fitted model.

- `data`: The data used to fit the model (as a list).

- `formulas`: The formulas used to create design matrices.

- `cols`: Important column names.

## See also

[`make_data()`](https://pinskylab.github.io/drmr/reference/make_data.md)

Other models:
[`fit_sdm()`](https://pinskylab.github.io/drmr/reference/fit_sdm.md)

## Author

lcgodoy

## Examples

``` r
if (instantiate::stan_cmdstan_exists()) {
  data(sum_fl)
  fit_drm(.data = sum_fl,
          y_col = "y",
          time_col = "year",
          site_col = "patch",
          seed = 2025)$stanfit$summary()
}
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 1 '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 312, column 2 to line 315, column 67)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration:  100 / 2000 [  5%]  (Warmup) 
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
#> Chain 1 finished in 3.3 seconds.
#> Chain 2 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Iteration:  100 / 2000 [  5%]  (Warmup) 
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
#> Chain 2 finished in 4.3 seconds.
#> Chain 3 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 3 '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 312, column 2 to line 315, column 67)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 3 '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 312, column 2 to line 315, column 67)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 3 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 3 '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 312, column 2 to line 315, column 67)
#> Chain 3 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 3 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 3 
#> Chain 3 Iteration:  100 / 2000 [  5%]  (Warmup) 
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
#> Chain 3 finished in 4.0 seconds.
#> Chain 4 Iteration:    1 / 2000 [  0%]  (Warmup) 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: Exception: gamma_lpdf: Inverse scale parameter[1] is 0, but must be positive finite! (in '/tmp/Rtmp77Yui3/pkg-lib1aba1e483534/drmr/bin/stan/utils/lpdfs.stan', line 97, column 4, included from
#> Chain 4 '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 2, column 0) (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 312, column 2 to line 315, column 67)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 4 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpJsUOE8/model-27b4495e770c.stan', line 309, column 4 to column 54)
#> Chain 4 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 4 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 4 
#> Chain 4 Iteration:  100 / 2000 [  5%]  (Warmup) 
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
#> Chain 4 finished in 2.7 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 3.6 seconds.
#> Total execution time: 14.6 seconds.
#> 
#> Warning: 137 of 4000 (3.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> # A tibble: 1,076 × 10
#>    variable        mean    median     sd    mad       q5      q95  rhat ess_bulk
#>    <chr>          <dbl>     <dbl>  <dbl>  <dbl>    <dbl>    <dbl> <dbl>    <dbl>
#>  1 lp__       -1363.     -1.36e+3 1.46   1.27   -1.37e+3 -1.36e+3  1.01     602.
#>  2 lxi[1]        -2.64   -2.31e+0 1.33   1.15   -5.31e+0 -1.14e+0  1.01     209.
#>  3 phi[1]         0.721   7.18e-1 0.0554 0.0546  6.32e-1  8.15e-1  1.01    1155.
#>  4 beta_r[1]      3.38    3.38e+0 0.0770 0.0777  3.26e+0  3.51e+0  1.01     753.
#>  5 beta_t[1]     -0.180  -2.63e-1 0.400  0.415  -6.96e-1  5.89e-1  1.02     194.
#>  6 xi[1]         -0.124  -9.94e-2 0.102  0.109  -3.20e-1 -4.93e-3  1.01     209.
#>  7 log_rec[1]     3.38    3.38e+0 0.0770 0.0777  3.26e+0  3.51e+0  1.01     753.
#>  8 log_rec[2]     3.38    3.38e+0 0.0770 0.0777  3.26e+0  3.51e+0  1.01     753.
#>  9 log_rec[3]     3.38    3.38e+0 0.0770 0.0777  3.26e+0  3.51e+0  1.01     753.
#> 10 log_rec[4]     3.38    3.38e+0 0.0770 0.0777  3.26e+0  3.51e+0  1.01     753.
#> # ℹ 1,066 more rows
#> # ℹ 1 more variable: ess_tail <dbl>
```
