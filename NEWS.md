# drmr 0.1.3

* Modified `examples.qmd` to prevent errors due to pathfinder failing.

* Updated github actions.

* Same problematic behaviour that was affecting `make_data` had to be dealed
  with in prior initialization.
  
* When `family = "gamma"` and a prior for `phi` is not provided, the default
  prior for this parameter is a Gamma with `shape = 2` and `rate = ybar *
  s2`. This prior has its theoretical mean at the method of moments estimator of
  `phi`.

* When `family = "lognormal"` and a prior for `phi` is not provided, the default
  prior for this parameter is a Gamma having its theoretical mean at the MLE of
  `phi`.

* fixed typo when defining log-logistic reparametrization

# drmr 0.1.2.9

* fixed weird behaviour with `make_data`

# drmr 0.1.2

* `lambda_drm` function added to recover the age-specific densities form a
  `fit_drm` call.
  
* fixed a bug on the `forecast.stan` script.

# drmr 0.1.1

* Vignettes updated so `R CMD check` "works".

* The `lambda`s are no longer returned by the fit and predict
  functions. Returning those values imply a significant decrease in the
  computational performance, especially as increasing the sample size.

* The AR term in the model was "corrected".

* Some bugs with prior initialization for the sdm were fixed.

# drmr 0.1.0

* Exporting `fix_linbeta`, `max_quad_x`, and `int_score` functions

* Fixing `check_between` function (not exported; This functions is just a
  helper)
   
* Changed prior on `zeta` and allowing for user to input the
  hyperparameters. Before, we had a standard normal prior on the logit of
  `zeta`. Now, we place a beta prior on `zeta`.

* Included support for different types of population dynamics initialization

* Fixed how pdf and random number generation from location-scale Student's t
  distribution

* `rlang` became a dependency.

* Included functions for simulating data from the DRM model and prior predictive
  checks.
   
* Trying to avoid overflow by making calculations on the log-scale whenever it's
  possible.

* Constraining `alpha` to $(0, 1)$.

* Prior on `phi` now is Gamma.

* `fit_drm` and `fit_sdm` returns changed. Now, the element `draws` is called
  `stanfit`. In addition, there is an additional element to the returned list
  called `formulas`. The `formulas` elements aims at making the `predict_*`
  functions less error prone.

* The QR parametrization toggles were completely removed.

* The `coef_*` parameters were converted to `beta_*`. That is, `coef_r` now is
  `beta_r`; while `coef_t` now is `beta_t`, and so on.

* The `pr_logsd_r_*` inputs converted to `pr_ltau_*`.

* Functions to plot effect of covariates on recruitment, survival, or absence
  probability were included.

* `est_mort` becomes `est_surv` (as it makes more sense with the text).

* `get-started` and `examples` vignette updated.

# drmr 0.0.24

* Initial values for population dynamics have been fixed.

* Vignettes to quarto

* New `between` function based on `data.table::between`.

* `int_score` function to calculate the interval score was also included. The
  interval score helps to assess interval predictions.

* `age_at_maturity` is replaced by `ages_movement`. The former can take either a
  single integer indicating the age at which individuals start to move, or a
  vector with 0s for age-groups that do not move and 1s for age-groups that are
  allowed to move.

# drmr 0.0.23

* The prior for $\alpha$ has been modified. In particular, instead of a pcp
  prior, now we put a Beta prior on $(\alpha + 1) / 2$. The hyperparameters of
  this Beta distribution are `pr_alpha_a` and `pr_alpha_b`, respectively.
   
* New functions to initialize the parameters from the prior have been
  introduced. Now the `init` parmeter from the `fit_sdm` and `fit_drm` functions
  can take three possible values:
  * "default": the standard initialization in `Stan` (For details see [`Stan's
    documentation`](https://mc-stan.org/docs/reference-manual/execution.html#random-initial-values)).
  * "prior": initialize the parameters using samples from their respective prior
    distributions;
  * "pathfinder": uses the
    [Pathfinder](https://mc-stan.org/docs/reference-manual/pathfinder.html)
    algorithm to initialize the parameters.

* Fixed some issues with documentation.

# drmr 0.0.22

* `fit_drm` and `fit_sdm` functions to make model fitting slightly simpler.

* `pr_phi_a` and `pr_phi_b` become `pr_phi_mu` and `pr_phi_sd`. The prior for
  $\log(\phi)$ is a Student's t with 3 degrees of freedom, mean `pr_phi_mu` and
  SD `pr_phi_sd`.
   
* Selectivity in `make_data` was fixed, thanks to Mark. Before, it was not being
  used when users input it.

# drmr 0.0.21

* `make_data_sdm` function (analogous to `make_data`) created for SDM.

* `predict_sdm` function (analogous to `predict_drm`) created for SDM.

# drmr 0.0.2

* Parameters in the code and documentation were properly matched.

* `p_error` toggle becomes `time_ar` toggle (more appropriate).

* `predict_drm` function created

* the `make_data` function now has a `family` argument indicating the
  probability distribution assumed for the response (given all the model
  parameters and latent variables)

# drmr 0.0.1

* First version of the package containing pre-compiled code for fitting the DRM,
  making forecasts and calculating density weighted centroids.
