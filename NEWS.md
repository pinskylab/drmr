# drmr (development version)

- New `areas` argument of `make_data()` (and so of `fit_drm()`) for the area
  of each site. The default is 1 for each site. A warning tells you when the
  model uses the areas but you did not give them.

- Movement now conserves abundance (density times area). Each moving age
  group goes through $\mathbf{A}^{-1} \mathbf{M}^\top \mathbf{A}$, where
  $\mathbf{M} = \zeta \mathbf{I} + (1 - \zeta) \mathbf{W}$, $\mathbf{W}$ is
  `adj_mat` with each row scaled to sum to 1, and $\mathbf{A}$ holds the site
  areas. Row $i$ of `adj_mat` now gives where individuals at site $i$ go.
  For an `adj_mat` that is not symmetric, this reverses the direction of
  earlier versions. Individuals at a site with no neighbors now stay there.
  Before, the share $1 - \zeta$ of them was lost every year.

- With `acc_dd = 1`, the regional SSB is now the area-weighted mean SSB
  density of the sites. Before, it was the sum of the SSB densities.

- Results for models with `movement = 1`, or with density dependence and
  `acc_dd = 1`, differ from those of earlier versions, even with equal
  areas. Refit these models. For fits from earlier versions, `predict()`,
  `ages_edens()`, and `get_adults()` still run, with an area of 1 for each
  site. But they apply the new movement and regional SSB to the old draws,
  so their results do not match the old fit.

- Changes to the ICAR random effects (`sp_re`). Results of models with
  `sp_re` change, so refit them.
  - Each pair of neighbors now enters the prior once. Before, a symmetric
    `adj_mat` listed each pair twice. So `sigma_s` was $\sqrt{2}$ times the
    typical marginal standard deviation of the spatial effect, and the prior
    pulled that SD toward values $\sqrt{2}$ times smaller than its
    hyperparameters said.
  - Graphs that are not connected now work, following Freni-Sterrantino et
    al. (2018). Each connected group of sites gets its own sum-to-zero
    constraint and scaling factor, and each site with no neighbors gets an
    independent effect. Before, the scaling factor was huge, so the spatial
    effect was close to zero. `gen_adj()` now accepts patches with no
    neighbors, such as islands.
  - The prior uses only which entries of `adj_mat` are positive, in either
    direction. Before, the scaling factor of a DRM also used the weights.
  - The effects now work when every site has the same number of neighbors,
    and when `adj_mat` has no row names. Before, they used wrong pairs of
    neighbors when the row names of `adj_mat` were not 1 to n in order.
  - The model parameter is now `w_s_raw`. `w_s` and `z_s` are transformed
    parameters.

- `fit_drm(init = "prior")` and `fit_sdm(init = "prior")` now work.
  `prior_sample()` drew initial values for random effects that the model did
  not have, used wrong names and dimensions, and failed for every DRM. The
  docs now say that, with `rho_mu = 1`, a chain can fail to start if the
  draw of `xi` is far in the tail of its prior.

- Clearer errors for wrong toggles and priors:
  - A misspelled toggle value, such as `rec_dd = "Ricker"`, or a value of
    `TRUE` or `FALSE`, now gives an error. Before, a misspelled value kept
    the default without any message.
  - An unknown toggle or prior name now gives an error that names it.
  - In an SDM, `ar_re`, `iid_re`, and `sp_re` must be 0 or 1.
  - Random effects on survival (`"surv"`) need `est_surv = 1`.
  - `adj_mat` must have finite, non-negative values when the model uses it.

- Removed `elpd()`. It could not run, because the package did not compile its
  Stan models.

- New `sim_drm()` simulates data from a DRM for known parameter values, for
  parameter-recovery and other simulation studies.

- drmr now requires R 4.3.0 or later.

# drmr 1.1.2

- Fixed bugs on data ordering and predictions.

- Using hyperbolic tangent transform for the temporal autocorrelation parameter.

# drmr 1.1.1

- Fixing some bugs in the parametrization of the loglogistic and
  truncated-Normal logliks.

# drmr 1.1.0

- Density dependence also support SSB for the whole region (as opposed to local
  SSB).

- added support for density dependence.

- Minor bug in the `log_lik` function was corrected.

- bug in the prediction function fixed.

# drmr 1.0.1

- Minor bug on the movement routine was corrected. The previous version was only
  smoothing the densities, and movement was not influencing the population
  dynamics.

- Speeding up movement through `Stan`'s [sparse matrix
  operations](https://mc-stan.org/docs/functions-reference/sparse_matrix_operations.html).

# drmr 1.0.0

* `ages_edens` now supports a `new_data` argument which allows for estimating
  out-of-sample age-specific densities.

* Fixed a bug on the `stan` function for the age-structured population dynamic.

* Changing the parametrization for survival (it does help with identifiability)

* `marg` method renamed to `effects_drm`.

* `plot` function for easily diagnosing MCMC convergence.

* A `mcmc_diag` function is also provided for assessing each parameters
  convergence. This function accompanies a `summary` method.
  
* A new `ppc` function is provided for some simple visual posterior predictive
  checks.

* Adding parameter description to the output of the `summary` method.

# drmr 0.5.0

* New vignette on "algorithms" detailing the methods that can be used to draw
  samples from the posterior.

* `summary` method for the output of the `ages_edens` output.

* `draws` method added to the package.

* `marg` method to .

* Fixed a bug in the `sdm.stan` AR(1) term.

* Fixed duplicated random effects in `drm.stan`

* Improved performance for `drm.stan` and `sdm.stan`

* Super class for DRMs and SDMs, so `print`, `summary`, and `draws` methods are
  simplified.

* A class `pred_drmr` was created for the output of the `predict` and `fitted`
  methods. It facilitates joining those outputs to the observed data.

* We had some functions called patches sites and vice-versa. Now everything is
  standardized using "site" for patches/sites (I know, it does sound confusing).

* `elpd` functions temporarily unavailable 

# drmr 0.4.0

* `predict` and `fitted` methods now have an argument called `type`, which
  allows to compute predictions based on the posterior predictive distribution
  (default and corresponds to old behavior), expected value of the posterior
  predictive distribution (i.e., its theoretical mean), and latent density
  (disconsidering observation error).

* Removing loglogistic distribution for now (need to do some work on the
  parametrization)

* Major bug in the fitted/predict functions fixed.

* Some improvements to documentation.

* Fixed bug in how we deal with zeros and NAs.

* A function called `summary_age_densities` to process the output of
  `ages_edens` is now available.

* Fixed (again) the `lambda_drm` function (to extract age specific densities)
  and changed its name to `ages_edens` (stands for "age-specific expected age
  densities".

* `log_lik` and `fitted` methods are now defined for `adrm` and `sdm` objects.

* Calculations for log-likelihood improved for easier maintenance.

* `predict_drm` and `predict_sdm` are now methods.

* Another increase in performance for the `predict` (and `elpd`)
  function. Especially relevant for larger datasets.

* `print` and `summary` methods are now available for `adrm` and `sdm` objects.

* Different inference algorithms are now alowed through the `algorithm`
  argument. A new `algo_args` argument is also introduced for specific inputs
  regarding the `algorithms`. For example, `iter_warmup` (specifying the number
  of warmup iterations), `iter_sampling` (number of sampling iterations),
  `chains` (number of MCMC chains), `parallel_chains` (number of chains to run
  in parallel), `adapt_delta` among others are supposed to be input as elements
  of a list passed as `algo_args` parameter.

# drmr 0.3.0

* New classes (and validators) created for the output of the `fit_drm` and
  `fit_sdm` functions. The classes are called `adrm` and `sdm`, respectively.

* An `update` method to facilitate model refitting is introduced.

# drmr 0.2.3

* Bug fixed when computing the log-lik from a truncated normal model.

* The `elpd` method has been included with the package and works with the
  outputs of `fit_drm` and `fit_sdm`. This method is helpful in calculating
  out-of-sample ELPD to perform leave-future-out cross validation (e.g., see
  https://mc-stan.org/loo/articles/loo2-lfo.html )

* The `predict_*` functions have been modified so they only return the predicted
  values at the new sites/times. Before this update, those functions also
  returned other quantities such as `mu_proj` (the expected value of the
  "non-zero" part of the density), `rho_proj` (probability of non-zero density),
  and `zt` (the AR random effect, if its switch is turned on).

* Improved the code to facilitate manteinance of the likelihood functions

* Fixing `lambda_gq`.

* Creating `fitted` and `log_likelihood` functions.

# drmr 0.2.2

* Ensuring there are no `NA`s in the covariates, `site`, or `time`.

* Added support for missingness on `y`.

* Functions for vectorzing the zero-inflated densities were moved to `R` (as
  opposed to `Stan`). Therefore, the indexes of the "zeros" in the response
  variables become an input from `R` to `Stan`.

* Spatial and iid random effects also included in the SDM functions.

* Typos in the "Parameters, priors, and toggles" vignette were fixed.

* `instantiate` version must be at least `0.2.3.9002`, so we can use `#include`
  for functions.

# drmr 0.2.1

* A new toggle called `rho_mu` introduced. Its default value is `1` and stands
  for explicitly relating the probability of observing a 0 (`rho`) to the latent
  density (`mu`). See [Yee 2014](https://doi.org/10.1016/j.csda.2013.01.012) and
  references therein.

* A new vignette with a comprehensive list of the model parameters and their
  priors is available.

* A new vignette detailing the initialization procedures is also included.

# drmr 0.2.0

* `init` now admits real numbers too.

* `init_data`: default now is to initialize the age-classes for all patches
  ranging from `.9` to `.01`.

* Allowing for choosing which process is correlated in time;

  * The variable `raw` became `w_t`;
  * `tau` became `sigma_t`
  * The flag `time_ar` became `ar_re` (stands for AR random effect). Now,
    `ar_re` can assume three values: 
    - "none", the default indicating no AR random effects.
    - "rec" (AR for recruitment)
    - "surv" (AR for survival)
    - "dens" (AR for density)

* Unstructured random effects.
  
  * The following variables are added to the model: `z_i` patch specific random
    effect and `sigma_i` the SD of the iid random effect.
  * Similarly to `ar_re`, the unstructured random effects can be enabled through
    the flat `iid_re`, which also admits the following entries:
    - "none" (default) indicating no IID random effects.
    - "rec" (IID for recruitment)
    - "surv" (IID for survival)
    - "dens" (iID for density)

* ICAR random effects
  
  * The following variables are added to the model: `w_s`, `z_s` patch specific
    random effect and `sigma_s` the approx marginal SD of the spatial random
    effect.
  * Similarly to the two structures mentioned above, the flag `sp_re` admits the
    following values: "none" (default), "rec" (recruitment), "surv" (survival),
    and "dens" (density).

* Functions for data simulation were temporarily removed.

* Fixing `init_data`

* Fixing bugs regarding the popdyn initialization.

# drmr 0.1.3

* Modified `examples.qmd` to prevent errors due to pathfinder failing.

* Updated github actions.

* Same problematic behaviour that was affecting `make_data` had to be dealed
  with prior initialization and predictions.
  
* Options for "data informed" phi priors through the floag `phi_hat`
  
  * When `family = "gamma"` and a prior for `phi` is not provided, the default
    prior for this parameter is a Gamma with `shape = 2` and `rate = ybar *
    s2`. This prior has its theoretical mean at the method of moments estimator
    of `phi`.

  * When `family = "lognormal"` and a prior for `phi` is not provided, the
    default prior for this parameter is a Gamma having its theoretical mean at
    the MLE of `phi`.

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
