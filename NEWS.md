# drmr 0.0.1

* First version of the package containing pre-compiled code for fitting the DRM,
  making forecasts and calculating density weighted centroids.

# drmr 0.0.2

* Parameters in the code and documentation were properly matched.

* `p_error` toggle becomes `time_ar` toggle (more appropriate).

* `predict_drm` function created

* the `make_data` function now has a `family` argument indicating the
  probability distribution assumed for the response (given all the model
  parameters and latent variables)

# drmr 0.0.21

* `make_data_sdm` function (analogous to `make_data`) created for SDM.

* `predict_sdm` function (analogous to `predict_drm`) created for SDM.

# drmr 0.0.22

* `fit_drm` and `fit_sdm` functions to make model fitting slightly simpler.

* `pr_phi_a` and `pr_phi_b` become `pr_phi_mu` and `pr_phi_sd`. The prior for
  $\log(\phi)$ is a Student's t with 3 degrees of freedom, mean `pr_phi_mu` and
  SD `pr_phi_sd`.
  
* Selectivity in `make_data` was fixed, thanks to Mark. Before, it was not being
  used when users input it.

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

# drmr 0.0.25

* Exporting `fix_linbeta`, `max_quad_x`, and `int_score` functions

* Fixing `check_between` function (not exported; This functions is just a
  helper)
  
* Changed prior on `zeta` and allowing for user to input the
  hyperparameters. Before, we had a standard normal prior on the logit of
  `zeta`. Now, we place a beta prior on `zeta`.

* Included support for different types of population dynamics initialization

* Included tools for prior sensitivity analysis

* Fixed how pdf and random number generation from location-scale Student's t
  distribution
