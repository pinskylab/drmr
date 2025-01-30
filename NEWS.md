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

* `make_data_drm` function (analogous to `make_data`) created for SDM.
