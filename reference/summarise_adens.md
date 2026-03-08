# Extract and Summarize Age-Specific Densities

Parses the output of `lambda_drm` to return a data frame of age-specific
densities, indexed by age, time, and patch.

## Usage

``` r
summarise_adens(lambda_obj, ages = NULL, probs = c(0.05, 0.5, 0.95))
```

## Arguments

- lambda_obj:

  A `CmdStanGQ` object returned by `lambda_drm`.

- ages:

  An optional vector of integers. If provided, the output is filtered to
  include only these specific ages.

- probs:

  A numeric vector of quantiles to calculate in the summary. Defaults to
  c(0.05, 0.5, 0.95).

## Value

A `data.frame` containing the summary statistics and parsed indices
(age, time, patch).
