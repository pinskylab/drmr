# Calculate the interval score

This function calculates the interval score for a given set of
observations, lower and upper bounds, and alpha parameter.

## Usage

``` r
int_score(y, l, u, alpha)
```

## Arguments

- y:

  A numeric vector of observations.

- l:

  A numeric vector of lower bounds for the prediction intervals.

- u:

  A numeric vector of upper bounds for the prediction intervals.

- alpha:

  A numeric value specifying the significance level (e.g., 0.05 for a
  95% interval).

## Value

A numeric vector of interval scores.

## Details

The interval score is a proper scoring rule that measures the accuracy
of interval predictions. It takes into account both the coverage and the
width of the prediction interval. A lower score indicates a better
prediction.
