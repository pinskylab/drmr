
<!-- README.md is generated from README.Rmd. Please edit that file -->

# drmr

## **D**ynamic **R**ange **M**odels in **R**

<!-- badges: start -->

[![R-CMD-check](https://github.com/pinskylab/drmr/workflows/check-cran/badge.svg)](https://github.com/pinskylab/drmr/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Installing the `drmr` package

The `drmr` package provides pre-compiled
[`cmdstanr`](https://mc-stan.org/cmdstanr/) models powered by the
[`instantiate`](https://github.com/wlandau/instantiate) package . These
pre-compiled models allow us to make inferences about the DRM (and SDM) using
`cmdstan` algorithms, such as the HMC NUTS Sampler.


The package is not on CRAN yet. To install the version hosted on GitHub, run:
```{r}
#| eval: false

remotes::install_github("pinskylab/drmr")
```
