## package config
if (!file.exists("cleanup"))
  instantiate::stan_package_configure()

## licensing
if (!file.exists("LICENSE.md"))
  usethis::use_gpl_license()

## Rmd readme
if (!file.exists("README.Rmd"))
  usethis::use_readme_rmd()

## badge for experimental
usethis::use_lifecycle_badge("experimental")

## later (building a webpage for the pkg)
## usethis::use_pkgdown()

## creating a vignette
if (!file.exists("vignettes/get-started.Rmd"))
  usethis::use_vignette(name = "get-started")

## build pkgdown website
usethis::use_pkgdown_github_pages()
