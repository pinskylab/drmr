# Default arguments for inference algorithm

Returns default for a given algorithm

## Usage

``` r
default_algo(algorithm = "nuts", algo_args = list())
```

## Arguments

- algorithm:

  a `character` specifying the algorithm used for inference. Default is
  `nuts` (the default MCMC in Stan). The remaining options are different
  flavors of variational bayes algorithms: "vb" (for ADVI), "pathfinder"
  (for Pathfinder), "laplace" (normal approximation centered at the mode
  of the posterior) or "optimize" for (penalized) MLEs.

- algo_args:

  a `list` with arguments for the sampling algorithms. For instance,
  `tol_rel_obj` for variational inference.

## Value

A `list` containing the default arguments for the specified algorithm.

## Author

lcgodoy
