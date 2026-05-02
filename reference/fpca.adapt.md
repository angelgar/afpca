# Adaptive Functional Principal Component Analysis

Decomposes functional observations using adaptive functional principal
component analysis. This method incorporates adaptive smoothness into
the estimated decomposition by using a penalized likelihood framework

## Usage

``` r
fpca.adapt(
  data,
  col_name = "name",
  basis = "trunc.poly",
  poly.degree = 2,
  knots = NA,
  nbs = 35,
  n.comp = 18,
  pve = 0.99,
  normalize.scores = TRUE,
  orthogonalize_scores = TRUE,
  orthogonalize_fpcs = TRUE,
  ntimes = 200,
  convergence.thresh = 1e-04
)
```

## Arguments

- data:

  A data frame, matrix, or a tidyfun "tfd.reg" object of functional
  observations.

- col_name:

  Column name of the data. Required if Y is a tibble. The package
  expects the column to be a "dft" object from tidyfun

- basis:

  Type of spline basis. "os" for B-spline basis and "trunc.poly" for
  truncated polynomials (default)

- poly.degree:

  The degree of the the spline basis

- knots:

  Number of knots used to generate a spline basis to estimate the mean
  function and functional principal components

- nbs:

  Number of spline basis to generate mean function and functional
  principal components

- n.comp:

  prespecified number of principal components

- pve:

  proportion of variance explained; used to choose the number of
  functional principal components

- normalize.scores:

  normalize scores at each iteration (default is true, recommended)

- orthogonalize_scores:

  orthogonalize scores at each iteration (default is true, recommended)

- orthogonalize_fpcs:

  orthogonalize FPCs at each iteration (default is true, recommended)

- ntimes:

  Maximum number of iterations

- convergence.thresh:

  Threshold for convergence of penalized likelihood

## Value

The function returns an object of class "afpca" that contains the
following:

- Y:

  The observed data. Pass either a matrix, a "dft" object from
  tidyverse, or a tibble. If a tibble is passed, then "col_name" is
  required

- Y_hat:

  The reconstructed fitted values

- mean:

  a vector of the estimated mean function

- fpcs:

  a tp (number of timepoints) by npc (number of fpcs) matrix of the
  estimated eigenfunctions

- scores:

  an n (number of subjects) by npc (number of fpcs) matrix of estimated
  scores

- eigenvalues:

  estimated eigenvalues i.e. variances of fpc scores

- Basis:

  list containing penalized (Theta_beta) and penalized (Theta_b) spline
  basis

- Theta:

  concatenation of (Theta_beta) and penalized (Theta_b) spline basis

- coef:

  estimated coefficients, coef %\*% theta = fpcs

- convergence:

  logical, algorithm converged or not

- no.iter:

  number of iterations needed for convergence

## References

Paper

## Examples

``` r

## Using simulated data
sim_data <- simulate_adaptive_functional_data(n.tp = 100)
afpca_output <- fpca.adapt(data = sim_data, nbs = 10)

## Using firing rate data (uses tidyfun, install beforehand)

if (FALSE) {
  data(firing_rates_data)
  afpca_firing_output <- fpca.adapt(data = firing_rates_data, col_name = "activation", nbs = 10)
}

```
