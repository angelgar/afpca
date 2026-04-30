# Simulate Adaptive Functional Data

Simulated functional data that includes with sharp changes over the
domain of the function (i.e. smoothness rapidly changes over the
domain). The data simulates data using sinusoidal functions with varying
period and amplitude.

## Usage

``` r
simulate_adaptive_functional_data(
  n.tp = 200,
  N.subj = 20,
  var.scores = c(4, 1),
  noise.var = 0.05,
  seed.num = 1
)
```

## Arguments

- n.tp:

  number of timepoints

- N.subj:

  number of subjects

- var.scores:

  variance of the scores (this is the value of the eigenvalues)

- noise.var:

  residual noise variance

- seed.num:

  seed number

## Value

An object of class "fpca_sim_data" that contains the following:

- data:

  simulated dataset

- data_true:

  simulated dataset without any noise

- Mu_true:

  The true data-generating mean-function

- Phi_true:

  The true data-generating functional principal components

- Scores_true:

  The true scores

## Examples

``` r
sim_data <- simulate_adaptive_functional_data(n.tp = 100, N.subj = 25)
```
