# Estimate likelihood of the data

Internal function to estimate the likelihood of the data, used to check
algorithm converges

## Usage

``` r
get_likelihood(
  data,
  coef,
  Theta,
  c_mat,
  N.Unp.Basis,
  nbs,
  n.comp,
  lambda,
  sigma
)
```

## Arguments

- data:

  input data

- coef:

  spline coefficients

- Theta:

  spline basis

- c_mat:

  matrix of scores

- N.Unp.Basis:

  number of unpenalized spline basis

- nbs:

  total number of spline basis

- n.comp:

  number of components being estimated

- lambda:

  penalty function

- sigma:

  residual variance
