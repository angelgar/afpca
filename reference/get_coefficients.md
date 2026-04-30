# Get the spline basis coefficients

Internal function to estimate the spline basis coefficients given the
scores and adaptive penalty

## Usage

``` r
get_coefficients(data, c_mat, Theta, Lambda, sigma, N.Unp.Basis, nbs)
```

## Arguments

- data:

  data input data

- c_mat:

  matrix of scores

- Theta:

  spline basis

- Lambda:

  penalty matrix

- sigma:

  residual variance

- N.Unp.Basis:

  number of unpelized spline basis

- nbs:

  total number of spline basis
