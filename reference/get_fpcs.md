# Estimate mean function and eigenfunctions given a spline basis and coefficients

Internal function to estimate mean function and fpcs given spline basis
and coefficients

## Usage

``` r
get_fpcs(Theta, coef, nbs, n.comp, orthogonalize_fpcs, N.Unp.Basis)
```

## Arguments

- Theta:

  spline basis

- coef:

  spline coefficients coefficients

- nbs:

  total number of spline basis

- n.comp:

  number of functional principal components being estimated

- orthogonalize_fpcs:

  wether or not to orthogonalize FPCs post estimation

- N.Unp.Basis:

  number of unpenalized spline basis
