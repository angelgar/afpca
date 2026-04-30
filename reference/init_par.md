# Get initial values for parameters (scores estimated from FPCA.FACE and coefficients from OLS)

Internal function to get scores given all other parameters

## Usage

``` r
init_par(data, n.comp, Theta, orthogonalize_fpcs, N.Unp.Basis)
```

## Arguments

- data:

  input data

- n.comp:

  number of fpcs being calculated

- Theta:

  spline basis matrix

- orthogonalize_fpcs:

  wether or not to orthogonalize estimated fpcs

- N.Unp.Basis:

  number of unpenalized spline basis
