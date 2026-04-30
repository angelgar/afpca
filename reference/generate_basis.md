# Generate Spline Basis to estimate FPC

This is an internal function to generate a spline basis used to
calculate the FPCs using asp2 from AdaptFitOS

## Usage

``` r
generate_basis(
  data,
  basis = "trunc.poly",
  poly.degree = 3,
  knots = NA,
  nbs = 40
)
```

## Arguments

- data:

  input data

- basis:

  Type of spline basis. "os" for B-spline basis and "trunc.poly" for
  truncated polynomials (default)

- poly.degree:

  Degree of spline basis

- knots:

  passed to AdaptFitOS

- nbs:

  number of spline basis
