# O'Sullivan spline decomposition into fixed (C) and random (Z) parts. Mirrors Predict.matrix.lme() / smooth.construct.os.smooth.spec() from asp-internal.r.

O'Sullivan spline decomposition into fixed (C) and random (Z) parts.
Mirrors Predict.matrix.lme() / smooth.construct.os.smooth.spec() from
asp-internal.r.

## Usage

``` r
os_design_matrices(x, nk, degree)
```

## Arguments

- x:

  predictor vector (length n)

- nk:

  number of interior knots (only the COUNT matters; knots are forced
  equidistant)

- degree:

  length-2 vector c(p, q) where p = B-spline degree, q = penalty order
