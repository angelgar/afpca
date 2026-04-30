# Main replacement for asp2()'s design matrix extraction. Internal function to generate knots

Main replacement for asp2()'s design matrix extraction. Internal
function to generate knots

## Usage

``` r
spline_design_matrices(x, knots, basis = "os", degree = 3)
```

## Arguments

- x:

  numeric predictor vector (your X_temp)

- knots:

  numeric knot vector (from default_knots() or user-supplied)

- basis:

  "trunc.poly" \| "tps" \| "os" (what you pass as 'basis' in f())

- degree:

  integer (or c(p,q) for "os") (your poly.degree)
