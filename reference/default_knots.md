# Direct replacement for AdaptFitOS::default.knots() / SemiPar::default.knots() Internal function to generate knots

Direct replacement for AdaptFitOS::default.knots() /
SemiPar::default.knots() Internal function to generate knots

## Usage

``` r
default_knots(x, num.knots, knotchoice = "quantiles")
```

## Arguments

- x:

  numeric vector (the predictor)

- num.knots:

  integer, number of knots; auto-computed when omitted

- knotchoice:

  "quantiles" (default) or "equidistant"
