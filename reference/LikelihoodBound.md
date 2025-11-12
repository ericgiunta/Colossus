# Generic likelihood boundary calculation function

`LikelihoodBound` Generic likelihood boundary calculation function

## Usage

``` r
LikelihoodBound(x, df, curve_control = list(), control = list(), ...)
```

## Arguments

- x:

  result object from a regression, class coxres or poisres

- df:

  a data.table containing the columns of interest

- curve_control:

  a list of control options for the likelihood boundary regression. See
  the Control_Options vignette for details.

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- ...:

  extended for other necessary parameters
