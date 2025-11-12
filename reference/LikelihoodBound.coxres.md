# Calculates the likelihood boundary for a completed cox model

`LikelihoodBound.coxres` solves the confidence interval for a cox model,
starting at the optimum point and iteratively optimizing end-points of
intervals.

## Usage

``` r
# S3 method for class 'coxres'
LikelihoodBound(x, df, curve_control = list(), control = list(), ...)
```

## Arguments

- x:

  result object from a regression, class coxres

- df:

  a data.table containing the columns of interest

- curve_control:

  a list of control options for the likelihood boundary regression. See
  the Control_Options vignette for details.

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- ...:

  can include the named entries for the curve_control list parameter

## Value

returns a list of the final results

## See also

Other Cox Wrapper Functions: [`CoxRun()`](CoxRun.md),
[`CoxRunMulti()`](CoxRunMulti.md),
[`RunCoxRegression_Omnibus()`](RunCoxRegression_Omnibus.md)
