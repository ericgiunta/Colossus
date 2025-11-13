# Calculates the likelihood boundary for a completed Poisson model

`LikelihoodBound.poisres` solves the confidence interval for a Poisson
model, starting at the optimum point and iteratively optimizing
end-points of intervals.

## Usage

``` r
# S3 method for class 'poisres'
LikelihoodBound(x, df, curve_control = list(), control = list(), ...)
```

## Arguments

- x:

  result object from a regression, class poisres

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

Other Poisson Wrapper Functions:
[`EventAssignment.poisres()`](EventAssignment.poisres.md),
[`EventAssignment.poisresbound()`](EventAssignment.poisresbound.md),
[`PoisRun()`](PoisRun.md), [`PoisRunJoint()`](PoisRunJoint.md),
[`PoisRunMulti()`](PoisRunMulti.md),
[`Residual.poisres()`](Residual.poisres.md),
[`RunPoissonRegression_Omnibus()`](RunPoissonRegression_Omnibus.md)
