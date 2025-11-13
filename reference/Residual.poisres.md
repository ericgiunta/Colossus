# Calculates the Residuals for a completed poisson model

`Residual.poisres` uses user provided data, person-year/event columns,
vectors specifying the model, and options to calculate residuals for a
solved Poisson regression

## Usage

``` r
# S3 method for class 'poisres'
Residual(
  x,
  df,
  control = list(),
  a_n = c(),
  pearson = FALSE,
  deviance = FALSE,
  ...
)
```

## Arguments

- x:

  result object from a regression, class poisres

- df:

  a data.table containing the columns of interest

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- pearson:

  boolean to calculate pearson residuals

- deviance:

  boolean to calculate deviance residuals

- ...:

  can include the named entries for the assign_control list parameter

## Value

returns a list of the final results

## See also

Other Poisson Wrapper Functions:
[`EventAssignment.poisres()`](EventAssignment.poisres.md),
[`EventAssignment.poisresbound()`](EventAssignment.poisresbound.md),
[`LikelihoodBound.poisres()`](LikelihoodBound.poisres.md),
[`PoisRun()`](PoisRun.md), [`PoisRunJoint()`](PoisRunJoint.md),
[`PoisRunMulti()`](PoisRunMulti.md),
[`RunPoissonRegression_Omnibus()`](RunPoissonRegression_Omnibus.md)
