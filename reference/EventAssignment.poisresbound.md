# Predicts how many events are due to baseline vs excess for a completed poisson likelihood boundary regression

`EventAssignment.poisresbound` uses user provided data,
person-year/event columns, vectors specifying the model, and options to
calculate background and excess events for a solved Poisson regression

## Usage

``` r
# S3 method for class 'poisresbound'
EventAssignment(
  x,
  df,
  assign_control = list(),
  control = list(),
  a_n = c(),
  ...
)
```

## Arguments

- x:

  result object from a regression, class poisres

- df:

  a data.table containing the columns of interest

- assign_control:

  control list for bounds calculated

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- ...:

  can include the named entries for the assign_control list parameter

## Value

returns a list of the final results

## See also

Other Poisson Wrapper Functions:
[`EventAssignment.poisres()`](EventAssignment.poisres.md),
[`LikelihoodBound.poisres()`](LikelihoodBound.poisres.md),
[`PoisRun()`](PoisRun.md), [`PoisRunJoint()`](PoisRunJoint.md),
[`PoisRunMulti()`](PoisRunMulti.md),
[`Residual.poisres()`](Residual.poisres.md)
