# Fully runs a joint poisson regression model, returning the model and results

`PoisRunJoint` uses a list of formula, data.table, and list of controls
to prepare and run a Colossus poisson regression function on a joint
dataset

## Usage

``` r
PoisRunJoint(
  model,
  df,
  a_n = list(c(0)),
  keep_constant = c(0),
  control = list(),
  gradient_control = list(),
  single = FALSE,
  observed_info = FALSE,
  cons_mat = as.matrix(c(0)),
  cons_vec = c(0),
  norm = "null",
  ...
)
```

## Arguments

- model:

  either a formula written for the get_form function, or the model
  result from the get_form function.

- df:

  a data.table containing the columns of interest

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- keep_constant:

  binary values to denote which parameters to change

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- gradient_control:

  a list of control options for the gradient descent algorithm. If any
  value is given, a gradient descent algorithm is used instead of
  Newton-Raphson. See the Control_Options vignette for details

- single:

  a boolean to denote that only the log-likelihood should be calculated
  and returned, no derivatives or iterations

- observed_info:

  a boolean to denote that the observed information matrix should be
  used to calculate the standard error for parameters, not the expected
  information matrix

- cons_mat:

  Matrix containing coefficients for a system of linear constraints,
  formatted as matrix

- cons_vec:

  Vector containing constants for a system of linear constraints,
  formatted as vector

- norm:

  methods used to normalize the covariates. Default is 'null' for no
  normalization. Other options include 'max' to normalize by the
  absolute maximum and 'mean' to normalize by the mean

- ...:

  can include the named entries for the control list parameter

## Value

returns a class fully describing the model and the regression results

## See also

Other Poisson Wrapper Functions:
[`EventAssignment.poisres()`](EventAssignment.poisres.md),
[`EventAssignment.poisresbound()`](EventAssignment.poisresbound.md),
[`LikelihoodBound.poisres()`](LikelihoodBound.poisres.md),
[`PoisRun()`](PoisRun.md), [`PoisRunMulti()`](PoisRunMulti.md),
[`Residual.poisres()`](Residual.poisres.md),
[`RunPoissonRegression_Omnibus()`](RunPoissonRegression_Omnibus.md)

## Examples

``` r
library(data.table)
df <- data.table::data.table(
  "UserID" = c(112, 114, 213, 214, 115, 116, 117),
  "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
  "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
  "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
  "Flu_Status" = c(0, 1, 0, 0, 1, 0, 1),
  "a" = c(0, 1, 1, 0, 1, 0, 1),
  "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
  "c" = c(10, 11, 10, 11, 12, 9, 11),
  "d" = c(0, 0, 0, 1, 1, 1, 1),
  "e" = c(0, 0, 1, 0, 0, 0, 1)
)
control <- list(
  "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
  "halfmax" = 1
)
formula_list <- list(Pois(Ending_Age, Cancer_Status) ~ plinear(d, 0),
  Pois(Ending_Age, Flu_Status) ~ loglinear(d, 0),
  "shared" = Pois(Ending_Age) ~ loglinear(a, b, c, 0)
)
res <- PoisRunJoint(formula_list, df, control = control)
```
