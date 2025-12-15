# Calculates hazard ratios for a reference vector

`coxres.RelativeRisk` uses a cox result object and data, to evaluate
relative risk in the data using the risk model from the result

## Usage

``` r
# S3 method for class 'coxres'
RelativeRisk(x, df, a_n = c(), ...)
```

## Arguments

- x:

  result object from a regression, class coxres

- df:

  a data.table containing the columns of interest

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- ...:

  extended to match any future parameters needed

## Value

returns a class fully describing the model and the regression results

## Examples

``` r
library(data.table)
df <- data.table::data.table(
  "UserID" = c(112, 114, 213, 214, 115, 116, 117),
  "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
  "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
  "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
  "a" = c(0, 1, 1, 0, 1, 0, 1),
  "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
  "c" = c(10, 11, 10, 11, 12, 9, 11),
  "d" = c(0, 0, 0, 1, 1, 1, 1),
  "e" = c(0, 0, 1, 0, 0, 0, 1)
)
control <- list(
  "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
  "halfmax" = 1
)
formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
  loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
res <- CoxRun(formula, df,
  a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)),
  control = control
)
res_risk <- RelativeRisk(res, df)
```
