# Performs Matched Case-Control Conditional Logistic Regression

`RunCaseControlRegression_Omnibus` uses user provided data, time/event
columns, vectors specifying the model, and options to control the
convergence and starting positions. Has additional options for starting
with several initial guesses, using stratification and/or matching by
time at risk, and calculation without derivatives

## Usage

``` r
RunCaseControlRegression_Omnibus(
  df,
  time1 = "%trunc%",
  time2 = "%trunc%",
  event0 = "event",
  names = c("CONST"),
  term_n = c(0),
  tform = "loglin",
  keep_constant = c(0),
  a_n = c(0),
  modelform = "M",
  control = list(),
  strat_col = "null",
  cens_weight = "null",
  model_control = list(),
  cons_mat = as.matrix(c(0)),
  cons_vec = c(0)
)
```

## Arguments

- df:

  a data.table containing the columns of interest

- time1:

  column used for time period starts

- time2:

  column used for time period end

- event0:

  column used for event status

- names:

  columns for elements of the model, used to identify data columns

- term_n:

  term numbers for each element of the model

- tform:

  list of string function identifiers, used for linear/step

- keep_constant:

  binary values to denote which parameters to change

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- modelform:

  string specifying the model type: M, ME, A, PA, PAE, GMIX, GMIX-R,
  GMIX-E

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- strat_col:

  column to stratify by if needed

- cens_weight:

  column containing the row weights

- model_control:

  controls which alternative model options are used, see the
  Control_Options vignette for further details

- cons_mat:

  Matrix containing coefficients for a system of linear constraints,
  formatted as matrix

- cons_vec:

  Vector containing constants for a system of linear constraints,
  formatted as vector

## Value

returns a list of the final results

## See also

Other Case Control Wrapper Functions:
[`CaseControlRun()`](CaseControlRun.md)

## Examples

``` r
library(data.table)
## basic example code reproduced from the starting-description vignette
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
# For the interval case
time1 <- "Starting_Age"
time2 <- "Ending_Age"
event <- "Cancer_Status"
names <- c("a", "b", "c", "d")
a_n <- list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4))
# used to test at a specific point
term_n <- c(0, 1, 1, 2)
tform <- c("loglin", "lin", "lin", "plin")
modelform <- "M"
keep_constant <- c(0, 0, 0, 0)
control <- list(
  "ncores" = 2, "lr" = 0.75, "maxiters" = c(5, 5, 5),
  "halfmax" = 5, "epsilon" = 1e-3, "deriv_epsilon" = 1e-3,
  "step_max" = 1.0, "thres_step_max" = 100.0,
  "verbose" = FALSE,
  "ties" = "breslow", "double_step" = 1
)
e <- RunCaseControlRegression_Omnibus(df, time1, time2, event,
  names, term_n, tform, keep_constant,
  a_n, modelform, control,
  model_control = list(
    "stata" = FALSE,
    "time_risk" = FALSE
  )
)
```
