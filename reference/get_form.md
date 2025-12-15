# Interprets a Colossus formula and makes necessary changes to data

`get_form` uses a formula and data.table, to fully describe the model
for a Colossus regression function.

## Usage

``` r
get_form(formula, df, nthreads = as.numeric(detectCores())/2)
```

## Arguments

- formula:

  a formula object, written in Colossus notation. See the Unified
  Equation Representation vignette for details.

- df:

  a data.table containing the columns of interest

- nthreads:

  number of threads to use, do not use more threads than available on
  your machine

## Value

returns a class fully describing the model and the updated data

## See also

Other Formula Interpretation: [`ColossusCoxSurv()`](ColossusCoxSurv.md),
[`ColossusLogitSurv()`](ColossusLogitSurv.md),
[`ColossusPoisSurv()`](ColossusPoisSurv.md),
[`get_form_joint()`](get_form_joint.md)

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
formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
  loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
model <- get_form(formula, df, 1)
```
