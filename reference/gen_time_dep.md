# Applies time dependence to parameters

`gen_time_dep` generates a new dataframe with time dependent covariates
by applying a grid in time

## Usage

``` r
gen_time_dep(
  df,
  time1,
  time2,
  event0,
  iscox,
  dt,
  new_names,
  dep_cols,
  func_form,
  fname,
  tform,
  nthreads = as.numeric(detectCores())
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

- iscox:

  boolean if rows not at event times should not be kept, rows are
  removed if true. a Cox proportional hazards model does not use rows
  with intervals not containing event times

- dt:

  spacing in time for new rows

- new_names:

  list of new names to use instead of default, default used if entry is
  ”"

- dep_cols:

  columns that are not needed in the new dataframe

- func_form:

  vector of functions to apply to each time-dependent covariate. Of the
  form func(df, time) returning a vector of the new column value

- fname:

  filename used for new dataframe

- tform:

  list of string function identifiers, used for linear/step

- nthreads:

  number of threads to use, do not use more threads than available on
  your machine

## Value

returns the updated dataframe

## See also

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md)

## Examples

``` r
library(data.table)
# Adapted from the tests
a <- c(20, 20, 5, 10, 15)
b <- c(1, 2, 1, 1, 2)
c <- c(0, 0, 1, 1, 1)
df <- data.table::data.table("a" = a, "b" = b, "c" = c)
time1 <- "%trunc%"
time2 <- "a"
event <- "c"
control <- list(
  "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9,
  "deriv_epsilon" = 1e-9, "step_max" = 1.0,
  "thres_step_max" = 100.0,
  "verbose" = FALSE, "ties" = "breslow", "double_step" = 1
)
grt_f <- function(df, time_col) {
  return((df[, "b"] * df[, get(time_col)])[[1]])
}
func_form <- c("lin")
df_new <- gen_time_dep(
  df, time1, time2, event, TRUE, 0.01, c("grt"), c(),
  c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 2
)
file.remove("test_new.csv")
#> [1] TRUE
```
