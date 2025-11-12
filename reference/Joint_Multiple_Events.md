# Automates creating data for a joint competing risks analysis

`Joint_Multiple_Events` generates input for a regression with multiple
non-independent events and models

## Usage

``` r
Joint_Multiple_Events(
  df,
  events,
  name_list,
  term_n_list = list(),
  tform_list = list(),
  keep_constant_list = list(),
  a_n_list = list()
)
```

## Arguments

- df:

  a data.table containing the columns of interest

- events:

  vector of event column names

- name_list:

  list of vectors for columns for event specific or shared model
  elements, required

- term_n_list:

  list of vectors for term numbers for event specific or shared model
  elements, defaults to term 0

- tform_list:

  list of vectors for subterm types for event specific or shared model
  elements, defaults to loglinear

- keep_constant_list:

  list of vectors for constant elements for event specific or shared
  model elements, defaults to free (0)

- a_n_list:

  list of vectors for parameter values for event specific or shared
  model elements, defaults to term 0

## Value

returns the updated dataframe and model inputs

## See also

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
a <- c(0, 0, 0, 1, 1, 1)
b <- c(1, 1, 1, 2, 2, 2)
c <- c(0, 1, 2, 2, 1, 0)
d <- c(1, 1, 0, 0, 1, 1)
e <- c(0, 1, 1, 1, 0, 0)
df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
time1 <- "t0"
time2 <- "t1"
df$pyr <- df$t1 - df$t0
pyr <- "pyr"
events <- c("e0", "e1")
names_e0 <- c("fac")
names_e1 <- c("fac")
names_shared <- c("t0", "t0")
term_n_e0 <- c(0)
term_n_e1 <- c(0)
term_n_shared <- c(0, 0)
tform_e0 <- c("loglin")
tform_e1 <- c("loglin")
tform_shared <- c("quad_slope", "loglin_top")
keep_constant_e0 <- c(0)
keep_constant_e1 <- c(0)
keep_constant_shared <- c(0, 0)
a_n_e0 <- c(-0.1)
a_n_e1 <- c(0.1)
a_n_shared <- c(0.001, -0.02)
name_list <- list("shared" = names_shared, "e0" = names_e0, "e1" = names_e1)
term_n_list <- list("shared" = term_n_shared, "e0" = term_n_e0, "e1" = term_n_e1)
tform_list <- list("shared" = tform_shared, "e0" = tform_e0, "e1" = tform_e1)
keep_constant_list <- list(
  "shared" = keep_constant_shared,
  "e0" = keep_constant_e0, "e1" = keep_constant_e1
)
a_n_list <- list("shared" = a_n_shared, "e0" = a_n_e0, "e1" = a_n_e1)
val <- Joint_Multiple_Events(
  df, events, name_list, term_n_list,
  tform_list, keep_constant_list, a_n_list
)
```
