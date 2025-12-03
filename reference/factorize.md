# Splits a parameter into factors

`factorize` uses user provided list of columns to define new parameter
for each unique value and update the data.table. Not for interaction
terms

## Usage

``` r
factorize(df, col_list, verbose = 0)
```

## Arguments

- df:

  a data.table containing the columns of interest

- col_list:

  an array of column names that should have factor terms defined

- verbose:

  integer valued 0-4 controlling what information is printed to the
  terminal. Each level includes the lower levels. 0: silent, 1: errors
  printed, 2: warnings printed, 3: notes printed, 4: debug information
  printed. Errors are situations that stop the regression, warnings are
  situations that assume default values that the user might not have
  intended, notes provide information on regression progress, and debug
  prints out C++ progress and intermediate results. The default level is
  2 and True/False is converted to 3/0.

## Value

returns a list with two named fields. df for the updated dataframe, and
cols for the new column names

## See also

Other Data Cleaning Functions: [`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
a <- c(0, 1, 2, 3, 4, 5, 6)
b <- c(1, 2, 3, 4, 5, 6, 7)
c <- c(0, 1, 2, 1, 0, 1, 0)
df <- data.table::data.table("a" = a, "b" = b, "c" = c)
col_list <- c("c")
val <- factorize(df, col_list)
df <- val$df
new_col <- val$cols
```
