# Automatically assigns missing values in listed columns

`Replace_Missing` checks each column and fills in NA values

## Usage

``` r
Replace_Missing(df, name_list, msv, verbose = FALSE)
```

## Arguments

- df:

  a data.table containing the columns of interest

- name_list:

  vector of string column names to check

- msv:

  value to replace na with, same used for every column used

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

returns a filled datatable

## See also

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
## basic example code reproduced from the starting-description vignette
df <- data.table::data.table(
  "UserID" = c(112, 114, 213, 214, 115, 116, 117),
  "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
  "Ending_Age" = c(30, 45, NA, 47, 36, NA, 55),
  "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0)
)
df <- Replace_Missing(df, c("Starting_Age", "Ending_Age"), 70)
```
