# uses a table, list of categories, and list of event summaries to generate person-count tables

`Event_Count_Gen` generates event-count tables

## Usage

``` r
Event_Count_Gen(table, categ, events, verbose = FALSE)
```

## Arguments

- table:

  dataframe with every category/event column needed

- categ:

  list with category columns and methods, methods can be either strings
  or lists of boundaries

- events:

  list of columns to summarize, supports counts and means and renaming
  the summary column

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

returns a grouped table and a list of category boundaries used

## See also

Other Data Cleaning Functions: [`Date_Shift()`](Date_Shift.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`factorize()`](factorize.md),
[`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
a <- c(0, 1, 2, 3, 4, 5, 6)
b <- c(1, 2, 3, 4, 5, 6, 7)
c <- c(0, 1, 0, 0, 0, 1, 0)
table <- data.table::data.table(
  "a" = a,
  "b" = b,
  "c" = c
)
categ <- list(
  "a" = "0/3/5]7",
  "b" = list(
    lower = c(-1, 3, 6),
    upper = c(3, 6, 10),
    name = c("low", "medium", "high")
  )
)
event <- list(
  "c" = "count AS cases",
  "a" = "mean", "b" = "mean"
)
e <- Event_Count_Gen(table, categ, event, T)
```
