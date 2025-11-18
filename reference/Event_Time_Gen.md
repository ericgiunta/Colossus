# uses a table, list of categories, list of summaries, list of events, and person-year information to generate person-time tables

`Event_Time_Gen` generates event-time tables

## Usage

``` r
Event_Time_Gen(table, pyr, categ, summaries, events, verbose = FALSE)
```

## Arguments

- table:

  dataframe with every category/event column needed

- pyr:

  list with entry and exit lists, containing day/month/year columns in
  the table

- categ:

  list with category columns and methods, methods can be either strings
  or lists of boundaries, includes a time category or entry/exit are
  both required for the pyr list

- summaries:

  list of columns to summarize, supports counts, means, and weighted
  means by person-year and renaming the summary column

- events:

  list of events or interests, checks if events are within each time
  interval

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

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Check_Strata_Model()`](Check_Strata_Model.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
a <- c(0, 1, 2, 3, 4, 5, 6)
b <- c(1, 2, 3, 4, 5, 6, 7)
c <- c(0, 1, 0, 0, 0, 1, 0)
d <- c(1, 2, 3, 4, 5, 6, 7)
e <- c(2, 3, 4, 5, 6, 7, 8)
f <- c(
  1900, 1900, 1900, 1900,
  1900, 1900, 1900
)
g <- c(1, 2, 3, 4, 5, 6, 7)
h <- c(2, 3, 4, 5, 6, 7, 8)
i <- c(
  1901, 1902, 1903, 1904,
  1905, 1906, 1907
)
table <- data.table::data.table(
  "a" = a, "b" = b, "c" = c,
  "d" = d, "e" = e, "f" = f,
  "g" = g, "h" = h, "i" = i
)
categ <- list(
  "a" = "-1/3/5]7",
  "b" = list(
    lower = c(-1, 3, 6), upper = c(3, 6, 10),
    name = c("low", "medium", "high")
  ),
  "time AS time" = list(
    "day" = c(1, 1, 1, 1, 1),
    "month" = c(1, 1, 1, 1, 1),
    "year" = c(1899, 1903, 1910)
  )
)
summary <- list(
  "c" = "count AS cases",
  "a" = "mean",
  "b" = "weighted_mean"
)
events <- list("c")
pyr <- list(
  entry = list(year = "f", month = "e", day = "d"),
  exit = list(year = "i", month = "h", day = "g"),
  unit = "years"
)
e <- Event_Time_Gen(table, pyr, categ, summary, events, T)
```
