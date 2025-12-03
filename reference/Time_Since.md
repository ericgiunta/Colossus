# Automates creating a date since a reference column

`Time_Since` generates a new dataframe with a column containing time
since a reference in a given unit

## Usage

``` r
Time_Since(df, dcol0, tref, col_name, units = "days")
```

## Arguments

- df:

  a data.table containing the columns of interest

- dcol0:

  list of ending month, day, and year

- tref:

  reference time in date format

- col_name:

  vector of new column names

- units:

  time unit to use

## Value

returns the updated dataframe

## See also

Other Data Cleaning Functions: [`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
m0 <- c(1, 1, 2, 2)
m1 <- c(2, 2, 3, 3)
d0 <- c(1, 2, 3, 4)
d1 <- c(6, 7, 8, 9)
y0 <- c(1990, 1991, 1997, 1998)
y1 <- c(2001, 2003, 2005, 2006)
df <- data.table::data.table(
  "m0" = m0, "m1" = m1,
  "d0" = d0, "d1" = d1,
  "y0" = y0, "y1" = y1
)
tref <- strptime("3-22-1997", format = "%m-%d-%Y", tz = "UTC")
df <- Time_Since(df, c("m1", "d1", "y1"), tref, "date_since")
```
