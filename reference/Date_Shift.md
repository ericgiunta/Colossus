# Automates creating a date difference column

`Date_Shift` generates a new dataframe with a column containing time
difference in a given unit

## Usage

``` r
Date_Shift(df, dcol0, dcol1, col_name, units = "days")
```

## Arguments

- df:

  a data.table containing the columns of interest

- dcol0:

  list of starting month, day, and year

- dcol1:

  list of ending month, day, and year

- col_name:

  vector of new column names

- units:

  time unit to use

## Value

returns the updated dataframe

## See also

Other Data Cleaning Functions:
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`factorize()`](factorize.md),
[`gen_time_dep()`](gen_time_dep.md)

## Examples

``` r
library(data.table)
m0 <- c(1, 1, 2, 2)
m1 <- c(2, 2, 3, 3)
d0 <- c(1, 2, 3, 4)
d1 <- c(6, 7, 8, 9)
y0 <- c(1990, 1991, 1997, 1998)
y1 <- c(2001, 2003, 2005, 2006)
df <- data.table::data.table("m0" = m0, "m1" = m1, "d0" = d0, "d1" = d1, "y0" = y0, "y1" = y1)
df <- Date_Shift(df, c("m0", "d0", "y0"), c("m1", "d1", "y1"), "date_since")
```
