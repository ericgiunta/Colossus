# Automatically applies a normalization to either an input or output

`apply_norm` applies a normalization factor

## Usage

``` r
apply_norm(df, norm, names, input, values, model_control)
```

## Arguments

- df:

  The data.table with columns to be normalized

- norm:

  The normalization option used, currently max or mean

- names:

  columns for elements of the model, used to identify data columns

- input:

  boolean if the normalization is being performed on the input values or
  on an output

- values:

  list of values using during normalization

- model_control:

  controls which alternative model options are used, see the
  Control_Options vignette for further details

## Value

returns list with the normalized values

## See also

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Check_Strata_Model()`](Check_Strata_Model.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`factorize()`](factorize.md),
[`gen_time_dep()`](gen_time_dep.md)
