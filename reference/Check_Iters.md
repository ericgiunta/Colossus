# Automatically checks the number of starting guesses

`Check_Iters` checks the number of iterations and number of guesses, and
corrects

## Usage

``` r
Check_Iters(control, a_n)
```

## Arguments

- control:

  list of parameters controlling the convergence, see the
  Control_Options vignette for details

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

## Value

returns a list with the corrected control list and a_n

## See also

Other Data Cleaning Functions: [`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)
