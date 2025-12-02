# Checks the default value for a given model, if every parameter were 0

`Check_Strata_Model` checks if a model is valid for stratified poisson

## Usage

``` r
Check_Strata_Model(term_n, tform, modelform, gmix_term, gmix_theta)
```

## Arguments

- term_n:

  term numbers for each element of the model

- tform:

  list of string function identifiers, used for linear/step

- modelform:

  string specifying the model type: M, ME, A, PA, PAE, GMIX, GMIX-R,
  GMIX-E

- gmix_term:

  binary vector to denote excess (1) and relative terms (0). Excess
  terms have 1 added before use in risk model

- gmix_theta:

  double, used in gmix model. Multiplicative combination is taken to
  power of theta, additive combination is taken to 1-theta.

## Value

TRUE if passed

## See also

Other Data Cleaning Functions: [`Check_Iters()`](Check_Iters.md),
[`Date_Shift()`](Date_Shift.md),
[`Event_Count_Gen()`](Event_Count_Gen.md),
[`Event_Time_Gen()`](Event_Time_Gen.md),
[`Joint_Multiple_Events()`](Joint_Multiple_Events.md),
[`Replace_Missing()`](Replace_Missing.md),
[`Time_Since()`](Time_Since.md), [`apply_norm()`](apply_norm.md),
[`factorize()`](factorize.md), [`gen_time_dep()`](gen_time_dep.md)
