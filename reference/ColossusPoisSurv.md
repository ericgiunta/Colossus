# Interprets basic poisson survival formula RHS

`ColossusPoisSurv` assigns and interprets interval columns for poisson
model. This functions is called using the arguments for Poisson or
Poisson_Strata in the right-hand side of the formula. Uses an
person-year column, number of events, and any strata columns. The first
two are expected to be in order or named: pyr and event. Anything beyond
the event name is assumed to be strata if Poisson_Strata is used.

## Usage

``` r
ColossusPoisSurv(...)
```

## Arguments

- ...:

  entries for a Poisson object with or without strata, pyr, event, and
  any strata columns. Either in order or named. The first two are
  assumed to be pyr and event, the rest assumed to be strata columns

## Value

returns list with duration, strata if used, and event

## See also

Other Formula Interpretation: [`ColossusCoxSurv()`](ColossusCoxSurv.md),
[`ColossusLogitSurv()`](ColossusLogitSurv.md),
[`get_form()`](get_form.md), [`get_form_joint()`](get_form_joint.md)
