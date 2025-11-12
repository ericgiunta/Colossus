# Interprets basic logistic survival formula RHS with no grouping

`ColossusLogitSurv` assigns and interprets columns for trials and events
in logistic model with no grouping.

## Usage

``` r
ColossusLogitSurv(...)
```

## Arguments

- ...:

  entries for a Logistic object, trials and events. trials not provided
  assumes one trial per row.

## Value

returns list with event

## See also

Other Formula Interpretation: [`ColossusCoxSurv()`](ColossusCoxSurv.md),
[`ColossusPoisSurv()`](ColossusPoisSurv.md),
[`get_form()`](get_form.md), [`get_form_joint()`](get_form_joint.md)
