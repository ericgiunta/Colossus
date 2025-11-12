# Interprets basic cox survival formula RHS

`ColossusCoxSurv` assigns and interprets interval columns for cox model.
This functions is called using the arguments for Cox in the right-hand
side of the formula. Uses an interval start time, end time, and event
status. These are expected to be in order or named: tstart, tend, and
event. The Fine-Gray and Stratified versions use strata and weight named
options or the last two entries.

## Usage

``` r
ColossusCoxSurv(...)
```

## Arguments

- ...:

  entries for a cox survival object, tstart, tend, and event. Either in
  order or named. If unnamed and two entries, tend and event are
  assumed.

## Value

returns list with interval endpoints and event

## See also

Other Formula Interpretation:
[`ColossusLogitSurv()`](ColossusLogitSurv.md),
[`ColossusPoisSurv()`](ColossusPoisSurv.md),
[`get_form()`](get_form.md), [`get_form_joint()`](get_form_joint.md)
