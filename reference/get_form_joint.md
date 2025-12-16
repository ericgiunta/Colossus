# Interprets a Poisson joint formula and makes necessary changes to data

`get_form_joint` uses two event formula, a shared formula, and
data.table, to fully describe the model for a joint Poisson model.

## Usage

``` r
get_form_joint(formula_list, df, nthreads = as.numeric(detectCores())/2)
```

## Arguments

- formula_list:

  a list of formula objects, each written in Colossus notation. See the
  Unified Equation Representation vignette for details. Each formula
  should include the elements specific to the specified event column.
  The list can include an entry named "shared" to denote shared terms.
  The person-year and strata columns should be the same.

- df:

  a data.table containing the columns of interest

- nthreads:

  number of threads to use, do not use more threads than available on
  your machine

## Value

returns a class fully describing the model and the updated data

## See also

Other Formula Interpretation: [`ColossusCoxSurv()`](ColossusCoxSurv.md),
[`ColossusLogitSurv()`](ColossusLogitSurv.md),
[`ColossusPoisSurv()`](ColossusPoisSurv.md), [`get_form()`](get_form.md)
