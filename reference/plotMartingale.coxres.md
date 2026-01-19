# Performs Cox Proportional Hazard model martingale residual plots

`plotMartingale.coxres` uses user provided data, time/event columns,
vectors specifying the model, and options to choose and save plots

## Usage

``` r
# S3 method for class 'coxres'
plotMartingale(x, df, plot_options, a_n = c(), ...)
```

## Arguments

- x:

  result object from a regression, class coxres

- df:

  a data.table containing the columns of interest

- plot_options:

  list of parameters controlling the plot options, see RunCoxPlots() for
  different options

- a_n:

  list of initial parameter values, used to determine the number of
  parameters. May be either a list of vectors or a single vector.

- ...:

  can include the named entries for the plot_options parameter

## Value

returns the data used for plots

## See also

Other Plotting Wrapper Functions: [`plot.coxres()`](plot.coxres.md),
[`plotRisk.coxres()`](plotRisk.coxres.md),
[`plotSchoenfeld.coxres()`](plotSchoenfeld.coxres.md),
[`plotSurvival.coxres()`](plotSurvival.coxres.md)
