# Generic Risk Plotting function, default option

`plotRisk.default` Generic Risk Plotting, by default nothing happens

## Usage

``` r
# Default S3 method
plotRisk(x, df, plot_options, a_n = c(), ...)
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
