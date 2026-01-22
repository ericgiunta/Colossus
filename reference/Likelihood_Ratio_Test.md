# Defines the likelihood ratio test

`Likelihood_Ratio_Test` uses two models and calculates the ratio

## Usage

``` r
Likelihood_Ratio_Test(alternative_model, null_model)
```

## Arguments

- alternative_model:

  the new model of interest in list form, output from a Poisson
  regression

- null_model:

  a model to compare against, in list form

## Value

returns the score statistic

## Examples

``` r
library(data.table)
# In an actual example, one would run two seperate RunCoxRegression regressions,
#    assigning the results to e0 and e1
a <- c(0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6)
b <- c(1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7)
c <- c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
d <- c(3, 4, 5, 6, 7, 8, 9, 1, 2, 1, 1, 2, 1, 2)
e <- c(1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2)
df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
keep_constant <- c(0)
a_n <- c(-0.1, 0.1, 0.1, 0.2)
control <- list("ncores" = 1, "maxiter" = 10, "verbose" = 0)
alternative_model <- CoxRun(Cox(a, b, c) ~ plinear(d * d, 0) + loglinear(factor(e)), df, control = control, a_n = a_n, norm = "max", keep_constant = c(0, 1, 0))
null_model <- CoxRun(Cox(a, b, c) ~ null(), df, control = control)
score <- Likelihood_Ratio_Test(alternative_model, null_model)
```
