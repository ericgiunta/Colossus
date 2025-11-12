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
e0 <- list("name" = "First Model", "LogLik" = -120)
e1 <- list("name" = "New Model", "LogLik" = -100)
score <- Likelihood_Ratio_Test(e1, e0)
```
