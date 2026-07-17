
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Colossus

<!-- badges: start -->

[![](https://img.shields.io/github/languages/code-size/ericgiunta/Colossus.svg)](https://github.com/ericgiunta/Colossus)
[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/ericgiunta/Colossus/graph/badge.svg?token=NMH5R502W8)](https://app.codecov.io/gh/ericgiunta/Colossus)
[![pkgdown](https://github.com/ericgiunta/Colossus/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ericgiunta/Colossus/actions/workflows/pkgdown.yaml)
[![OS_Checks](https://github.com/ericgiunta/Colossus/actions/workflows/OS_TEST.yml/badge.svg?branch=main)](https://github.com/ericgiunta/Colossus/actions/workflows/OS_TEST.yml)
[![](https://cranlogs.r-pkg.org/badges/grand-total/Colossus)](https://CRAN.R-project.org/package=Colossus)
<!-- badges: end -->

The goal of `Colossus` is to provide an open-source means of performing
survival analysis on big data with complex risk formulas. `Colossus` is
designed to perform Cox Proportional Hazard regressions and Poisson
regressions on datasets loaded as data.tables or data.frames. The risk
models allowed are sums or products of linear, log-linear, or several
other radiation dose response formulas highlighted in the vignettes.
Additional plotting capabilities are available.

By default, a fully portable version of the code is compiled, which does
not support OpenMP on every system. Note that `Colossus` requires OpenMP
support to perform parallel calculations. The environment variable
“R_COLOSSUS_NOT_CRAN” is checked to determine if OpenMP should be
disabled for linux compiling with clang. The number of cores is set to 1
if the environment variable is empty, the operating system is detected
as linux, and the default compiler or R compiler is clang. `Colossus`
testing checks for the “NOT_CRAN” variable to determine if additional
tests should be run. Setting “NOT_CRAN” to “false” will disable the
longer tests. Currently, OpenMP support is not configured for linux
compiling with clang.

Note: From versions 1.3.1 to 1.4.1 the expected inputs changed.
Regressions are now run with [CoxRun()](reference/CoxRun.html) and
[PoisRun()](reference/PoisRun.html) and formula inputs. Please see the
[Unified Equation Representation
vignette](articles/Equation_Expression.html) for more details.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(data.table)
library(parallel)
library(Colossus)
## basic example code reproduced from the starting-description vignette
set.seed(3742)
df <- data.table(
  "UserID" = 1:100,
  "Starting_Age" = c(rep(10, 25), rep(15, 25), rep(20, 25), rep(25, 25)),
  "Ending_Age" = c(rep(16, 25), rep(21, 25), rep(26, 25), rep(31, 25)),
  "a" = rbinom(100, 3, 0.35),
  "b" = rbinom(100, 5, 0.25)
)

df$Cancer_Status <- vapply(df$a, function(x) rbinom(1, 1, x / 4), numeric(1))

model <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, 0) + plinear(b, 0) + multiplicative()

a_n <- c(0.1, 0.1)

control <- list(
  "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-9,
  "deriv_epsilon" = 1e-9, "step_max" = 1.0,
  "verbose" = 2, "ties" = "breslow"
)

e <- CoxRun(model, df, a_n = a_n, control = control)
print(e)
#> |--------------------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Central Estimate Standard Error 95% Confidence Interval
#>       <char>  <char>            <num>          <num>                  <char>
#> 1:         a  loglin            0.901         0.2560        (0.399 - 1.4029)
#> 2:         b    plin           -0.155         0.0976       (-0.347 - 0.0358)
#>    2-tail p-value
#>             <num>
#> 1:       0.000432
#> 2:       0.111185
#> |- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
#> 
#> Cox Model Used
#> Entry Age Column was: 'Starting_Age', Survival Age Column was: 'Ending_Age', Outcome Column was: 'Cancer_Status'
#> Risk Groups Used: 4
#> |- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
#> -2*Log-Likelihood: 182.74,  AIC: 186.74
#> Iterations run: 5
#> maximum step size: 1.133e-03, maximum first derivative: 5.375e-03
#> Last iteration improved the log-likelihood by: 1.500e-05
#> Analysis converged
#> Records Used: 100, Records Removed: 0
#> Run finished in 0.224 seconds
#> |--------------------------------------------------------------------------------|
```
