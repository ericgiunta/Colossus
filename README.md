
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Colossus

<!-- badges: start -->

[![](https://img.shields.io/github/languages/code-size/ericgiunta/Colossus.svg)](https://github.com/ericgiunta/Colossus)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/ericgiunta/Colossus/graph/badge.svg?token=NMH5R502W8)](https://app.codecov.io/gh/ericgiunta/Colossus)
[![pkgdown](https://github.com/ericgiunta/Colossus/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ericgiunta/Colossus/actions/workflows/pkgdown.yaml)
[![OS_Checks](https://github.com/ericgiunta/Colossus/actions/workflows/OS_TEST.yml/badge.svg?branch=main)](https://github.com/ericgiunta/Colossus/actions/workflows/OS_TEST.yml)
[![](https://cranlogs.r-pkg.org/badges/grand-total/Colossus)](https://CRAN.R-project.org/package=Colossus)
<!-- badges: end -->

The goal of Colossus is to provide an open-source means of performing
survival analysis on big data with complex risk formulas. Colossus is
designed to perform Cox Proportional Hazard regressions and Poisson
regressions on datasets loaded as data.tables or data.frames. The risk
models allowed are sums or products of linear, log-linear, or several
other radiation dose response formulas highlighted in the vignettes.
Additional plotting capabilities are available.

By default, a fully portable version of the code is compiled, which does
not support OpenMP on every system. Note that Colossus requires OpenMP
support to perform parallel calculations. The environment variable
“R_COLOSSUS_NOT_CRAN” is checked to determine if OpenMP should be
disabled for linux compiling with clang. The number of cores is set to 1
if the environment variable is empty, the operating system is detected
as linux, and the default compiler or R compiler is clang. Colossus
testing checks for the “NOT_CRAN” variable to determine if additional
tests should be run. Setting “NOT_CRAN” to “false” will disable the
longer tests. Currently, OpenMP support is not configured for linux
compiling with clang.

Note: From versions 1.3.1 to 1.4.1 the expected inputs changed.
Regressions are now run with CoxRun and PoisRun and formula inputs.
Please see the “Unified Equation Representation” vignette for more
details.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(data.table)
library(parallel)
library(Colossus)
## basic example code reproduced from the starting-description vignette

df <- data.table(
  "UserID" = c(112, 114, 213, 214, 115, 116, 117),
  "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
  "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
  "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
  "a" = c(0, 1, 1, 0, 1, 0, 1),
  "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
  "c" = c(10, 11, 10, 11, 12, 9, 11),
  "d" = c(0, 0, 0, 1, 1, 1, 1)
)

model <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, 0) + linear(b, c, 1) + plinear(d, 2) + multiplicative()

a_n <- c(0.1, 0.1, 0.1, 0.1)

keep_constant <- c(0, 0, 0, 0)

control <- list(
  "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-9,
  "deriv_epsilon" = 1e-9, "step_max" = 1.0,
  "verbose" = 2, "ties" = "breslow"
)

e <- CoxRun(model, df, a_n = a_n, control = control)
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:         a  loglin           0         21.67085            NaN            NaN
#> 2:         b     lin           1          0.10000            NaN            NaN
#> 3:         c     lin           1          0.10000            NaN            NaN
#> 4:         d    plin           2          0.10000            Inf              1
#> 
#> Cox Model Used
#> -2*Log-Likelihood: 2.64,  AIC: 10.64
#> Iterations run: 27
#> maximum step size: 7.50e-01, maximum first derivative: 5.49e-10
#> Analysis converged
#> Run finished in 0.04 seconds
#> |-------------------------------------------------------------------|
```
