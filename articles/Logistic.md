# Logistic Regression

``` r
library(Colossus)
library(data.table)
library(survival)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:data.table':
#> 
#>     between, first, last
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Logistic Modeling for Binomial Odds and Binary Data

One alternative regression method in Colossus for predicting the
probability of an event in independent trials is logistic regression.
The theory is presented in the following section.

### General Theory

We start by describing the structure of the data. The logistic
regression models in Colossus assume the data can be represented in
terms of trials, events, and risk factors. Similar to a Poisson model
using person-years and number of events in each row, a logistic model
requires some measure of the number of independent trials and events for
each group. A special case is an analysis in which each row represents
one trial, and either 0 or 1 events are observed. The likelihood ($L$)
for a logistic model is a function of the trials ($N_{i}$), events
($y_{i}$), and the probability of a success ($p_{i}$) for each row, $i$.
In Colossus, the log-likelihood ($Ll$) is optimized.

$$\begin{array}{r}
{L = \prod\limits_{i}\left\lbrack p_{i}^{y_{i}} \times \left( 1 - p_{i} \right)^{N_{i} - y_{i}} \right\rbrack} \\
{Ll = \sum\limits_{i}\left\lbrack y_{i}\ln\left( p_{i} \right) + \left( N_{i} - y_{i} \right)\ln\left( 1 - p_{i} \right) \right\rbrack} \\

\end{array}$$

Similar to other survival models in Colossus, a model is written as a
function of the various risk factors
($f\left( \overset{\rightarrow}{\beta},\overset{\rightarrow}{x} \right)$).
This function is related to the probability through a linking function.
Three options are currently available: odds, identity, and complementary
log linking, based on the options present in 32-bit Epicure.

$$\begin{array}{r}
{p_{odds} = \frac{f\left( \overset{\rightarrow}{\beta},\overset{\rightarrow}{x} \right)}{1 + f\left( \overset{\rightarrow}{\beta},\overset{\rightarrow}{x} \right)}} \\
{p_{id} = f\left( \overset{\rightarrow}{\beta},\overset{\rightarrow}{x} \right)} \\
{p_{comp} = e^{- f{(\overset{\rightarrow}{\beta},\overset{\rightarrow}{x})}}}
\end{array}$$

Similar to other regression models in Colossus, the first derivative
vector and second derivative matrix are calculated to optimize the model
parameters.

$$\begin{array}{r}
{\frac{\partial Ll}{\partial\mu} = \sum\limits_{i}\left\lbrack y_{i}\frac{\partial p_{i}/\partial\mu}{p_{i}} - \left( N_{i} - y_{i} \right)\frac{\partial p_{i}/\partial\mu}{1 - p_{i}} \right\rbrack} \\
{\frac{\partial^{2}Ll}{\partial\mu\partial\nu} = \sum\limits_{i}\left\lbrack y_{i}\left( \frac{\partial^{2}p_{i}/\partial\mu\partial\nu}{p_{i}} - \frac{\partial p_{i}/\partial\mu}{p_{i}}\frac{\partial p_{i}/\partial\nu}{p_{i}} \right) - \left( N_{i} - y_{i} \right)\left( \frac{\partial^{2}p_{i}/\partial\mu\partial\nu}{1 - p_{i}} + \frac{\partial p_{i}/\partial\mu}{1 - p_{i}}\frac{\partial p_{i}/\partial\nu}{1 - p_{i}} \right) \right\rbrack} \\

\end{array}$$

### Examples of Use

The use of a logistic model in Colossus closely follows the other
regression options. We start by loading and preparing the data. In this
example, we will be using the Veterans’ Administration Lung Cancer study
data from the survival package. Every row has 1 trial and a binary event
status column. Once again, our model uses a trial column, an event
column, and the “risk” equation that will be converted to a probability.
In this case, the number of trials is 1, so the column can be omitted,
but to illustrate, we will artificially create a trial column named
“trial”.

``` r
#
veteran %>% setDT()
df <- copy(veteran)
# Make the same adjustments as Epicure example 6.5
karno <- df$karno
karno[93] <- 20
df$karno <- karno
df$trt <- df$trt - 1
df$trt <- as.integer(df$trt == 0)
cell_lvl <- c("large", "squamous", "smallcell", "adeno")
df$cell <- as.integer(factor(df$celltype, level = cell_lvl)) - 1
df$karno50 <- df$karno - 50

df$trial <- 1
```

We will start with testing the effects of cell type. By default,
Colossus will always remove the first level of a factored column. We can
force Colossus to include every level by defining the lowest level to an
unused value, -1 in this case. This gives the same final model as if we
removed the baseline and added an intercept column (CONST). In this
case, the p-values are for the null hypothesis that the cell type has no
effect. If we removed the baseline and added an intercept, the p-values
would be for the null hypothesis that the cell type has no effect above
the baseline.

``` r
df$cell <- factor(df$cell, levels = c(-1, 0, 1, 2, 3))

model <- Logit(trial, status) ~ loglinear(cell)
model <- Logit(status) ~ loglinear(cell)

control <- list(verbose = 0, ncores = 2)
e <- LogisticRun(model, df, control = control)
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Constant Central Estimate Standard Error
#>       <char>  <char>       <int>   <lgcl>           <char>         <char>
#> 1:    cell_0  loglin           0    FALSE             3.26          1.019
#> 2:    cell_1  loglin           0    FALSE             2.05          0.531
#> 3:    cell_2  loglin           0    FALSE             2.71          0.596
#> 4:    cell_3  loglin           0    FALSE             3.26          1.019
#>    2-tail p-value
#>            <char>
#> 1:       1.39e-03
#> 2:       1.16e-04
#> 3:       5.58e-06
#> 4:       1.39e-03
#> 
#> Logisitic Model Used
#> -2*Log-Likelihood: 64.429,  Deviation: 64.429,  AIC: 72.429,  BIC: 84.109
#> Iterations run: 11
#> maximum step size: 1.075e-04, maximum first derivative: 3.454e-05
#> Analysis converged
#> Run finished in 0.021 seconds
#> |-------------------------------------------------------------------|
```

In this analysis, we were using the default of the odds ratio. The
linking function can be changed by using the “link” option. The
following examples show the impact of using each linking function. Note
that the probability must be between 0-1, so the starting values are set
to ensure that the initial probability estimates are valid. Because we
are using categorical risk factors, every model converges to the same
score and event probabilities; however, the standard error and p-values
are affected. If we were using continuous variables, the linking
functions would give a way to define further complicated risk models.

``` r
a_n <- c(0.1, 0.1, 0.1, 0.1)
e <- LogisticRun(model, df, control = control, a_n = a_n, link = "odds")
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Constant Central Estimate Standard Error
#>       <char>  <char>       <int>   <lgcl>           <char>         <char>
#> 1:    cell_0  loglin           0    FALSE             3.26          1.019
#> 2:    cell_1  loglin           0    FALSE             2.05          0.531
#> 3:    cell_2  loglin           0    FALSE             2.71          0.596
#> 4:    cell_3  loglin           0    FALSE             3.26          1.019
#>    2-tail p-value
#>            <char>
#> 1:       1.39e-03
#> 2:       1.16e-04
#> 3:       5.58e-06
#> 4:       1.39e-03
#> 
#> Logisitic Model Used
#> -2*Log-Likelihood: 64.429,  Deviation: 64.429,  AIC: 72.429,  BIC: 84.109
#> Iterations run: 11
#> maximum step size: 9.569e-05, maximum first derivative: 3.072e-05
#> Analysis converged
#> Run finished in 0.011 seconds
#> |-------------------------------------------------------------------|

a_n <- c(-1, -1, -1, -1)
e <- LogisticRun(model, df, control = control, a_n = a_n, link = "ident")
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Constant Central Estimate Standard Error
#>       <char>  <char>       <int>   <lgcl>           <char>         <char>
#> 1:    cell_0  loglin           0    FALSE          -0.0377         0.0377
#> 2:    cell_1  loglin           0    FALSE          -0.1214         0.0607
#> 3:    cell_2  loglin           0    FALSE          -0.0645         0.0373
#> 4:    cell_3  loglin           0    FALSE          -0.0377         0.0377
#>    2-tail p-value
#>            <char>
#> 1:         0.3174
#> 2:         0.0456
#> 3:         0.0833
#> 4:         0.3174
#> 
#> Logisitic Model Used
#> -2*Log-Likelihood: 64.429,  Deviation: 64.429,  AIC: 72.429,  BIC: 84.109
#> Iterations run: 11
#> maximum step size: 4.859e-05, maximum first derivative: 1.145e-02
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.011 seconds
#> |-------------------------------------------------------------------|

a_n <- c(0.1, 0.1, 0.1, 0.1)
e <- LogisticRun(model, df, control = control, a_n = a_n, link = "loglink")
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Constant Central Estimate Standard Error
#>       <char>  <char>       <int>   <lgcl>           <char>         <char>
#> 1:    cell_0  loglin           0    FALSE            -3.28          1.000
#> 2:    cell_1  loglin           0    FALSE            -2.11          0.500
#> 3:    cell_2  loglin           0    FALSE            -2.74          0.577
#> 4:    cell_3  loglin           0    FALSE            -3.28          1.000
#>    2-tail p-value
#>            <char>
#> 1:       1.05e-03
#> 2:       2.49e-05
#> 3:       2.08e-06
#> 4:       1.05e-03
#> 
#> Logisitic Model Used
#> -2*Log-Likelihood: 64.429,  Deviation: 64.429,  AIC: 72.429,  BIC: 84.109
#> Iterations run: 12
#> maximum step size: 1.263e-04, maximum first derivative: 4.514e-05
#> Analysis converged
#> Run finished in 0.01 seconds
#> |-------------------------------------------------------------------|
```
