# Matched Case-Control Logistic Regression

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

## Matched Case-Control Modeling

One alternative regression method in Colossus is matched case-control
logistic regression. The theory is presented in the following section.

### General Theory

#### Conditional Logistic Regression

Suppose we have matched case-control data and divide our data into each
matched set. Each set has $m$ cases and $n$ records. We denote the
relative risk for individual $i$ in the set by $r_{i}$. We can calculate
the probability of case exposures conditional on all exposures in the
set by taking the ratio of the product of relative risks in the cases to
the sum of the product of relative risks for every way of selecting $m$
individuals from the $n$ at risk.

$$\begin{array}{r}
\frac{\prod\limits_{i}^{m}r_{i}}{\sum\limits_{c \in R}\left( \prod\limits_{j = 1}^{m}r_{c_{j}} \right)} \\
{L = \sum\limits_{i = 1}^{m}log\left( r_{i} \right) - log\left( \sum\limits_{c \in R}\left( \prod\limits_{j = 1}^{m}r_{c_{j}} \right) \right)}
\end{array}$$

Using the methods presented in Gail et al. (1981) we can calculate the
combination of all $n!/m!(n - m)!$ ways to select $m$ items with a more
manageable recursive formula $B(m,n)$.

$$\begin{array}{r}
{B(m,n) = \sum\limits_{c \in R}\left( \prod\limits_{j = 1}^{m}r_{c_{j}} \right)} \\
{B(m,n) = B(m,n - 1) + r_{n}B(m - 1,n - 1)} \\
{B(m,n) = \begin{cases}
{\sum\limits_{j}^{n}r_{j}} & {m = 1} \\
0 & {m > n}
\end{cases}}
\end{array}$$

We can then directly solve for the first and second derivatives and
their recursive formula.

$$\begin{array}{r}
{\frac{\partial r_{i}}{\partial\beta_{\mu}} = :r_{i}^{\mu}} \\
{\frac{\partial B(m,n)}{\partial\beta_{\mu}} = B^{\mu}(m,n) = \sum\limits_{c \in R}\left\lbrack \left( \sum\limits_{j = 1}^{m}\frac{r_{c_{j}}^{\mu}}{r_{c_{j}}} \right)\prod\limits_{j = 1}^{m}r_{c_{j}} \right\rbrack} \\
{B^{\mu}(m,n) = B^{\mu}(m,n - 1) + r_{n}B^{\mu}(m - 1,n - 1) + r_{n}^{\mu}B(m - 1,n - 1)} \\
{B^{\mu}(m,n) = \begin{cases}
{\sum\limits_{j}^{n}r_{j}^{\mu}} & {m = 1} \\
0 & {m > n}
\end{cases}}
\end{array}$$

$$\begin{array}{r}
{\frac{\partial^{2}r_{i}}{\partial\beta_{\mu}\partial\beta_{\nu}} = :r_{i}^{\mu,\nu}} \\
{\frac{\partial^{2}B(m,n)}{\partial\beta_{\mu}\partial\beta_{\nu}} = B^{\mu,\nu}(m,n) = \sum\limits_{c \in R}\left\lbrack \left( \sum\limits_{j = 1}^{m}\frac{r_{c_{j}}^{\mu,\nu}}{r_{c_{j}}} + \left( \sum\limits_{j = 1}^{m}\frac{r_{c_{j}}^{\mu}}{r_{c_{j}}} \right)\left( \sum\limits_{j = 1}^{m}\frac{r_{c_{j}}^{\nu}}{r_{c_{j}}} \right) - \sum\limits_{j = 1}^{m}\frac{r_{c_{j}}^{\mu}}{r_{c_{j}}}\frac{r_{c_{j}}^{\nu}}{r_{c_{j}}} \right)\prod\limits_{j = 1}^{m}r_{c_{j}} \right\rbrack} \\
{B^{\mu,\nu}(m,n) = B^{\mu,\nu}(m,n - 1) + r_{n}^{\mu,\nu}B(m - 1,n - 1) + r_{n}^{\nu}B^{\mu}(m - 1,n - 1) + r_{n}^{\mu}B^{\nu}(m - 1,n - 1) + r_{n}B^{\mu,\nu}(m - 1,n - 1)} \\
{B^{\mu,\nu}(m,n) = \begin{cases}
{\sum\limits_{j}^{n}r_{j}^{\mu,\nu}} & {m = 1} \\
0 & {m > n}
\end{cases}}
\end{array}$$

Finally, these expressions for $B(m,n)$ can be substituted into the
equations for the contribution of Log-Likelihood and its derivatives
from each matched set. The model is then optimized via the same methods
as the other regression models.

$$\begin{array}{r}
{L_{set} = \sum\limits_{i = 1}^{m}log\left( r_{i} \right) - log\left( B(m,n) \right)} \\
{L_{set}^{\mu} = \sum\limits_{i = 1}^{m}\frac{r_{i}^{\mu}}{r_{i}} - \frac{B^{\mu}(m,n)}{B(m,n)}} \\
{L_{set}^{\mu,\nu} = \sum\limits_{i = 1}^{m}\left( \frac{r_{i}^{\mu,\nu}}{r_{i}} - \frac{r_{i}^{\mu}}{r_{i}}\frac{r_{i}^{nu}}{r_{i}} \right) - \left( \frac{B^{\mu,\nu}(m,n)}{B(m,n)} - \frac{B^{\mu}(m,n)}{B(m,n)}\frac{B^{\nu}(m,n)}{B(m,n)} \right)}
\end{array}$$

#### Unconditional Logistic Regression

It is important to note that the recursive formula calculation can
quickly become time-consuming, particularly if there is a large number
of cases. To make the matched case-control method generally applicable,
the likelihood function can be changed to a logistic regression model in
matched sets with a large number of cases. In general, the matched
case-control regression function adds an item to the model control list,
“cond_thres”, to set the threshold to switch to a logistic regression
model.

The logistic loglikelihood is defined by treating the matched
case-control data as single trial data. The likelihood is a function of
the event status ($\theta_{i}$) and odds ratio ($O_{i}$). The odds ratio
for any row is calculated as the product of the odds ratio for the
matched set ($O_{s}$) and the relative risk for the row ($r_{i}$).
Behind the scenes, Colossus optimizes both the model parameters
($\overset{\rightarrow}{\beta}$) for the relative risk as well as a
logistic model for the matched set odds ratios
($O_{s} = e^{\alpha_{s}}$).

$$\begin{array}{r}
{L_{i} = y_{i}ln\left( O_{s}r_{i} \right) - ln\left( 1 + O_{s}r_{i} \right)} \\
{\frac{\partial L_{i}}{\partial\beta_{j}} = y_{i}\frac{r_{i}^{j}}{r_{i}} - O_{s}\frac{r_{i}^{j}}{1 + O_{s}r_{i}}} \\
{\frac{\partial L_{i}}{\partial\alpha_{s}} = \left( y_{i} - 1 \right) + \frac{1}{1 + O_{s}r_{i}}} \\
{\frac{\partial^{2}L_{i}}{\partial\beta_{j}\partial\beta_{k}} = y_{i}\left( \frac{r_{i}^{j,k}}{r_{i}} - \frac{r_{i}^{j}}{r_{i}}\frac{r_{i}^{k}}{r_{i}} \right) - O_{s}\left( \frac{r_{i}^{j,k}}{1 + O_{s}r_{i}} - O_{s}\frac{r_{i}^{j}}{1 + O_{s}r_{i}}\frac{r_{i}^{k}}{1 + O_{s}r_{i}} \right)} \\
{\frac{\partial^{2}L_{i}}{\partial\beta_{j}\partial\alpha_{s}} = \frac{- O_{s}r_{i}^{j}}{1 + O_{s}r_{i}}} \\
{\frac{\partial^{2}L_{i}}{\partial\alpha_{s}^{2}} = \frac{- O_{s}r_{i}}{1 + O_{s}r_{i}}}
\end{array}$$

### Examples of Use

The following provides a basic example of how to use the matched
case-control regression function. The data used is from a lung cancer
study, included in the R survival library. Slight adjustments were made
to make the data line up with data included with 32-bit Epicure, for
comparison. In short, we want to model the effects of treatment and
Karnofsky performance score on sets matched by cancer cell type.

``` r
#
data(cancer, package = "survival")
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
```

In the following examples, we are using matching by strata and changing
the conditional threshold. In the first case, every matched set uses the
recursive formula. In the second case, one set uses the simplified
formula. Finally in the third case, every set uses the simplified
formula. In all cases, the regression returns the typical output and can
be summarized similarly to other regression function outputs.

``` r
model <- CaseControl_Strata(status, cell) ~ loglinear(karno50, trt)


control <- list(verbose = 2, maxiters = c(25, 25), ncores = 2)
e0 <- CaseControlRun(model, df, control = control, conditional_threshold = 100)
e1 <- CaseControlRun(model, df, control = control, conditional_threshold = 40)
e2 <- CaseControlRun(model, df, control = control, conditional_threshold = 0)


print(e0)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:   karno50  loglin           0             0.01          0.017          0.556
#> 2:       trt  loglin           0             0.01          0.700          0.989
#> 
#> Matched Case-Control Model Used
#> Model stratified by 'cell'
#> Deviance: 57.029
#> 0 out of 4 matched sets used Unconditional Likelihood
#> Iterations run: 3
#> maximum step size: 3.052e-05, maximum first derivative: 1.465e+02
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.022 seconds
#> |-------------------------------------------------------------------|
print(e1)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:   karno50  loglin           0             0.01         0.0171          0.558
#> 2:       trt  loglin           0             0.01         0.5726          0.986
#> 
#> Matched Case-Control Model Used
#> Model stratified by 'cell'
#> Deviance: 59.956
#> 1 out of 4 matched sets used Unconditional Likelihood
#> Iterations run: 3
#> maximum step size: 3.052e-05, maximum first derivative: 1.463e+02
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.011 seconds
#> |-------------------------------------------------------------------|
print(e2)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:   karno50  loglin           0             0.01         0.0172          0.561
#> 2:       trt  loglin           0             0.01         0.5164          0.985
#> 
#> Matched Case-Control Model Used
#> Model stratified by 'cell'
#> Deviance: 67.073
#> 4 out of 4 matched sets used Unconditional Likelihood
#> Iterations run: 3
#> maximum step size: 3.052e-05, maximum first derivative: 1.498e+02
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.011 seconds
#> |-------------------------------------------------------------------|
```
