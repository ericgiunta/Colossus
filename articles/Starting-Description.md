# Colossus Description

``` r
Sys.setenv("OMP_THREAD_LIMIT" = 1) # Reducing core use, to avoid accidental use of too many cores
library(Colossus)
library(data.table)
library(parallel)
```

## Colossus Description

## Model Structure

At its full potential, Colossus can analyze a vast number of possible
risk/rate models. Specifically, Colossus is designed to allow for a
combination of linear and non-linear models to estimate Cox proportional
hazard ratios, Poisson model rates, and Fine-Gray Competing risk ratios.
The simplest of which is an exponential model. The simplest model
generally used for hazard ratios in survival analysis is the exponential
of a linear function of covariates.

$$\begin{array}{r}
{R\left( \overset{\rightarrow}{\beta},\overset{\rightarrow}{x} \right) = \exp\left( \overset{\rightarrow}{\beta} \cdot \overset{\rightarrow}{x} \right)}
\end{array}$$

In Colossus, this general model is extended by abstracting to the sum
and product of terms and subterms. In this vignette, the value
calculated to quantify the risk of an event will be denoted as the risk,
R, however the exact meaning differs. In Cox proportional hazard
modeling, the risk calculated is a hazard ratio. In Poisson modeling the
risk calculated is an estimated number of events per person-year. The
risk ($R$) has a set formula dependent on terms ($T$). Each term has a
formula dependent on the product of subterms ($S$). Each subterm is a
function of a covariates ($x$) and parameters ($\alpha,\beta$).

There are currently five types of risk models available. The risk can be
expressed as an additive model ($R_{A}$), product additive model
($R_{PA}$), Product additive excess model ($R_{PAE}$), multiplicative
relative model ($R_{M}$), multiplicative excess model ($R_{ME}$), or the
geometric mixture model with relative risks or excess risks
($R_{GMIX}$).

$$\begin{aligned}
R_{A} & {= \sum\limits_{i = 0}^{n}T_{i}} \\
R_{PA} & {= T_{0} \times \sum\limits_{i = 1}^{n}T_{i}} \\
R_{PAE} & {= T_{0} \times \left( 1 + \sum\limits_{i = 1}^{n}T_{i} \right)} \\
R_{M} & {= T_{0} \times \prod\limits_{i = 1}^{n}\left( T_{i} \right)} \\
R_{ME} & {= T_{0} \times \prod\limits_{i = 1}^{n}\left( 1 + T_{i} \right)} \\
R_{GMIX} & {= T_{0} \times \left( \prod\limits_{i = 1}^{n}\left( T_{i}^{*} \right) \right)^{\theta} \times \left( 1 + \sum\limits_{i = 1}^{n}\left( T_{i}^{*} - 1 \right) \right)^{1 - \theta}} \\
T_{i}^{*} & {= \begin{cases}
T_{i} & \text{Relative Risk} \\
{T_{i} + 1} & \text{Excess Risk}
\end{cases}} \\
 & 
\end{aligned}$$

Each term is composed of a combination of 4 types of subterms. Every
covariate is part of a log-linear subterm, a linear subterm, a
product-linear subterm, or a general non-linear term. The log-linear
subterm ($S_{LL}$) is the exponential of a linear combination of
covariates. The linear subterm ($S_{L}$) is a linear combination of
covariates. The product-linear subterm ($S_{PL}$) is one plus a linear
combination of covariates. The general non-linear term ($S_{NL}$) is the
sum of exponential terms, quadratic terms, linear-threshold terms
($F_{LT}$), step function terms ($F_{STP}$), linear-quadratic terms
($F_{LQ}$), and linear-exponential terms ($F_{LEXP}$). Each term is the
product of the non-empty subterms.

$$\begin{aligned}
S_{LL} & {= \prod\limits_{i}\left( \exp\left( x_{i} \cdot \beta_{i} \right) \right)} \\
S_{L} & {= \sum\limits_{i}\left( x_{i} \cdot \beta_{i} \right)} \\
S_{PL} & {= 1 + \sum\left( x_{i} \cdot \beta_{i} \right)} \\
S_{NL} & {= \sum\limits_{i}\left( \alpha_{i} \times \exp\left( x_{i} \cdot \beta_{i} \right) \right) + \sum\limits_{i}\left( \beta_{i} \cdot \left( x_{i} \right)^{2} \right)} \\
 & {+ \sum\limits_{i}F_{LT} + \sum\limits_{i}F_{STP} + \sum\limits_{i}F_{LQ} + \sum\limits_{i}F_{LEXP}} \\
F_{LT} & {= \begin{cases}
{\alpha_{i} \cdot \left( x - \beta_{i} \right)} & \left( x > \beta_{i} \right) \\
0 & \text{else}
\end{cases}} \\
F_{STP} & {= \begin{cases}
\alpha_{i} & \left( x > \beta_{i} \right) \\
0 & \text{else}
\end{cases}} \\
F_{LQ} & {= \begin{cases}
{\beta_{i} \cdot x} & \left( x > \alpha_{i} \right) \\
{\lambda_{i} \cdot x^{2} + \nu_{i}} & \text{else}
\end{cases}} \\
F_{LEXP} & {= \begin{cases}
{\beta_{i} \cdot x} & \left( x > \alpha_{i} \right) \\
{\lambda_{i} - \exp\left( \nu_{i} + \mu \cdot x \right)} & \text{else}
\end{cases}} \\
T_{j} & {= S_{LL,j} \times S_{L,j} \times S_{PL,j} \times S_{NL,j}}
\end{aligned}$$

In short, every element of the risk model has a subterm type, term
number, covariate, and parameter value.

## Using The Standard Model

In Colossus, every equation is defined in a very similar fashion.
Elements of the risk equation can be viewed as rows of a table: the
columns store covariate names, term numbers, subterms types, and a
starting point. For more complex regression most parameters can also be
set to remain constant over the regression. Doing so forces the
parameter value to be constant but not remove the element from the risk
calculations. Constant elements are used for risk calculation but are
not used for calculating steps or standard deviations. For a
multiplicative excess risk model, the following table and equation are
equivalent.

| Term number | subterm type | covariate |
|:-----------:|:------------:|:---------:|
|      0      |  LogLinear   |     a     |
|      1      |    Linear    |     b     |
|      1      |    Linear    |     c     |

$$\begin{array}{r}
{R = \exp\left( \beta_{a} \cdot x_{a} \right) \times \left( 1 + \beta_{b} \cdot x_{b} + \beta_{c} \cdot x_{c} \right)}
\end{array}$$

If the user wanted to update the model to include a new term with a
product-linear subterm, then the table would only need to be updated
with a new row.

| Term number |  subterm type  | covariate |
|:-----------:|:--------------:|:---------:|
|      0      |   LogLinear    |     a     |
|      1      |     Linear     |     b     |
|      1      |     Linear     |     c     |
|      2      | Product-Linear |     d     |

$$\begin{array}{r}
{R = \exp\left( \beta_{a} \cdot x_{a} \right) \times \left( 1 + \left( \beta_{b} \cdot x_{b} + \beta_{c} \cdot x_{c} \right) \right) \times \left( 2 + \beta_{d} \cdot x_{d} \right)}
\end{array}$$

Note that the multiplicative excess model is taking the product of
terms, each of which is the product of subterms. The same equation may
be written with subterms distributed over different terms in the
multiplicative model without changing the final value. The only
exception is subterms moved from or into the “first” term in the
multiplicative excess model. The same equation could be expressed in
ways with different computational complexity.

Models in Colossus are passed in the form of a formula. Each subterm of
the model is expressed as an element of the formula. The model
identifiers (A, M, PA, PAE, GMIX) are included at the end of the
formula. Every model except the additive and multiplicative models have
a distinction between the first and remaining terms. Colossus defaults
to term 0 being the unique term.

The final table is also equivalent to the following code.

``` r
names <- c("a", "b", "c", "d")
term_n <- c(0, 1, 1, 2)
tform <- c("loglin", "lin", "lin", "plin")
modelform <- "M"

a_n <- c(0.1, 0.1, 0.1, 0.1)

model <- outcome ~ loglinear(a, 0) + linear(b, c, 1) +
  plinear(d, 2) + ME()
```

Interior colossus functions convert between the vectors and formula
representations. *names* denotes the column names used, *term_n* denotes
the term numbers, *tform* denotes the subterm formula, *modelform*
denotes the term formula, and *a_n* denotes the initial guesses for the
parameter values. All of which are assumed to be in the same order. A
function is called before regression that reorders the inputs in order
of terms, subterm types, etc.

## Survival Time and Event Data

Colossus performs survival analysis via either a Cox Proportional
Hazards or a Poisson model regression. In both cases, the user specifies
which columns contain time duration and events of interest. For the
Poisson model regression, these would be which column contains the
person-years and number of events for each row of data. For the Cox
Proportional Hazards regression, the user identifies which columns
provide starting and ending times, and what column gives the event
status. The event status is assumed to be a binary covariate which is 1
for intervals containing an event. Colossus supports left-censored data,
right-censored data, and interval-censored data. The data is interval
censored by default, so left and right censoring are handled by defining
interval endpoints outside the minimum or maximum event times. For
Poisson model regression the user only provides columns for the
person-years per row and number of events. Poisson model regression
supports any non-negative number of events.

On the user side, the names of columns containing the time and events
need to be given. Let us assume we have a dataframe organized as
follows:

| UserID | Starting_Age | Ending_Age | Cancer_Status |
|:------:|:------------:|:----------:|:-------------:|
|  112   |      18      |     30     |       0       |
|  114   |      20      |     45     |       0       |
|  213   |      18      |     57     |       1       |
|   …    |      …       |     …      |       …       |

For Cox proportional hazard, we need to provide three column names: the
starting age, the ending age, and if an event happened during the
interval. In this case that would be “Starting_Age”, “Ending_Age”, and
“Cancer_Status”. In this example, we have interval data but that may not
always be the case. Colossus is designed to allow the user to omit
interval columns to represent truncation. Colossus assumes that this
means that the person has the missing endpoint outside of the available
data range and creates a dummy column to reference as the interval
endpoint. In the formula, the interval endpoints can be named “tstart”
and “tend”. So in code, these variables would look like this:

``` r
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
# For the interval case
tstart <- "Starting_Age"
tend <- "Ending_Age"
event <- "Cancer_Status"
```

Suppose instead we were interested in a Poisson model regression, the
only difference is instead of providing interval endpoints we give a
column of interval lengths. Suppose we are now working with the
following table:

| UserID | Starting_Age | Ending_Age | Person_Years | Cancer_Status |
|:------:|:------------:|:----------:|:------------:|:-------------:|
|  112   |      18      |     30     |      12      |       0       |
|  114   |      20      |     45     |      25      |       0       |
|  213   |      18      |     57     |      39      |       1       |
|   …    |      …       |     …      |      …       |       …       |

We would just provide the duration column and event column, or
“Person_Years” and “Cancer_Status”. Which in code looks like the
following:

``` r
df$Person_Years <- df$Ending_Age - df$Starting_Age
pyr <- "Person_Years"
event <- "Cancer_Status"
```

The regression types are added to the formula on the left-hand side.
Similar to the Surv function in the survival library, the left-hand side
of the formula specifies the type of model and which time/event columns
are used.

``` r
# For the interval case
RHS <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ risk_factors

# Supposing we had left truncated data the following would change
RHS <- Cox(tstart = Starting_Age, event = Cancer_Status) ~ risk_factors

# and with right truncated data the following is used
RHS <- Cox(Ending_Age, Cancer_Status) ~ risk_factors
```

## Control and Verbosity

Colossus offers several options that control how the input data is used
and returned. This can be split into two general categories. The first
category is convergence parameters. These cover the number of iterations
used, limits on parameter changes, and stopping criteria. The second
category is parameters used for debugging and additional information.

The common convergence parameters are as follows:

|   Parameter   |                                                                                                                                                                                                                                                               Description                                                                                                                                                                                                                                                               |
|:-------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|      lr       |                                                                                                                                                                      learning rate applied to steps taken, smaller values avoid overshooting at the risk of slower convergence. Applied to the solution for setting the first derivative of Log-Likelihood to zero                                                                                                                                                                      |
|    maxiter    |                                                                                                                                                                                                                                                    maximum number of iterations run                                                                                                                                                                                                                                                     |
|    halfmax    |                                                                                                                                                                                                                                            Maximum number of half-steps taken each iteration                                                                                                                                                                                                                                            |
|    epsilon    |                                                                                                                                                                                                                      Minimum parameter change above convergence. If the change is below this, the regression exits                                                                                                                                                                                                                      |
|   dbeta_max   |                                                                                                                                                                                                                      learning rate applied to limit on the steps taken. Applied to solution to zero Log-Likelihood                                                                                                                                                                                                                      |
| deriv_epsilon |                                                                                                                                                                                                     Minimum Log-Likelihood first derivative above convergence. If the largest first derivative is below this, the regression exits                                                                                                                                                                                                      |
|   step_max    |                                                                                                                                                                                          Largest starting step size. Code reduces it every time halfmax is hit without finding an improvement. When step_max is below epsilon, the code exits                                                                                                                                                                                           |
|    verbose    | integer valued 0-4 controlling what information is printed to the terminal. Each level includes the lower levels. 0: silent, 1: errors printed, 2: warnings printed, 3: notes printed, 4: debug information printed. Errors are situations that stop the regression, warnings are situations that assume default values that the user might not have intended, notes provide information on regression progress, and debug prints out C++ progress and intermediate results. The default level is 2 and True/False is converted to 3/0. |
|    ncores     |                                                                                                                                                                                                                         Number of physical cores assigned, cannot be greater than the number of cores available                                                                                                                                                                                                                         |

The previous options are all contained within a list of parameters,
however, there are additional standalone parameters that control the
convergence. A common example is *keep_constant* which can be used to
hold parameters constant. *keep_constant* can be assigned a vector of
0/1 for each element to denote free (0) and constant (1) parameters.

Going back to our example, the code may look like the following:

``` r
keep_constant <- c(0, 0, 0, 0)

control <- list(
  "ncores" = 1, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
  "dbeta_max" = 0.5, "deriv_epsilon" = 1e-3, "step_max" = 1.0,
  "verbose" = 2
)
```

In this example, no parameters are held constant and 100 iterations with
5 half-steps are used. There are no errors for unused control parameters
being used, they will simply be ignored. Every commonly required control
parameter has a default value which is automatically used if not
declared.

## Running the Regression and Output

Finally, the user calls the regression function they need. The names
used throughout this vignette are the defaults assumed. Both the Cox PH
regression and Poisson model regression are called directly and return a
list of results. Colossus contains a suite of additional checks it runs
on the inputs before starting any calculations, which are designed to
output explicit details on any issues. Printing error details may
require the verbose option to be TRUE. In code, the functions are called
as follows:

``` r
# assuming the table of covariates is stored in a data.table "df"

model_cox <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, 0) +
  linear(b, c, 1) + plinear(d, 2) + multiplicative()

e <- CoxRun(model_cox, df, a_n = a_n, control = control)
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:         a  loglin           0             31.1        1079824          1.000
#> 2:         b     lin           1             14.4            NaN            NaN
#> 3:         c     lin           1             14.4            NaN            NaN
#> 4:         d    plin           2             31.1            477          0.948
#> 
#> Cox Model Used
#> -2*Log-Likelihood: 1.374,  AIC: 9.374
#> Iterations run: 30
#> maximum step size: 1.000e+00, maximum first derivative: 9.764e-04
#> Analysis converged
#> Run finished in 0.024 seconds
#> |-------------------------------------------------------------------|

# or a Poisson model regression
model_pois <- Poisson(Person_Years, Cancer_Status)  ~ loglinear(a, 0) + linear(b, c, 1) + plinear(d, 2) + multiplicative()
e <- PoisRun(model_pois, df, a_n = a_n, control = control)
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Term Number Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>       <int>            <num>          <num>          <num>
#> 1:         a  loglin           0         -0.29681        1.10226       7.88e-01
#> 2:         b     lin           1         -0.01166        0.04784       8.07e-01
#> 3:         c     lin           1          0.00725        0.00993       4.65e-01
#> 4:         d    plin           2         -0.84726        0.19126       9.42e-06
#> 
#> Poisson Model Used
#> -2*Log-Likelihood: 11.992,  Deviation: 7.992,  AIC: 15.992,  BIC: 19.776
#> Iterations run: 5
#> maximum step size: 9.766e-04, maximum first derivative: 4.346e+02
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.012 seconds
#> |-------------------------------------------------------------------|
```

The following are output in the regression output lists:

|        Item        |                                                           Description                                                           |
|:------------------:|:-------------------------------------------------------------------------------------------------------------------------------:|
|       LogLik       |                                                      Final Log-Likelihood                                                       |
|     First_Der      |                                      First derivative of Log-Likelihood at the final point                                      |
|     Second_Der     |                                Matrix of second derivatives of Log-Likelihood at the final point                                |
|       beta_0       |                                                    The final parameter list                                                     |
| Standard_Deviation |                                       Inverse-Hessian estimate of the standard deviation                                        |
|        AIC         |                                       Akaike information criterion for the final solution                                       |
|     Deviation      |   Deviance calculated using the difference in log-likelihood contributions between the saturated model and the current model    |
|  Parameter_Lists   |                            the term numbers and term formulas used, to compare against the expected                             |
|    Control_List    | a list containing the number of iterations, the maximum step at the final iteration, and the highest first derivative magnitude |
|     Converged      |                                             TRUE/FALSE if the regression converged                                              |

In this case, the poisson regression did not converge. Given that the
examples were arbitrary this is not unexpected. The Poisson model
regression exited due to the step limit being below the threshold, which
means that the score was still changing significantly even at low step
sizes. To find a better solution the user would likely have to either
change the model equation or provide a new starting guess with a lower
step size limit. This is an example of how to interpret the results of a
Colossus run.

There are additional variants of these functions which are described in
greater detail in the other vignettes.
