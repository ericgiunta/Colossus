# Alternative Regression Options

``` r
Sys.setenv("OMP_THREAD_LIMIT" = 1) # Reducing core use, to avoid accidental use of too many cores
library(Colossus)
library(data.table)
```

## General Options

Both Cox proportional hazards and Poisson model regressions have
additional functions that can account for specific situations. There
general situations are as follows:

|               Option               |                                                                                                          Description                                                                                                          | Cox PH | Poisson |
|:----------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:------:|:-------:|
|           Stratification           | In Cox PH, the stratification is applied in the risks groups to compare only like rows. In Poisson model regression the stratification is an additional term in the model to account for the effects of stratified covariates |   x    |    x    |
|          Simplified Model          |                                                         For a multiplicative Log-Linear model, simplifications can be made to the code. These give a faster version.                                                          |   x    |         |
|     Non-Derivative Calculation     |                                                                If a single iteration is needed without derivatives, there is time and memory that can be saved                                                                |   x    |    x    |
|          Competing Risks           |                                    If there is a competing event, rows are weighted by an estimate of the censoring rate to approximate a dataset without the effect of a competing event                                     |   x    |         |
|           Joint Analysis           |              A Poisson method which allows a model to solve for multiple outcomes, by creating multiple copies of rows with multiple events and using a factor variable to add covariates for specific outcomes               |        |    x    |
| Distributed Initial Parameter Sets |                                  The user provides distribution for each starting parameter, Colossus runs low iteration regressions at random points, and then a full run at the best guess                                  |   x    |    x    |

The following sections review the math behind the basic functions and
how each option changes that.

## Stratification

### Cox Proportional Hazards

In Cox Proportional Hazards, the Log-Likelihood is calculated by taking
the ratio of the hazard ratio at each event to the sum of the hazard
ratios of every other row at risk. This defines the risk group for every
event time to be the intervals containing the event time. Intervals are
assumed to be open on the left and closed on the right, and events are
assumed to take place at the right endpoint. This gives the following
common equation for the Log-Likelihood:

$$\begin{array}{r}
{Ll = \prod\limits_{i}^{n}\left( \frac{r_{i}}{\sum\limits_{j:t_{j} \in R_{i}}r_{j}} \right)^{\delta_{i}}}
\end{array}$$

In which r denotes hazard ratios, the denominator is the sum of hazard
ratios of intervals containing the event time, and each term is raised
to the power of 1 if the interval has an event and 0 otherwise.
Different tie methods modify the denominator based on how the order of
events is assumed, but the general form still stands. The goal is to
compare each event interval to intervals within a similar time span.
Stratification adds a condition: if the user stratifies over a covariate
“F” then each risk group is split into subgroups with the same value of
“F”. So the goal becomes to compare event intervals to intervals with
similar strata AND time. This is done to remove the influence of the
stratification variables from the calculations.

In code, this is done by adding an additional parameter for the
stratification column and using a different response term in the model.
Multiple strata columns can be provided by including a vector of
columns.

``` r
Strat_Col <- "s0"
e <- CoxRun(Cox_Strata(time1, time2, event, s0) ~ loglinear(dose), df,
  a_n = a_n, control = control
)
Strat_Cols <- c("s0", "s1", "s2")
e <- CoxRun(Cox_Strata(time1, time2, event, c(s0, s1, s2)) ~ loglinear(dose), df,
  a_n = a_n, control = control
)
```

### Poisson Regression

Poisson model regression does not have risk groups to account for, but
the theory is the same. To remove the influence of stratification
covariate a new term is added to account for their effects. In Colossus,
this is a Log-Linear term. So the following may be the model used
without stratification:

$$\begin{array}{r}
{R = \sum\limits_{i}\left( x_{i} \cdot \beta_{i} \right)}
\end{array}$$

Then the stratified model may look like this:

$$\begin{array}{r}
{R = \left( \sum\limits_{i}\left( x_{i} \cdot \beta_{i} \right) \right) \times \left( \exp\left( \sum\limits_{strata}\left( x_{strata} \cdot \beta_{strata} \right) \right) \right)}
\end{array}$$

Only results associated with the non-stratified parameters are returned
by default. The strata effect is calculated by taking the ratio of the
events in each strata to the sum of person-years multiplied by risk in
each strata. The derivatives are then written, substituting the strata
effect.

$$\begin{array}{r}
{\beta_{strata} = \left( \sum\limits_{i}{event}_{i} \right)/\left( \sum\limits_{i}P_{i} \cdot R_{i} \right)}
\end{array}$$

In code, this is done by adding in a list of stratification columns and
using a different response term in the model. Colossus combines the list
of stratification columns into a single interaction, so it does not
matter if the user provides a list or combines them themselves.

``` r
Strat_Col <- c("e")
e <- PoisRun(Poisson_Strata(pyr, event, e) ~ loglinear(dose), df,
  a_n = a_n, control = control
)
```

## Non-Derivative Calculation

Colossus uses Newton’s method to perform the regression, which can
become computationally intensive as the model becomes more complicated.
So, Colossus contains functions to calculate only the scores for a
parameter set and skip the derivative calculations. These results could
be used to perform a bisection method of regression or plot the
dependence between the score and parameter values. This ability is left
for the user’s convenience.

The code is similar to previous examples:

``` r
e <- CoxRun(Cox(time1, time2, event) ~ loglinear(dose), df,
  a_n = a_n, control = control, single = TRUE
)
```

## Competing Risks

### Fine-Gray

In Cox PH there is an assumption that every individual is recorded until
they have an event or they are naturally censored. Censoring is assumed
to be statistically independent with regard to the event being studied.
These assumptions are commonly violated when there is a competing event
occurring. In sensitivity analysis, two extremes are generally tested.
Either every person with a competing event is treated as having the
event of interest instead, or they are assumed to never experience the
event of interest. However, there are methods to find a more realistic
alternative. Colossus applies the Fine-Gray model for competing risks,
which instead weights the contribution of competing event intervals in
future intervals by the probability they would not have been censored.

As previously established, the risk groups are formed to measure the
probability that an individual survived up to the time and experienced
the event of interest, given that they did not experience an event up to
that time:

$$\begin{array}{r}
{\lambda(t) = \lim\limits_{\Delta t\rightarrow 0}\frac{\left( P(t \leq T \leq t + \Delta t){\mspace{6mu}\text{and}\mspace{6mu}}(k = 1),T \geq t \right)}{\Delta t}}
\end{array}$$

The competing risks model adjusts this to be the probability that an
individual survived up to the time and experienced the event of
interest, given that they did not experience an event up to that time OR
they survived up to a previous time and experienced a different event:

$$\begin{array}{r}
{\lambda(t) = \lim\limits_{\Delta t\rightarrow 0}\frac{\left( P(t \leq T \leq t + \Delta t){\mspace{6mu}\text{and}\mspace{6mu}}(k = 1),(T \geq t){\mspace{6mu}\text{or}\mspace{6mu}}\left( (T < t){\mspace{6mu}\text{and}\mspace{6mu}}(k \neq 1) \right) \right)}{\Delta t}}
\end{array}$$

This means that risk groups would contain both intervals actually at
risk and intervals with competing events treated as if they were at
risk. If we assume that no time-dependent covariates are being used,
then all that remains is to weigh the contribution of these competing
events. In general, a weighting column is provided before regression,
the values can be solved from a survival curve or by using the
“finegray” function of the survival library.

In code, the call is similar to the standard function. The process for
finding a weighting is included below. In this example, we assume the
event column, lung, to contain a 0 for no events, 1 for the primary
event, and 2 for the competing event.

``` r
pdata <- finegray(Surv(time2, event) ~ ., data = df)

e <- CoxRun(FineGray(fgstart, fgstop, fgstatus, fgwt) ~ loglinear(dose), pdata,
  a_n = a_n, control = control
)
```

### Poisson Joint Analysis

For a dataset with multiple outcomes, there are often two ways that
Poisson models are fit. Either multiple independent models are fit or
the events are combined and one model is fit to several events. These
methods limit models to be completely independent or identical. The true
model may be a combination of shared terms and terms specific to each
event, which cannot be modeled with these methods. To fit this type of
model a joint analysis method (Cologne, 2019) is available in Colossus.

Suppose one has a table of person-years, a covariate, and counts for two
events. Assume we are fitting the event rate ($\lambda$) for each event
($y,z$) and have reason to believe that the background rate ($\beta$) is
the same for each event.

|  Time   |    a    |    y    |    z    |
|:-------:|:-------:|:-------:|:-------:|
| $t_{1}$ | $a_{1}$ | $y_{1}$ | $z_{1}$ |
| $t_{2}$ | $a_{2}$ | $y_{2}$ | $z_{2}$ |

$$\begin{array}{r}
{\lambda_{y}(a) = \beta*\exp\left( \mu_{y}*a \right)} \\
{\lambda_{z}(a) = \beta*\exp\left( \mu_{z}*a \right)}
\end{array}$$

If one were to solve the equations separately the table would be split
into two tables, each with one event column. The premise of a joint
analysis is to write a model that can be applied to every event,
including a factor covariate to select which event is being solved for.
Then the split tables can be recombined and multiple events can be
solved at once allowing for shared and event-specific parameters.

|  Time   |    a    | $\alpha_{y}$ | $\alpha_{z}$ | events  |
|:-------:|:-------:|:------------:|:------------:|:-------:|
| $t_{1}$ | $a_{1}$ |      1       |      0       | $y_{1}$ |
| $t_{1}$ | $a_{1}$ |      0       |      1       | $z_{1}$ |
| $t_{2}$ | $a_{2}$ |      1       |      0       | $y_{2}$ |
| $t_{2}$ | $a_{2}$ |      0       |      1       | $z_{2}$ |

$$\begin{array}{r}
{\lambda\left( a,\alpha_{y},\alpha_{z} \right) = \beta*\exp\left( \mu_{y}*\left( a*\alpha_{y} \right) + \mu_{z}*\left( a*\alpha_{z} \right) \right)}
\end{array}$$

Colossus includes several functions to apply this method, a function
that produces input for a joint analysis regression and a function that
also runs the regression. The most general converts a table and lists of
formula into the input for a regression.

``` r
a <- c(0, 0, 0, 1, 1, 1)
b <- c(1, 1, 1, 2, 2, 2)
c <- c(0, 1, 2, 2, 1, 0)
d <- c(1, 1, 0, 0, 1, 1)
e <- c(0, 1, 1, 1, 0, 0)
df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
time1 <- "t0"
time2 <- "t1"
df$pyr <- df$t1 - df$t0
pyr <- "pyr"
events <- c("e0", "e1")
```

Colossus generally accepts a formula to describe the elements of the
model. For a joint analysis, Colossus instead expects a list of formula.
The formula are expect to cover the model elements specific to each
event, and potentially a model for shared elements. Colossus expects the
name of the shared model to be “shared”. The left hand side of the
shared model is not used.

``` r
model_1 <- Pois(pyr, e0) ~ loglin(fac, 0)
model_2 <- Pois(pyr, e1) ~ loglin(fac, 0)
model_s <- Pois(pyr) ~ plinear(t0, 0)
formula_list <- list(model_1, model_2, "shared" = model_s)
```

The function returns a list containing the combined table and model
object.

``` r
get_form_joint(formula_list, df, nthreads = 1)
#> $model
#> $person_year
#> [1] "pyr"
#> 
#> $event
#> [1] "events"
#> 
#> $strata
#> [1] "NONE"
#> 
#> $null
#> [1] FALSE
#> 
#> $term_n
#> [1] 0 0 0
#> 
#> $tform
#> [1] "plin"   "loglin" "loglin"
#> 
#> $names
#> [1] "t0"     "fac_e0" "fac_e1"
#> 
#> $a_n
#> [1] 0 0 0
#> 
#> $keep_constant
#> [1] 0 0 0
#> 
#> $modelform
#> [1] "ME"
#> 
#> $gmix_term
#> [1] 1
#> 
#> $gmix_theta
#> [1] 0
#> 
#> $expres_calls
#> list()
#> 
#> attr(,"class")
#> [1] "poismodel"
#> 
#> $data
#>        t0    t1 events    e0    e1   fac   pyr fac_e0 fac_e1
#>     <num> <num>  <num> <num> <num> <num> <num>  <num>  <num>
#>  1:     0     1      0     1     0     0     1      0      0
#>  2:     0     1      1     1     0     1     1      1      0
#>  3:     0     1      2     1     0     1     1      1      0
#>  4:     1     2      2     1     0     1     1      1      0
#>  5:     1     2      1     1     0     0     1      0      0
#>  6:     1     2      0     1     0     0     1      0      0
#>  7:     0     1      1     0     1     0     1      0      0
#>  8:     0     1      1     0     1     1     1      0      1
#>  9:     0     1      0     0     1     1     1      0      1
#> 10:     1     2      0     0     1     1     1      0      1
#> 11:     1     2      1     0     1     0     1      0      0
#> 12:     1     2      1     0     1     0     1      0      0
```

These results could be used as input for any of the available regression
functions. Colossus also includes a wrapper function to directly call
the Poisson model regression function.

``` r
control <- list(
  "ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 5, "epsilon" = 1e-6,
  "deriv_epsilon" = 1e-6, "verbose" = 2
)
e <- PoisRunJoint(formula_list, df, control = control)
print(e)
#> |-------------------------------------------------------------------|
#> Final Results
#>    Covariate Subterm Central Estimate Standard Error 2-tail p-value
#>       <char>  <char>            <num>          <num>          <num>
#> 1:        t0    plin           -0.184          0.385          0.631
#> 2:    fac_e0  loglin            0.574          0.468          0.219
#> 3:    fac_e1  loglin           -1.035          1.009          0.305
#> 
#> Poisson Model Used
#> -2*Log-Likelihood: 20.891,  Deviation: 6.436,  AIC: 12.436,  BIC: 28.346
#> Iterations run: 10
#> maximum step size: 4.775e-06, maximum first derivative: 1.595e-06
#> Analysis did not converge, check convergence criteria or run further
#> Run finished in 0.017 seconds
#> |-------------------------------------------------------------------|
```
