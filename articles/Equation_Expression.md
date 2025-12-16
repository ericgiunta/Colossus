# Unified Equation Representation

``` r
Sys.setenv("OMP_THREAD_LIMIT" = 1) # Reducing core use, to avoid accidental use of too many cores
library(Colossus)
library(data.table)
```

The underlying functions in Colossus specify a survival model using a
series of vectors, which can be overly complicated at times. Colossus
includes functions to interpret formulas, and then call the underlying
functions.

## Equation Expression

Generally, a survival model is defined by using a Surv object,
initialized with the interval and event columns, and a list of columns.
Because Colossus can handle much more complicated models, the definition
has to be more complicated. Instead of a Surv object, the model is
specified with ‘cox’, ‘poisson’, or ‘finegray’, each also including a
strata option. The right side of the equation is listed similarly.
Columns are added to subterm types (loglinear, linear, plinear, etc.).
Each subterm has a list of columns and a term number. Finally, the
method to combine terms is added as a separate item (multiplicative,
additive, etc.).

$$\begin{array}{r}
{{\text{Surv(interval, event)}\mspace{6mu}} \sim \text{name + ... + factor(name)}} \\
{{\text{’Survival(interval, event)}\mspace{6mu}} \sim {\mspace{6mu}\text{subterm(name, factor(name), term\_number) + ... + term\_model()’}}}
\end{array}$$

Based on the number of entries, the format of the time interval is
assumed. Cox models use an entry time, an exit time, and an event
status. These can be either named ‘start’, ‘tend’, and ‘event’ or
included in order. If none are named, a left truncated interval is
assumed. A stratified cox model includes the named entry ‘strata’ or the
right-most entry is assumed to be the stratification column. A Fine-Gray
model adds a weighting column to the standard Cox interval
representation, assumed to be an entry named ‘weight’ or the right-most
entry. In total, the assumed order is: entry time, exit time, event
status, strata column, and weighting column.

$$\begin{array}{r}
{{\text{’Cox(tstart, tend, event)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}} \\
{{\text{’Cox\_Strata(tstart, tend, event, strata)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}} \\
{{\text{’FineGray(tstart, tend, event, weight)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}} \\
{{\text{’FineGray\_Strata(tstart, tend, event, strata, weight)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}}
\end{array}$$

Poisson models use a measure of duration (named ‘pyr’), the number of
events (named ‘event’), and any strata columns. The assumed order is
‘pyr’, ‘event’, and the strata columns listed sequentially. The duration
and event columns can be named, but the naming of the strata columns is
currently not supported.

$$\begin{array}{r}
{{\text{’Pois(pyr, event)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}} \\
{{\text{’Pois\_Strata(pyr, event, strata\_0)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}} \\
{{\text{’Pois\_Strata(pyr, event, strata\_0, strata\_1, ...)}\mspace{6mu}} \sim {\mspace{6mu}\text{...’}}}
\end{array}$$

The right hand side adds elements to subterms. Each subterm should have
the included columns, and optionally the term number. Subterms are
assumed to be in term 0 if a term number is not provided, and the model
is assumed to be multiplicative excess. The factor option automatically
removes the first level as a reference. Note that unused levels can be
listed first to use all levels.

The following expressions are equivalent.

``` r
term_n <- c(0, 0, 1)
tform <- c("loglin", "loglin", "lin")
names <- c("dose0", "dose1", "dose2")
modelform <- "M"
tstart <- "t0"
tend <- "t1"
event <- "lung"

Model_Eq <- Cox(t0, t1, lung) ~ loglinear(dose0, dose1, 0) +
  linear(dose2, 1) + ME()
Model_Eq <- Cox(t0, t1, lung) ~ loglinear(dose0, dose1) +
  linear(dose2, 1)

df <- data.table(
  "dose0" = 1:4, "dose1" = 2:5, "dose2" = 3:6,
  "t0" = c(0, 0, 1, 0), "t1" = 5:8, "lung" = c(1, 0, 0, 1)
)
res <- get_form(Model_Eq, df, nthreads = 1)
formula <- res$model
new_data <- res$data
```

The following tables cover the subterms and term models that are
currently implemented. Every option listed has multiple equivalent
aliases that all refer to the same subterm or model formula.

|                      Survival Type                      |                                                                     Equivalent Aliases                                                                     |
|:-------------------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------:|
|                           Cox                           |                                                                       “cox”, “coxph”                                                                       |
|                     Cox stratified                      |                                                                “cox_strata”, “coxph_strata”                                                                |
|                         Poisson                         |                                                                     “poisson”, “pois”                                                                      |
|                   Poisson stratified                    |                                                              “poisson_strata”, “pois_strata”                                                               |
|                        Fine-Gray                        |                                                                      “finegray”, “fg”                                                                      |
|                  Fine-Gray stratified                   |                                                               “finegray_strata”, “fg_strata”                                                               |
|                        Logistic                         |                                                                    “logit”, “logistic”                                                                     |
|           Matched Case-Control with one group           |                                                          “casecon”, “casecontrol”, “case_control”                                                          |
|      Matched Case-Control grouped by time at risk       |                                                  “casecon_time”, “casecontrol_time”, “case_control_time”                                                   |
|         Matched Case-Control grouped by strata          |                                               “casecon_strata”, “casecontrol_strata”, “case_control_strata”                                                |
| Matched Case-Control grouped by time at risk and strata | “casecon_strata_time”, “casecontrol_strata_time”, “case_control_strata_time”, “casecon_time_strata”, “casecontrol_time_strata”, “case_control_time_strata” |

|                Subterm Type                 |                            Equivalent Aliases                             |
|:-------------------------------------------:|:-------------------------------------------------------------------------:|
|                    plin                     |                    “plin”, “plinear”, “product-linear”                    |
|                     lin                     |                              “lin”, “linear”                              |
|                   loglin                    |        “loglin”, “loglinear”, “log-linear” , “exponential”, “exp”         |
|           loglin_slope/loglin_top           |            “loglin-dose”, “loglinear-dose”, “log-linear-dose”             |
|              lin_slope/lin_top              |               “lin-dose”, “linear-dose”, “linear-piecewise”               |
|                 quad_slope                  |            “quadratic”, “quad”, “quad-dose”, “quadratic-dose”             |
|             step_slope/step_int             |                       “step-dose”, “step-piecewise”                       |
|         lin_quad_slope/lin_quad_int         |  “lin-quad-dose”, “linear-quadratic-dose”, “linear-quadratic-piecewise”   |
| lin_exp_slope/lin_exp_int/lin_exp_exp_slope | “lin-exp-dose”, “linear-exponential-dose”, “linear-exponential-piecewise” |

|                   Term Type                   |          Equivalent Aliases           |
|:---------------------------------------------:|:-------------------------------------:|
|                Multiplicative                 |         “m”, “multiplicative”         |
|             Multiplicative Excess             |     “me”, “multiplicative-excess”     |
|                   Additive                    |            “a”, “additive”            |
|               Product Additive                |       “pa”, “product-additive”        |
|            Product Additive Excess            |   “pae”, “product-additive-excess”    |
|            Geometric Mixture Model            |      “gmix”,“geometric-mixture”       |
| Geometric Mixture Model, relative error terms | “gmix-r”,“relative-geometric-mixture” |
|  Geometric Mixture Model, excess error terms  |  “gmix-e”,“excess-geometric-mixture”  |
