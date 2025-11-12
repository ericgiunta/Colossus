# Unified Equation Representation

``` r
library(Colossus)
#> Note: From versions 1.3.1 to 1.4.1 the expected inputs changed. Regressions are now run with CoxRun and PoisRun and formula inputs. Please see the 'Unified Equation Representation' vignette for more details.
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
\text{Surv(interval, event) ~ name + ... + factor(name)} \\
\text{’Survival(interval, event) ~ subterm(name, factor(name), term\_number) + ... + term\_model()’}
\end{array}$$

Based on the number of entries, the format of the time interval is
assumed. Cox models used an entry, exit, and event time. These can be
either named ‘start’, ‘tend’, and ‘event’ or included in order. If none
are named, a left truncated interval is assumed. A stratified cox model
includes the named entry ‘strata’ or the right-most entry is assumed to
be the stratification column. A Fine-Gray model adds a weighting column
to the standard Cox interval representation, assumed to be an entry
named ‘weight’ or the right-most entry. In total, the assumed order is:
entry time, exit time, strata column, and weighting column.

$$\begin{array}{r}
\text{’Cox(tstart, tend, event) ~ ...’} \\
\text{’Cox\_Strata(tstart, tend, event, strata) ~ ...’} \\
\text{’FineGray(tstart, tend, event, weight) ~ ...’} \\
\text{’FineGray\_Strata(tstart, tend, event, strata, weight) ~ ...’}
\end{array}$$

Poisson models use a measure of duration (named ‘pyr’), the event column
(named ‘event’), and any strata columns. The assumed order is ‘pyr’,
‘event’, and the strata columns listed sequentially. The duration and
event columns can be named, but the naming of the strata columns is
currently not supported.

$$\begin{array}{r}
\text{’Pois(pyr, event) ~ ...’} \\
\text{’Pois\_Strata(pyr, event, strata\_0) ~ ...’} \\
\text{’Pois\_Strata(pyr, event, strata\_0, strata\_1, ...) ~ ...’}
\end{array}$$

The right hand side adds elements to subterms. Each subterm should have
the included columns, and optionally the term number. Subterms are
assumed to be in term 0 is not provided, and the model is assumed to be
multiplicative. The factor option automatically removes the first level
as a reference.

The following expressions are equivalent.

``` r
term_n <- c(0, 0, 1)
tform <- c("loglin", "loglin", "lin")
names <- c("dose0", "dose1", "dose2")
modelform <- "M"
tstart <- "t0"
tend <- "t1"
event <- "lung"

Model_Eq <- Cox(t0, t1, lung) ~ loglinear(dose0, dose1, 0) + linear(dose2, 1) + multiplicative()
Model_Eq <- Cox(t0, t1, lung) ~ loglinear(dose0, dose1) + linear(dose2, 1)

df <- data.table(
  "dose0" = 1:4, "dose1" = 2:5, "dose2" = 3:6,
  "t0" = c(0, 0, 1, 0), "t1" = 5:8, "lung" = c(1, 0, 0, 1)
)
get_form(Model_Eq, df)
#> $model
#> $start_age
#> [1] "t0"
#> 
#> $end_age
#> [1] "t1"
#> 
#> $event
#> [1] "lung"
#> 
#> $strata
#> [1] "NONE"
#> 
#> $weight
#> [1] "NONE"
#> 
#> $null
#> [1] FALSE
#> 
#> $term_n
#> [1] 0 0 1
#> 
#> $tform
#> [1] "loglin" "loglin" "lin"   
#> 
#> $names
#> [1] "dose0" "dose1" "dose2"
#> 
#> $a_n
#> [1] 0.01 0.01 0.01
#> 
#> $keep_constant
#> [1] 0 0 0
#> 
#> $modelform
#> [1] "M"
#> 
#> $gmix_term
#> NULL
#> 
#> $gmix_theta
#> [1] 0
#> 
#> $expres_calls
#> list()
#> 
#> attr(,"class")
#> [1] "coxmodel"
#> 
#> $data
#>    dose0 dose1 dose2    t0    t1  lung
#>    <int> <int> <int> <num> <int> <num>
#> 1:     1     2     3     0     5     1
#> 2:     2     3     4     0     6     0
#> 3:     3     4     5     1     7     0
#> 4:     4     5     6     0     8     1
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
|                   loglin                    |                    “loglin”, “loglinear”, “log-linear”                    |
|           loglin_slope/loglin_top           |            “loglin-dose”, “loglinear-dose”, “log-linear-dose”             |
|              lin_slope/lin_top              |               “lin-dose”, “linear-dose”, “linear-piecewise”               |
|                 quad_slope                  |            “quadratic”, “quad”, “quad-dose”, “quadratic-dose”             |
|             step_slope/step_int             |                       “step-dose”, “step-piecewise”                       |
|         lin_quad_slope/lin_quad_int         |  “lin-quad-dose”, “linear-quadratic-dose”, “linear-quadratic-piecewise”   |
| lin_exp_slope/lin_exp_int/lin_exp_exp_slope | “lin-exp-dose”, “linear-exponential-dose”, “linear-exponential-piecewise” |

| Term Type |                  Equivalent Aliases                  |
|:---------:|:----------------------------------------------------:|
|     M     | “m”, “me”, “multiplicative”, “multiplicative-excess” |
|     A     |                   “a”, “additive”                    |
|    PA     |               “pa”, “product-additive”               |
|    PAE    |           “pae”, “product-additive-excess”           |
|   GMIX    |              “gmix”,“geometric-mixture”              |
|  GMIX-R   |        “gmix-r”,“relative-geometric-mixture”         |
|  GMIX-E   |         “gmix-e”,“excess-geometric-mixture”          |
