# Functions for Plotting and Analysis

``` r
Sys.setenv("OMP_THREAD_LIMIT" = 1) # Reducing core use, to avoid accidental use of too many cores
library(Colossus)
library(data.table)
if (system.file(package = "survival") != "") {
  library(survival)
}
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
if (system.file(package = "ggplot2") != "") {
  library(ggplot2)
}
```

## Example Setup

We will use analysis of the lung dataset in the survival package to
visualize the different plotting methods available in Colossus.

``` r
if (system.file(package = "survival") != "") {
  data(cancer, package = "survival")
  cancer %>% setDT()
  df <- copy(cancer)
} else {
  status <- c(2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1)
  sex <- c(1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 1, 2, 1, 2)
  time <- c(306, 455, 1010, 210, 883, 1022, 310, 361, 218, 166, 170, 654, 728, 71, 567, 144, 613, 707, 61, 88, 301, 81, 624, 371, 394, 520, 574, 118, 390, 12, 473, 26, 533, 107, 53, 122, 814, 965, 93, 731, 460, 153, 433, 145, 583, 95, 303, 519, 643, 765, 735, 189, 53, 246, 689, 65, 5, 132, 687, 345, 444, 223, 175, 60, 163, 65, 208, 821, 428, 230, 840, 305, 11, 132, 226, 426, 705, 363, 11, 176, 791, 95, 196, 167, 806, 284, 641, 147, 740, 163, 655, 239, 88, 245, 588, 30, 179, 310, 477, 166, 559, 450, 364, 107, 177, 156, 529, 11, 429, 351, 15, 181, 283, 201, 524, 13, 212, 524, 288, 363, 442, 199, 550, 54, 558, 207, 92, 60, 551, 543, 293, 202, 353, 511, 267, 511, 371, 387, 457, 337, 201, 404, 222, 62, 458, 356, 353, 163, 31, 340, 229, 444, 315, 182, 156, 329, 364, 291, 179, 376, 384, 268, 292, 142, 413, 266, 194, 320, 181, 285, 301, 348, 197, 382, 303, 296, 180, 186, 145, 269, 300, 284, 350, 272, 292, 332, 285, 259, 110, 286, 270, 81, 131, 225, 269, 225, 243, 279, 276, 135, 79, 59, 240, 202, 235, 105, 224, 239, 237, 173, 252, 221, 185, 92, 13, 222, 192, 183, 211, 175, 197, 203, 116, 188, 191, 105, 174, 177)
  age <- c(74, 68, 56, 57, 60, 74, 68, 71, 53, 61, 57, 68, 68, 60, 57, 67, 70, 63, 56, 57, 67, 49, 50, 58, 72, 70, 60, 70, 53, 74, 69, 73, 48, 60, 61, 62, 65, 66, 74, 64, 70, 73, 59, 60, 68, 76, 74, 63, 74, 50, 72, 63, 68, 58, 59, 62, 65, 57, 58, 64, 75, 48, 73, 65, 69, 68, 67, 64, 68, 67, 63, 48, 74, 40, 53, 71, 51, 56, 81, 73, 59, 55, 42, 44, 44, 71, 62, 61, 44, 72, 63, 70, 66, 57, 69, 72, 69, 71, 64, 70, 58, 69, 56, 63, 59, 66, 54, 67, 55, 75, 69, 44, 80, 75, 54, 76, 49, 68, 66, 80, 75, 60, 69, 72, 70, 66, 50, 64, 77, 48, 59, 53, 47, 55, 67, 74, 58, 56, 54, 56, 73, 74, 76, 65, 57, 53, 71, 54, 82, 59, 70, 60, 62, 53, 55, 69, 68, 62, 63, 56, 62, 44, 69, 63, 64, 57, 60, 46, 61, 65, 61, 58, 56, 43, 53, 59, 56, 55, 53, 74, 60, 39, 66, 65, 51, 45, 72, 58, 64, 53, 72, 52, 50, 64, 71, 70, 63, 64, 52, 60, 64, 73, 63, 50, 63, 62, 55, 50, 69, 59, 60, 67, 69, 64, 65, 65, 41, 76, 70, 57, 67, 71, 76, 77, 39, 75, 66, 58)

  df <- data.table(
    "status" = status,
    "sex" = sex,
    "time" = time,
    "age" = age
  )
}

df$UserID <- seq_len(nrow(df))

df$status <- df$status - 1
df$sex <- df$sex - 1

control <- list(ncore = 1)
a_n <- c(0.01701289, -0.51256478)
coxres <- CoxRun(Cox(time, status) ~ loglinear(age, sex, 0),
  df,
  control = control, a_n = a_n
)
```

## Survival Function Approximation

After fitting a Cox proportional hazards model, one may be interested in
what the baseline survival rate is. One method of doing so is by
weighting the number of events at each event time by the total risk. In
the absence of excess risk, the hazard at any point in time is equal to
the event rate and the hazard ratio of every row is equal to one, so
this assumption holds for the case with negligible excess risk.

$$\begin{array}{r}
{\lambda(t) = \frac{d\lbrack t\rbrack}{n\lbrack t\rbrack} = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}(1)} = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}\left( r_{i} \right)}}
\end{array}$$

Suppose every row at risk at a time was twice as likely to experience an
event than the baseline. Then we would expect there to be twice as many
events as there would be at baseline. The same logic applies to the case
with every row at risk being half as likely. This generalizes to any
average risk. Colossus allows the instantaneous hazard to be
approximated for both stratified and non-stratified models.

$$\begin{array}{r}
{\lambda(t) = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}\left( r_{i} \right)} = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}(2)} = \frac{d\lbrack t\rbrack/2}{n\lbrack t\rbrack}} \\
{\lambda(t) = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}\left( r_{i} \right)} = \frac{d\lbrack t\rbrack}{\sum\limits_{i = 1}^{n}(0.5)} = \frac{d\lbrack t\rbrack*2}{n\lbrack t\rbrack}}
\end{array}$$

Once the instantaneous hazard is approximated, the cumulative hazard can
be approximated. Then the surviving fraction is approximately equal to
the exponential of the negative cumulative hazard at each event time.

$$\begin{array}{r}
{\Lambda(t) = \int_{0}^{t}\left( \lambda(x)dx \right)} \\
{S(t) = \exp\left( - \Lambda(t) \right)}
\end{array}$$

``` r
plot_options <- list(
  "fname" = paste(tempfile(), "run", sep = ""), "studyid" = "UserID",
  "verbose" = 2, "surv_curv" = T, "martingale" = F, "strat_haz" = F, "km" = F
)

e <- plotSurvival(coxres, df, plot_options)

norm_surv <- e[["standard"]]
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(norm_surv, aes(x = .data$t, y = .data$h)) +
    geom_point(color = "black") +
    labs(x = "age", y = "Instantaneous Hazard")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-3-1.png)

``` r
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(norm_surv, aes(x = .data$t, y = .data$ch)) +
    geom_line(color = "black", alpha = 1) +
    labs(x = "age", y = "Cumulative Hazard")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-3-2.png)

``` r
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(norm_surv, aes(x = .data$t, y = .data$surv)) +
    geom_line(color = "black", alpha = 1) +
    labs(x = "age", y = "Surviving Fraction")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-3-3.png)

``` r
plot_options <- list(
  "fname" = paste(tempfile(), "run", sep = ""), "studyid" = "UserID",
  "verbose" = 2, "surv_curv" = F, "martingale" = F, "strat_haz" = F, "km" = T
)

e <- plotSurvival(coxres, df, plot_options)

km <- e[["kaplin-meier"]]
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(km, aes(x = .data$t_t, y = .data$n_t)) +
    geom_line(color = "black", alpha = 1) +
    labs(x = "age", y = "KM Survival")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-3-4.png)

## Cox Proportional Hazards model assumptions

The Cox proportional hazards model by definition assumes proportional
hazards independent of time. If this is violated, then the results of a
regression may be misleading. There are two checks that Colossus
provides that can be used to test this assumption, the Schoenfeld
residuals and Martingale residuals. In both cases, the premise is that
if the hazard ratio were independent of time, then the residuals should
also be independent of time.

### Schoenfeld Residuals

Schoenfeld residuals compare the average covariate value of rows with
events with the risk-weighted average of the covariate in rows at risk.
Consistently high or low residuals may be due to the hazard being much
higher or lower than the model predicts. If these residuals are also
correlated with event time, then the hazard ratio may be dependent on
event time. There is also the option to scale the residuals by the
standard deviation.

$$\begin{array}{r}
{s(t,x) = \frac{\sum\limits_{i \in events}x_{i}}{d\lbrack t\rbrack} - \frac{\sum\limits_{i}^{n}x_{i}*r_{i}}{\sum\limits_{i}^{n}r_{i}}}
\end{array}$$

``` r
plot_options <- list(
  "fname" = paste(tempfile(), "run", sep = ""),
  "studyid" = "UserID", "verbose" = 2
)

res_all <- plotSchoenfeld(coxres, df, plot_options)

res_age <- res_all[["age"]]
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_age, aes(x = .data$time, y = .data$y)) +
    geom_point(color = "black") +
    labs(
      x = paste("Survival Time", sep = ""),
      y = paste("Schoenfeld Residual (age)", sep = " ")
    )
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-4-1.png)

``` r
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_age, aes(x = .data$time, y = .data$y_scale)) +
    geom_point(color = "black") +
    labs(
      x = paste("Survival Time", sep = ""),
      y = paste("Schoenfeld Residual Scaled (age)", sep = " ")
    )
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-4-2.png)

``` r

res_sex <- res_all[["sex"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_sex, aes(x = .data$time, y = .data$y)) +
    geom_point(color = "black") +
    labs(
      x = paste("Survival Time", sep = ""),
      y = paste("Schoenfeld Residual (sex)", sep = " ")
    )
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-4-3.png)

``` r
if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_sex, aes(x = .data$time, y = .data$y_scale)) +
    geom_point(color = "black") +
    labs(
      x = paste("Survival Time", sep = ""),
      y = paste("Schoenfeld Residual Scaled (sex)", sep = " ")
    )
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-4-4.png)

### Martingale Residuals

Martingale residuals compare the event status and cumulative hazard for
each subject. The residual is bound in the open interval
$( - \infty,1)$. Negative residuals correspond with subjects that
survive despite high cumulative hazard and residuals near unity
correspond with subjects that experienced an event despite low
cumulative hazard. The distribution of residuals with event time can
indicate if the model is over-predicting or under-predicting dependent
on time.

$$\begin{array}{r}
{m_{j} = \delta_{j} - \int_{t_{0}}^{t_{1}}\lambda(t)*r_{j}dt}
\end{array}$$

``` r
plot_options <- list(
  "fname" = paste(tempfile(), "run", sep = ""),
  "studyid" = "UserID", "verbose" = 2, "surv_curv" = F,
  "martingale" = T, "strat_haz" = F, "km" = F, "cov_cols" = c("age", "sex")
)
res_all <- plotMartingale(coxres, df, plot_options)

res_age <- res_all[["age"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot() +
    geom_point(
      data = res_age,
      aes(x = .data$cov_max, y = .data$res_sum, group = .data$event, color = .data$event)
    )
  g <- g + labs(x = "Max Age", y = "Martingale Residuals")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-5-1.png)

``` r
res_sex <- res_all[["sex"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot() +
    geom_point(
      data = res_sex,
      aes(x = .data$cov_max, y = .data$res_sum, group = .data$event, color = .data$event)
    )
  g <- g + labs(x = "Sex", y = "Martingale Residuals")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-5-2.png)

``` r
res_surv <- res_all[["survival_time"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot() +
    geom_point(
      data = res_surv,
      aes(x = .data$time_max, y = .data$res_sum, group = .data$event, color = .data$event)
    )
  g <- g + labs(x = "Survival Time", y = "Martingale Residuals")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-5-3.png)

## General Evaluation of Risk

Colossus also offers scripts that plot the relative risk by covariate
value for each covariate in the model, assuming every other covariate is
zero. Results for specific covariates can be requested by passing an
optional argument, otherwise every covariate used in the risk model is
included. Risk at Wald confidence interval boundaries can also be
calculated by including a ‘boundary’ option set to a z-score. In this
example, 95 percent confidence intervals are included and plotted.

``` r
plot_options <- list(
  "fname" = paste(tempfile(), "run", sep = ""), "studyid" = "UserID",
  "verbose" = 2, "cov_cols" = c("age", "sex"), boundary = 1.96
)
res_all <- plotRisk(coxres, df, plot_options)

res_age <- res_all[["age"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_age, aes(x = .data$x, y = .data$y)) +
    geom_errorbar(aes(ymin = .data$"y:lower", ymax = .data$"y:upper"), color = "black") +
    geom_point(color = "black") +
    labs(x = "Age", y = "Relative Risk")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-6-1.png)

``` r


res_sex <- res_all[["sex"]]

if (system.file(package = "ggplot2") != "") {
  g <- ggplot(res_sex, aes(x = .data$x, y = .data$y)) +
    geom_errorbar(aes(ymin = .data$"y:lower", ymax = .data$"y:upper"), color = "black") +
    geom_point(color = "black") +
    labs(x = "Sex", y = "Relative Risk")
} else {
  g <- message("ggplot2 wasn't detected. Please install to see the plot")
}
g
```

![](Plotting_And_Analysis_files/figure-html/unnamed-chunk-6-2.png)
