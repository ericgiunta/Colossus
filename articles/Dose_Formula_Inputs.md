# Dose Response Formula Terms

``` r
library(Colossus)
library(data.table)
```

## Dose Response Formula

As previously discussed, Colossus features a term composed of the sum of
multiple linear and non-linear elements which can be used to define many
dose-response curves used in radiation epidemiology. These terms are
referred to as dose-response terms, but there is nothing prohibiting
them from being used for non-dose covariates. The following formulae are
available, reproduced from the starting description vignette.

$$\begin{array}{r}
{S_{NL} = \sum\limits_{i}\left( \alpha_{i} \times \exp\left( x_{i} \cdot \beta_{i} \right) \right) + \sum\limits_{i}\left( \beta_{i} \cdot \left( x_{i} \right)^{2} \right) + \sum\limits_{i}F_{LT} + \sum\limits_{i}F_{STP} + \sum\limits_{i}F_{LQ} + \sum\limits_{i}F_{LEXP}} \\
{F_{LT} = \begin{cases}
{\alpha_{i} \cdot \left( x - \beta_{i} \right)} & \left( x > \beta_{i} \right) \\
0 & \text{else}
\end{cases}} \\
{F_{STP} = \begin{cases}
\alpha_{i} & \left( x > \beta_{i} \right) \\
0 & \text{else}
\end{cases}} \\
{F_{LQ} = \begin{cases}
{\beta_{i} \cdot x} & \left( x > \alpha_{i} \right) \\
{\lambda_{i} \cdot x^{2} + \nu_{i}} & \text{else}
\end{cases}} \\
{F_{LEXP} = \begin{cases}
{\beta_{i} \cdot x} & \left( x > \alpha_{i} \right) \\
{\lambda_{i} - \exp\left( \nu_{i} + \mu \cdot x \right)} & \text{else}
\end{cases}} \\
{T_{j} = S_{LL,j} \times S_{L,j} \times S_{PL,j} \times S_{NL,j}}
\end{array}$$

For every subterm type, there are between 1 and 3 parameters that fully
define the curve. The Linear-Quadratic and Linear-Exponential curves are
continuously differentiable, so there are only 2-3 parameters that can
be set.

$$\begin{array}{r}
{\lambda_{LQ} = \beta_{LQ}/\left( 2\alpha_{LQ} \right)} \\
{\nu_{LQ} = \left( \beta_{LQ}*\alpha_{LQ} \right)/2} \\
{\nu_{LEXP} = \ln\left( \beta_{LEXP} \right) - \ln\left( \mu_{LEXP} \right) + \mu_{LEXP}*\alpha_{LEXP}} \\
{\lambda_{LEXP} = \beta_{LEXP}*\alpha_{LEXP} + exp\left( \nu_{LEXP} - \mu_{LEXP}*\alpha_{LEXP} \right)}
\end{array}$$

## Using The Different subterms

|      Subterm       |    tform Entry    |                                    Description                                     |
|:------------------:|:-----------------:|:----------------------------------------------------------------------------------:|
|    Exponential     |    loglin_top     |                 parameter in the exponent of the term, $\beta_{i}$                 |
|    Exponential     |   loglin_slope    | parameter multiplied by the exponential assumed to be 1 if not given, $\alpha_{i}$ |
|  Linear Threshold  |     lin_slope     |                      slope for the linear term, $\alpha_{i}$                       |
|  Linear Threshold  |      lin_int      |                     intercept for the linear term, $\beta_{i}$                     |
|   Step Function    |    step_slope     |                         step function value, $\alpha_{i}$                          |
|   Step Function    |     step_int      |                        step function intercept, $\beta_{i}$                        |
|     Quadratic      |    quad_slope     |               parameter multiplied by the squared value, $\beta_{i}$               |
| Linear-Exponential |   lin_exp_slope   |                           Linear slope term, $\beta_{i}$                           |
| Linear-Exponential |    lin_exp_int    |               Intercept between linear to exponential, $\alpha_{i}$                |
| Linear-Exponential | lin_exp_exp_slope |                      Slope term in the exponential, $\mu_{i}$                      |
|  Linear-Quadratic  |  lin_quad_slope   |                           Linear slope term, $\beta_{i}$                           |
|  Linear-Quadratic  |   lin_quad_int    |                Intercept between linear to quadratic, $\alpha_{i}$                 |

The linear-exponential and linear-quadratic curves must be either
completely fixed or completely free. In contrast, the exponential,
linear threshold, and step-function curves can be partially fixed. The
exponential term can be provided with only the covariate in the exponent
and assume the magnitude to be 1. The linear threshold and step
functions can be provided a fixed threshold covariate, which can be used
to define a linear-no-threshold model or a combination of linear and
step functions with a known threshold.
