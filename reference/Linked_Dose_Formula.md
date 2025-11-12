# Calculates Full Parameter list for Special Dose Formula

`Linked_Dose_Formula` Calculates all parameters for linear-quadratic and
linear-exponential linked formulas

## Usage

``` r
Linked_Dose_Formula(tforms, paras, verbose = 0)
```

## Arguments

- tforms:

  list of formula types

- paras:

  list of formula parameters

- verbose:

  integer valued 0-4 controlling what information is printed to the
  terminal. Each level includes the lower levels. 0: silent, 1: errors
  printed, 2: warnings printed, 3: notes printed, 4: debug information
  printed. Errors are situations that stop the regression, warnings are
  situations that assume default values that the user might not have
  intended, notes provide information on regression progress, and debug
  prints out C++ progress and intermediate results. The default level is
  2 and True/False is converted to 3/0.

## Value

returns list of full parameters

## Examples

``` r
library(data.table)
tforms <- list("cov_0" = "quad", "cov_1" = "exp")
paras <- list("cov_0" = c(1, 3.45), "cov_1" = c(1.2, 4.5, 0.1))
full_paras <- Linked_Dose_Formula(tforms, paras)
```
