# Calculates The Additional Parameter For a linear-exponential formula with known maximum

`Linked_Lin_Exp_Para` Calculates what the additional parameter would be
for a desired maximum

## Usage

``` r
Linked_Lin_Exp_Para(y, a0, a1_goal, verbose = 0)
```

## Arguments

- y:

  point formula switch

- a0:

  linear slope

- a1_goal:

  exponential maximum desired

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

returns parameter used by Colossus

## Examples

``` r
library(data.table)
y <- 7.6
a0 <- 1.2
a1_goal <- 15
full_paras <- Linked_Lin_Exp_Para(y, a0, a1_goal)
```
