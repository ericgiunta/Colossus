# Prints a poisson likelihood boundary regression output clearly

`print.poisresbound` uses the list output from a regression, prints off
a table of results and summarizes the score and convergence.

## Usage

``` r
# S3 method for class 'poisresbound'
print(x, ...)
```

## Arguments

- x:

  result object from a regression, class poisresbound

- ...:

  can include the number of digits, named digit, or an unnamed integer
  entry assumed to be digits

## Value

return nothing, prints the results to console

## See also

Other Output and Information Functions:
[`System_Version()`](System_Version.md),
[`print.caseconres()`](print.caseconres.md),
[`print.coxres()`](print.coxres.md),
[`print.coxresbound()`](print.coxresbound.md),
[`print.logitres()`](print.logitres.md),
[`print.poisres()`](print.poisres.md)
