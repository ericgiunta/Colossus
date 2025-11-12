# Prints a cox regression output clearly

`print.coxres` uses the list output from a regression, prints off a
table of results and summarizes the score and convergence.

## Usage

``` r
# S3 method for class 'coxres'
print(x, ...)
```

## Arguments

- x:

  result object from a regression, class coxres

- ...:

  can include the number of digits, named digit, or an unnamed integer
  entry assumed to be digits

## Value

return nothing, prints the results to console

## See also

Other Output and Information Functions:
[`System_Version()`](System_Version.md),
[`print.caseconres()`](print.caseconres.md),
[`print.coxresbound()`](print.coxresbound.md),
[`print.logitres()`](print.logitres.md),
[`print.poisres()`](print.poisres.md),
[`print.poisresbound()`](print.poisresbound.md)
