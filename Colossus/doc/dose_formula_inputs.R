## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Colossus)
library(data.table)

## -----------------------------------------------------------------------------
Term_n <- c(0,1,1,0,0)
tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
keep_constant <- c(0,0,0,1,0)
a_n <- c(1,2,3,4,5)
names <- c("a","a","a","a","a")
val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names)
Term_n <- val$Term_n
tform <- val$tform
keep_constant <- val$keep_constant
a_n <- val$a_n
der_iden <- val$der_iden
names <- val$names

