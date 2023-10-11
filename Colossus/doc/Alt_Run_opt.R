## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Colossus)
library(data.table)


## ----eval=FALSE---------------------------------------------------------------
#  Strat_Col <- "e"
#  e <- RunCoxRegression_STRATA(df, time1, time2, event, names, Term_n, tform, keep_constant,
#                               a_n, modelform, fir, der_iden, control,Strat_Col)

## ----eval=FALSE---------------------------------------------------------------
#  Strat_Col <- c("e")
#  e <-RunPoissonRegression_STRATA(df, pyr, event, names, Term_n, tform, keep_constant,
#                                  a_n, modelform, fir, der_iden, control,Strat_Col)

## ----eval=FALSE---------------------------------------------------------------
#  e <- RunCoxRegression_Basic(df, time1, time2, event, names,
#                              keep_constant, a_n, der_iden, control)

## ----eval=FALSE---------------------------------------------------------------
#  e <- RunCoxRegression_Single(df, time1, time2, event, names, Term_n, tform,
#                               a_n, modelform, fir, control)
#  
#  e <- RunPoissonRegression_Single(df, pyr, event, names, Term_n, tform,
#                                   a_n, modelform, fir, control)

## ----eval=FALSE---------------------------------------------------------------
#  df$censor <- (df$lung==0) #censoring column made
#  event <- "censor" #event type switched to censoring
#  
#  plot_options <- list("name"="run_2","verbose"=FALSE,"studyID"="studyID","age_unit"="years")
#  #modified plotting function used to get censoring weights
#  dft <- GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant,
#                       a_n, modelform, fir, control, plot_options) #generates a survival curve
#  t_ref <- dft$t
#  surv_ref <- dft$surv
#  t_c <- df$t1
#  cens_weight <- approx(t_ref, surv_ref, t_c,rule=2)$y
#  #the surviving proportions used as censoring weight
#  event <- "lung" #event switched back
#  
#  e <- RunCoxRegression_CR(df, time1, time2, event, names, Term_n, tform, keep_constant,
#                           a_n, modelform, fir, der_iden, control,cens_weight)

