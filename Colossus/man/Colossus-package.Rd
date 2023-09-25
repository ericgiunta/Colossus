\name{Colossus-package}
\alias{Colossus-package}
\alias{Colossus}
\docType{package}
\title{
  Risk Model Regression and Analysis with Complex Non-Linear Models
}
\description{
  Performs survival analysis using general non-linear models. Risk models can be the sum or product of terms. Each term is the product of exponential/linear functions of covariates. Additionally sub-terms can be defined as a sum of exponential, linear threshold, and step functions. Cox Proportional hazards, Poisson, and Fine-Grey competing risks regression are supported.
}
\details{
  Colossus requires a data.table with columns defining time duration, covariates used in the risk model, and a column with event status/counts. Cox Proportional Hazards and Fine-Grey models require two columns to set the interval start and end points, one may be specified as truncated to set one-sided intervals. Poisson models require a single column for person-years. The event column is assumed to contain 0/1 for Cox Proportional Hazards models for censored/event, 0/1/2 for Fine-Grey models for censored/main event/competing event, and a count of events for Poisson models.
  Risk models in Colossus are composed of subterms multiplied to calculate terms, then the risk is calculated as the product and/or sum of terms. Please see the starting_description vignette for more detail. The model is passed to Colossus as a collection of vectors which provide the column, subterm type, term number, parameter value, and fixed/free status for every model element.
  Colossus also uses a series of control lists to control different aspects of the regression/calculation. Generally each function has a control list that provides hardware/convergence values, control_guess to control distributed starting point generation if needed, and model_control to control stratification/competing risks/etc additional model options. Please read the Alt_Distrib_Starts and Alt_Run_opt vignettes for more details.
  Generally all of the Cox Proportional Hazards functions are wrappers which call RunCoxRegression_Omnibus with different options, similarly all of the Poisson functions call RunPoissonRegression_Omnibus. All of the plotting and assumptions tests are accessed through the RunCoxPlots function.
}
\author{
Eric Giunta, egiunta@ksu.edu

Maintainer: Eric Giunta <egiunta@ksu.edu>
}
\keyword{ package }
\examples{
  library(data.table)
  ## basic example code reproduced from the starting-description vignette

  df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
              "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
                "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
            "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
                        "a"=c(0,   1,   1,   0,   1,   0,   1),
                        "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
                        "c"=c(10,  11,  10,  11,  12,  9,   11),
                        "d"=c(0,   0,   0,   1,   1,   1,   1),
                        "e"=c(0,   0,   1,   0,   0,   0,   1))
  # For the interval case
  time1 <- "Starting_Age"
  time2 <- "Ending_Age"
  pyr <- "Ending_Age"
  event <- "Cancer_Status"
  names <- c('a','b','c','d')
  a_n <- list(c(1.1, -0.1, 0.2, 0.5),c(1.6, -0.12, 0.3, 0.4))
  #used to test at a specific point
  Term_n <- c(0,1,1,2)
  tform <- c("loglin","lin","lin","plin")
  modelform <- "M"
  fir <- 0

  keep_constant <- c(0,0,0,0)
  der_iden <- 0

  control <- list("Ncores"=2,'lr' = 0.75,'maxiters' = c(5,5,5),
      'halfmax' = 5,'epsilon' = 1e-3, 'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
      'abs_max'=1.0,'change_all'=TRUE, 'dose_abs_max'=100.0,'verbose'=FALSE,
      'ties'='breslow','double_step'=1, "guesses"=2)

  e <- RunCoxRegression_Omnibus(df, time1, time2, event,
                                names, Term_n, tform, keep_constant,
                                a_n, modelform, fir, der_iden, control,
                                model_control=list("single"=FALSE,
                                "basic"=FALSE, "CR"=FALSE, 'null'=FALSE))
  Strat_Col <- 'e'
  e <- RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant,
                                    a_n, modelform, fir, der_iden, control,Strat_Col)
  plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")), "studyID"="UserID",
                    'verbose'=FALSE)

  RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant,
              a_n, modelform, fir, control, plot_options)
}
