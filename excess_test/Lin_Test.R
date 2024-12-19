library(Colossus)
library(data.table)
library(parallel)

if (identical(Sys.getenv("TESTTHAT"), "true")){
} else {
  fname <- "tests/testthat/base_example.csv"
  df <- fread(fname)

  time1 <- "entry"
  time2 <- "exit"
  event <- "event"
  names <- c("dose0", "dose1", "dose0")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "lin")
  keep_constant <- c(0, 0, 0)
  # a_n <- c(0.2462, 5.020, -0.5909)
  a_n <- c(-1.493177, 5.020007, 1.438377)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  #
  alpha <- 0.005
  alpha_list <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 4, "ties" = "breslow", "double_step" = 1, "guesses" = 10)
  for (alpha_i in 1:length(alpha_list)) {
    alpha <- alpha_list[alpha_i]
    a_n <- c(-1.493177, 5.020007, 1.438377)
    model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = TRUE)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "nan", model_control = model_control)
  }
}
