test_that("Coxph lin_err check no error", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  df$lung <- (df$lung > 0)
  names <- c("dose", "rand", "rand0", "rand1", "rand2")
  term_n <- c(0, 0, 0, 0, 0)
  tform <- c("plin", "loglin", "loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- rep(0.1, length(names))
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  #
  event <- "lung"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  model_control <- list("linear_err" = T)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = "M", fir = 0, der_iden = der_iden, control = control, model_control = model_control))
})
test_that("Coxph lin_err check value", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  df$lung <- (df$lung > 0)
  names <- c("dose", "rand", "rand0", "rand1", "rand2")
  term_n <- c(0, 0, 0, 0, 0)
  tform <- c("plin", "loglin", "loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- rep(0.1, length(names))
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  #
  event <- "lung"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  for (i in 1:5) {
      a_n <- runif(length(names))
      tform <- c("loglin", "loglin", "loglin", "loglin", "loglin")
      tform[i] <- "plin"
      model_control <- list("linear_err" = T)
      e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = "M", fir = 0, der_iden = der_iden, control = control, model_control = model_control)
      model_control <- list("linear_err" = F)
      e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = "M", fir = 0, der_iden = der_iden, control = control, model_control = model_control)
      expect_equal(e0$LogLik, e1$LogLik)
      expect_equal(e0$First_Der, e1$First_Der)
      expect_equal(e0$Second_Der, e1$Second_Der)
  }
})
test_that("Coxph basic check value", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  df$lung <- (df$lung > 0)
  names <- c("dose", "rand", "rand0", "rand1", "rand2")
  term_n <- c(0, 0, 0, 0, 0)
  tform <- c("plin", "loglin", "loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- rep(0.1, length(names))
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  #
  event <- "lung"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  for (i in 1:5) {
      a_n <- runif(length(names))
      tform <- c("loglin", "loglin", "loglin", "loglin", "loglin")
      model_control <- list("basic" = T)
      e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = "M", fir = 0, der_iden = der_iden, control = control, model_control = model_control)
      model_control <- list("basic" = F)
      e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = "M", fir = 0, der_iden = der_iden, control = control, model_control = model_control)
      expect_equal(e0$LogLik, e1$LogLik)
      expect_equal(e0$First_Der, e1$First_Der)
      expect_equal(e0$Second_Der, e1$Second_Der)
  }
})
