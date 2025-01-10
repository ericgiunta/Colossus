test_that("Pois_tier_guess various_fixes", {
  fname <- "MULTI_COV.csv"
  colTypes <- c("double", "double", "integer", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  event <- "lung"
  names <- c("t0", "a")
  term_n <- c(1, 2)
  tform <- c("loglin", "loglin")
  keep_constant <- c(0, 0)
  a_n <- c(0.01, -15)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 2, "guesses" = 2, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "strata" = TRUE, "term_initial" = c(0, 1), "verbose" = 0)
  strat_col <- c("b")
  options(warn = -1)
  expect_no_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
  keep_constant <- c(1, 1)
  expect_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
  keep_constant <- c(0, 0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 2, "guesses" = 10, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "strata" = TRUE, "term_initial" = c(0, 1), "verbose" = 0)
  expect_no_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
})
test_that("Poisson_basic_guess_cpp various_fixes", {
  fname <- "MULTI_COV.csv"
  colTypes <- c("double", "double", "integer", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  pyr <- "t1"
  event <- "lung"
  names <- c("a", "b")
  term_n <- c(0, 1)
  tform <- c("loglin", "loglin")
  strat_col <- "a"
  keep_constant <- c(0, 0)
  a_n <- c(0.01, -15)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 1, "guesses" = 1, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "strata" = FALSE, "term_initial" = c(0, 1), "verbose" = TRUE)
  expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
  model_control <- list("strata" = TRUE)
  expect_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  guesses_control <- list("iterations" = 1, "guesses" = 1, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "term_initial" = c(0, 1), "verbose" = TRUE)
  model_control <- list("strata" = FALSE)
  expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  keep_constant <- c(1, 1)
  expect_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  keep_constant <- c(0, 0)
  names <- c("a", "CONST")
  expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
})

test_that("Cox_tier_guess combinations", {
  fname <- "MULTI_COV.csv"
  colTypes <- c("double", "double", "integer", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("a", "b")
  term_n <- c(0, 1)
  tform <- c("loglin", "loglin")
  keep_constant <- c(0, 0)
  a_n <- c(0.01, -15)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 1, "guesses" = 1, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", strata = FALSE, term_initial = c(0))
  strat_col <- "a"
  keep_constant <- c(1, 1)
  options(warn = -1)
  expect_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
  names <- c("a", "CONST")
  keep_constant <- c(0, 0)
  guesses_control <- list("iterations" = 1, "guesses" = 1, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", strata = FALSE, term_initial = c(0), rmin = c(1, 1, 1, 1), rmax = c(1, 1))
  expect_no_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
  #
  keep_constant <- c(0, 0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 2, "guesses" = 10, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "strata" = TRUE, "term_initial" = c(0, 1), "verbose" = 0)
  expect_no_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col))
})

test_that("Cox_basic_guess_cpp combinations", {
  fname <- "MULTI_COV.csv"
  colTypes <- c("double", "double", "integer", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("a", "b")
  term_n <- c(0, 1)
  tform <- c("loglin", "loglin")
  keep_constant <- c(0, 0)
  a_n <- c(0.01, -15)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  strat_col <- "a"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  guesses_control <- list("iterations" = 2, "guesses" = 2, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "strata" = FALSE, "term_initial" = c(0, 1), "verbose" = TRUE)
  model_control <- list("strata" = TRUE)
  expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  guesses_control <- list("iterations" = 2, "guesses" = 2, "lin_min" = 0.001, "lin_max" = 1, "loglin_min" = -1, "loglin_max" = 1, "lin_method" = "uniform", "loglin_method" = "uniform", "term_initial" = c(0, 1), "verbose" = TRUE)
  names <- c("a", "CONST")
  expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  keep_constant <- c(1, 1)
  expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  keep_constant <- c(0, 0)
  model_control <- list("strata" = FALSE)
  lung_temp <- df$lung
  df$lung <- rep(0, nrow(df))
  expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, strat_col, model_control))
  model_control <- list("strata" = TRUE)
})
