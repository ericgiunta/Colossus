test_that("Coxph multidose", {
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
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  j_iterate <- 1
  LL_comp_1 <- c(-450.7215, -450.7215, -382.8276, -382.8276)
  LL_comp_2 <- c(-449.5319, -449.5319, -381.6798, -381.6798)
  LL_comp_3 <- c(-450.8742, -450.8742, -382.8966, -382.8966)
  k <- 1
  for (i in c(FALSE, TRUE)) {
    for (j in c(FALSE, TRUE)) {
      model_control <- list("strata" = i, "basic" = j)
      if (verbose) {
        print(model_control)
      }
      a_n <- c(-0.1, -0.1)
      # expect_equal(0,0)
      control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
      e <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")
      val <- e$LogLik
      expect_equal(LL_comp_1[k], val[1], tolerance = 1e-4)
      expect_equal(LL_comp_2[k], val[2], tolerance = 1e-4)
      expect_equal(LL_comp_3[k], val[3], tolerance = 1e-4)
      k <- k + 1
    }
  }
})
test_that("Coxph multidose negative shift check", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- -1 * df$rand1
  time1 <- "t0"
  time2 <- "t1"
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "plin")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(0, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)
  model_control <- list()
  # expect_equal(0,0)
  control <- list("ncores" = 2, "lr" = 0.95, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  e <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")
  val <- e$LogLik
  expect_equal(c(-450.7240, -449.4633, -449.4633), val, tolerance = 1e-4)
})
test_that("Coxph multidose, extra warnings and checks", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  #
  df$weighting <- floor(runif(nrow(df), min = 0, max = 2))
  time1 <- "t0"
  time2 <- "t1"
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  options(warn = -1)
  verbose <- FALSE
  j_iterate <- 1
  # LL_comp <- c(-69.51585, -69.51585, -77.97632, -77.97632, -59.95167, -60.05273, -75.34028, -75.3691, -69.51585, -69.51585, -77.97632, -77.97632, -59.95167, -60.05273, -75.34028, -75.3691, -111.3009, -111.3009, -119.9814, -119.9814, -100.8329, -101.007, -117.0147, -117.0539, -111.3009, -111.3009, -119.9814, -119.9814, -100.8329, -101.007, -117.0147, -117.0539)
  model_control <- list("cr" = TRUE)
  a_n <- c(-0.1, -0.1)
  # expect_equal(0,0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  expect_no_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "weighting"))
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "bad_weighting"))
  #
  keep_constant <- c(1, 1)
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "bad_weighting"))
  keep_constant <- c(0, 0)
  #
  names <- c("CONST", "rand")
  expect_no_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "weighting"))
  names <- c("dose", "rand")
  #
  model_control <- list("strata" = TRUE)
  d <- df[1, ]
  d$fac <- 101
  d$lung <- 0
  df <- rbind(df, d)
  expect_no_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "weighting"))
})
test_that("Coxph multidose failures", {
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
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  j_iterate <- 1

  model_control <- list("strata" = FALSE, "basic" = FALSE)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand", "more_rand")
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
  realization_columns <- matrix(c("rand0", "rand1", "bad"), nrow = 1)
  realization_index <- c("rand")
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
  realization_columns <- matrix(c("rand0", "rand1"), nrow = 1)
  realization_index <- c("rand")
  names <- c("bad", "rand")
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
  names <- c("dose", "rand")
  df$lung <- 0
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
})
test_that("Coxph multidose model failures", {
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
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  j_iterate <- 1

  model_control <- list("strata" = FALSE, "basic" = FALSE)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  model_control <- list("single" = TRUE, "basic" = FALSE)
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
  model_control <- list("null" = TRUE, "basic" = FALSE)
  expect_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
  model_control <- list("gradient" = TRUE, "basic" = FALSE)
  expect_no_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null"))
})
test_that("Coxph multidose MCML repeated column", {
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
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  e0 <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")

  realization_columns <- matrix(c("rand0", "rand0"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e1 <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")

  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
})
test_that("Coxph multidose MCML swapped columns", {
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
  # df$censor <- (df$lung==0)
  df$lung <- (df$lung > 0)
  # event <- "censor"
  names <- c("dose", "rand")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  realization_columns <- matrix(c("rand1", "rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  e0 <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")

  realization_columns <- matrix(c("rand0", "rand1"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e1 <- RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "fac", model_control = model_control, cens_weight = "null")

  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
})
