test_that("Logistic multioutcome", {
  # generating data
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  realization_columns <- c("event0", "event1", "event2")
  #
  a_n <- c(0, 0)
  keep_constant <- c(0, 0)
  #
  # testing for basic run errors and warnings
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model <- logit(event0) ~ loglinear(CONST, dose, 0)
  res <- get_form(model, df)
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, ncores = 1))
  expect_warning(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, ncores = 1, ncores = 1))
  expect_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "bad"))
  expect_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = "bad"))
  expect_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, bad = "wrong"))
  expect_error(LogisticRunMulti("bad", df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  e <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control)
  val <- e$LogLik ## -628.2178 -607.3623 -608.5099
  #
  # testing against general logistic
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2)
  e <- LogisticRun(res$model, df, a_n = a_n, keep_constant = keep_constant, control = control)
  ll_comp_1 <- e$LogLik ## -628.2261
  model <- logit(event1) ~ loglinear(CONST, dose, 0)
  res <- get_form(model, df)
  e <- LogisticRun(res$model, df, a_n = a_n, keep_constant = keep_constant, control = control)
  ll_comp_2 <- e$LogLik ## -607.3674
  model <- logit(event2) ~ loglinear(CONST, dose, 0)
  res <- get_form(model, df)
  e <- LogisticRun(res$model, df, a_n = a_n, keep_constant = keep_constant, control = control)
  ll_comp_3 <- e$LogLik ## -608.5099
  expect_equal(val[1], ll_comp_1, tolerance = 1e-4)
  expect_equal(val[2], ll_comp_2, tolerance = 1e-4)
  expect_equal(val[3], ll_comp_3, tolerance = 1e-4)
})
test_that("Logistic multioutcome extra checks", {
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  realization_columns <- c("event0", "event1", "event2")
  #
  a_n <- c(0, 0)
  keep_constant <- c(0, 0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model <- logit(event0) ~ loglinear(CONST, dose, 0)
  res <- get_form(model, df)
  #
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "odds", observed_info = TRUE))
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "ident", observed_info = TRUE))
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "loglink", observed_info = TRUE))
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "probit", observed_info = TRUE))
  #
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = c(control, constraint = TRUE)))
  expect_no_error(LogisticRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = c(control, constraint = TRUE, gradient = TRUE)))
  #
  expect_error(LogisticRunMulti(logit(event0) ~ linear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, link = "ident"))
})
test_that("Logistic multioutcome failures", {
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  realization_columns <- "event0"
  #
  a_n <- c(0, 0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  keep_constant <- c(1, 1)
  expect_error(LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  keep_constant <- c(0, 0)
  #
  expect_error(LogisticRunMulti(logit(event0) ~ null(), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  #
  df$event0[1] <- -1
  expect_error(LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  df$event0 <- rep(0, times = nrow(df))
  expect_error(LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  #
  realization_columns <- "event3"
  expect_error(LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
  df$event3 <- rep("NAN", times = nrow(df))
  expect_error(LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control))
})
test_that("Logistic multioutcome repeated column", {
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  #
  a_n <- c(-0.1, -0.2)
  keep_constant <- c(1, 0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  # with gradient control
  realization_columns <- "event0"
  e0 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, gradient_control = list())
  realization_columns <- c("event0", "event0")
  e1 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, gradient_control = list())
  expect_equal(rep(e0$LogLik, 2), e1$LogLik, tolerance = 1e-4) ## -645.6296
  #
  # without gradient control
  realization_columns <- "event0"
  e0 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control)
  realization_columns <- c("event0", "event0")
  e1 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control)
  expect_equal(rep(e0$LogLik, 2), e1$LogLik, tolerance = 1e-4) ## -644.7275
})
test_that("Logistic multioutcome swapped columns", {
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  #
  a_n <- c(0, 0)
  keep_constant <- c(1, 0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  realization_columns <- c("event0", "event1")
  e0 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, gradient_control = list())
  realization_columns <- c("event1", "event0")
  e1 <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, gradient_control = list())
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4) ## -693.1472
})
test_that("Logistic multioutcome model, single check", {
  fname <- "ll_comp_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, select = 1:3, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$event0 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event1 <- rbinom(nrow(df), size = 1, prob = 0.3)
  df$event2 <- rbinom(nrow(df), size = 1, prob = 0.3)
  realization_columns <- c("event0", "event1", "event2")
  #
  a_n <- c(-0.1, 0)
  keep_constant <- c(0, 0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 1, "halfmax" = 1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 2, "ties" = "breslow")
  #
  e <- LogisticRunMulti(logit(event0) ~ loglinear(CONST, dose, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, control = control, single = TRUE, observed_info = TRUE)
  expect_equal(e$LogLik, c(-676.5967, -673.9967, -674.1967), tolerance = 1e-4)
})
