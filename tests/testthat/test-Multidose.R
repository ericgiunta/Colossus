test_that("Coxph multidose", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  j_iterate <- 1
  LL_comp_1 <- c(-450.7215, -382.8276)
  LL_comp_2 <- c(-449.5319, -381.6798)
  LL_comp_3 <- c(-450.8742, -382.8966)
  k <- 1
  a_n <- c(-0.1, -0.1)
  # expect_equal(0,0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  ##
  expect_error(CoxRunMulti("bad", df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, bad = "wrong"))
  expect_no_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, fma = TRUE, ncores = 1))
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = "wrong"))
  ##
  e <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  #
  val <- e$LogLik
  expect_equal(LL_comp_1[k], val[1], tolerance = 1e-4)
  expect_equal(LL_comp_2[k], val[2], tolerance = 1e-4)
  expect_equal(LL_comp_3[k], val[3], tolerance = 1e-4)
  k <- k + 1
  e <- CoxRunMulti(Cox_Strata(t0, t1, lung, fac) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  val <- e$LogLik
  expect_equal(LL_comp_1[k], val[1], tolerance = 1e-4)
  expect_equal(LL_comp_2[k], val[2], tolerance = 1e-4)
  expect_equal(LL_comp_3[k], val[3], tolerance = 1e-4)
})
test_that("Pois multidose", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  j_iterate <- 1
  LL_comp_1 <- c(-329.7453, -397.8750)
  LL_comp_2 <- c(-328.7705, -387.2309)
  LL_comp_3 <- c(-329.9031, -405.6820)
  k <- 1
  a_n <- c(-0.02, -0.1, -0.1)
  # expect_equal(0,0)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model <- Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0)
  res <- get_form(model, df)
  expect_no_error(PoisRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, fma = TRUE))
  expect_no_error(PoisRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, ncores = 1))
  expect_error(PoisRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = "bad"))
  expect_error(PoisRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, bad = "wrong"))
  expect_error(PoisRunMulti("bad", df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  e <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  val <- e$LogLik
  expect_equal(LL_comp_1[k], val[1], tolerance = 1e-4)
  expect_equal(LL_comp_2[k], val[2], tolerance = 1e-4)
  expect_equal(LL_comp_3[k], val[3], tolerance = 1e-4)
  k <- k + 1
  e <- PoisRunMulti(Pois_Strata(t1, lung, fac) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  val <- e$LogLik
  expect_equal(LL_comp_1[k], val[1], tolerance = 1e-4)
  expect_equal(LL_comp_2[k], val[2], tolerance = 1e-4)
  expect_equal(LL_comp_3[k], val[3], tolerance = 1e-4)
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(0, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)
  model_control <- list()
  # expect_equal(0,0)
  control <- list("ncores" = 1, "lr" = 0.95, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, 0) + plinear(rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  val <- e$LogLik
  expect_equal(c(-450.7240, -449.4633, -449.4633), val, tolerance = 1e-4)
})
test_that("Pois multidose negative shift check", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- -1 * df$rand1
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(0, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.02, -0.1, -0.1)
  keep_constant <- c(0, 0, 0)
  model_control <- list()
  # expect_equal(0,0)
  control <- list("ncores" = 1, "lr" = 0.95, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, 0) + plinear(rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control)
  val <- e$LogLik
  expect_equal(c(-329.7473, -328.71485, -328.71485), val, tolerance = 1e-3)
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  options(warn = -1)
  verbose <- FALSE
  j_iterate <- 1
  model_control <- list("cr" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model <- FineGray(t0, t1, lung, weighting) ~ loglinear(dose, rand, 0)
  res <- get_form(model, df)
  expect_no_error(CoxRunMulti(res$model, df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, ncores = 1))
  expect_error(CoxRunMulti(FineGray(t0, t1, lung, bad_weighting) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  #
  keep_constant <- c(1, 1)
  expect_error(CoxRunMulti(FineGray(t0, t1, lung, weighting) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  keep_constant <- c(0, 0)
  #
  names <- c("CONST", "rand")
  expect_no_error(CoxRunMulti(FineGray(t0, t1, lung, weighting) ~ loglinear(CONST, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  names <- c("dose", "rand")
  #
  model_control <- list("strata" = TRUE)
  d <- df[1, ]
  d$fac <- 101
  d$lung <- 0
  df <- rbind(df, d)
  expect_no_error(e <- CoxRunMulti(Cox_Strata(t0, t1, lung, fac) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  options(warn = 0)
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  j_iterate <- 1

  model_control <- list("strata" = FALSE, "basic" = FALSE)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand", "more_rand")
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  realization_columns <- matrix(c("rand0", "rand1", "bad"), nrow = 1)
  realization_index <- c("rand")
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  realization_columns <- matrix(c("rand0", "rand1"), nrow = 1)
  realization_index <- c("rand")
  names <- c("bad", "rand")
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(bad, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  names <- c("dose", "rand")
  df$lung <- 0
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  j_iterate <- 1

  model_control <- list("strata" = FALSE, "basic" = FALSE)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model_control <- list("single" = TRUE, "basic" = FALSE)
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, single = TRUE))
  model_control <- list("null" = TRUE, "basic" = FALSE)
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ null(), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control))
  model_control <- list("gradient" = TRUE, "basic" = FALSE)
  expect_no_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list()))
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e0 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand0"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e1 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)

  realization_columns <- matrix(c("rand0"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e0 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand0"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e1 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, mcml = TRUE)
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
})
test_that("Poisson multidose MCML repeated column", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.02, -0.1, -0.1)
  keep_constant <- c(0, 0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.02, -0.1, -0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e0 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand0"), nrow = 1)
  a_n <- c(-0.02, -0.1, -0.1)
  e1 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
  #
  realization_columns <- matrix(c("rand0"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e0 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand0"), nrow = 1)
  a_n <- c(-0.02, -0.1, -0.1)
  e1 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, mcml = TRUE)
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
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand1", "rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e0 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand1"), nrow = 1)
  a_n <- c(-0.1, -0.1)
  e1 <- CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
  expect_error(CoxRunMulti(Cox(t0, t1, lung) ~ loglinear(dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE, fma = TRUE))
})
test_that("Pois multidose MCML swapped columns", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
  df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
  df$lung <- (df$lung > 0)
  realization_columns <- matrix(c("rand1", "rand0"), nrow = 1)
  realization_index <- c("rand")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  cens_weight <- c(0)
  #
  event <- "lung"
  a_n <- c(-0.02, -0.1, -0.1)
  keep_constant <- c(0, 0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  model_control <- list("mcml" = TRUE)
  a_n <- c(-0.02, -0.1, -0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  realization_columns <- matrix(c("rand1", "rand0"), nrow = 1)
  e0 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  realization_columns <- matrix(c("rand0", "rand1"), nrow = 1)
  a_n <- c(-0.02, -0.1, -0.1)
  e1 <- PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE)
  expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-4)
  expect_error(PoisRunMulti(Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0), df, a_n = a_n, keep_constant = keep_constant, realization_columns = realization_columns, realization_index = realization_index, control = control, gradient_control = list(), mcml = TRUE, fma = TRUE))
})
