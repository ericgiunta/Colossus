test_that("Coxph time column missing", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox(a_bad, b, c) ~ loglinear(d, 0) + multiplicative(), df, a_n = a_n, control = control))
})
test_that("Coxph no events", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  keep_constant <- c(0, 0)
  a_n <- c(-0.1, 0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox(a, b, c) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
})
#
test_that("Coxph_strata time column missing", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 0, 0, 1, 0, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox_Strata(a_bad, b, c, e) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
})
test_that("Coxph_strata no events", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 0, 0, 1, 0, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox_Strata(a, b, c, e) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
  #
})
test_that("Coxph_strata no events in strata", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 0, 1, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(2, 1, 2, 0, 1, 0, 2)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  keep_constant <- c(0, 0)
  a_n <- c(-0.1, 0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  options(warn = -1)
  expect_no_error(CoxRun(Cox_Strata(a, b, c, e) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
  options(warn = 0)
})
test_that("Coxph_strata no strata", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 0, 0, 1, 0, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox_Strata(a, b, c, e_bad) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
})
test_that("Coxph_strata strata with no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 0, 0, 1, 0, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 5, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(e0 <- CoxRun(Cox_Strata(a, b, c, e) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
  expect_no_error(e1 <- CoxRun(Cox_Strata(a, b, c, c(e)) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
  expect_no_error(e2 <- CoxRun(Cox_Strata(a, b, c, c(e, e)) ~ loglinear(d, 0) + loglinear(d, 1) + multiplicative(), df, a_n = a_n, control = control))
  #
})

test_that("Coxph null time column missing", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox(a_bad, b, c) ~ null(), df, a_n = a_n, control = control))
})
test_that("Coxph null no events", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(CoxRun(Cox(a, b, c) ~ null(), df, a_n = a_n, control = control))
})
#
test_that("Coxph dose list", {
  fname <- "dose.csv"
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  a_n <- c(1.0, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(t0, t1, lung) ~ loglin_dose(dose, 0) + lin_dose(dose, 0) + quad(dose, 0) + step_dose(dose, 0) + lin_quad_dose(dose, 0) + lin_exp_dose(dose, 0) + multiplicative(), df, a_n = a_n, keep_constant = keep_constant, control = control)
  expect_equal(e$beta_0, c(1.00000000, -0.10370531, -0.10093897, 0.98640244, -0.06027164, 1.12393188, 2.82997092, 0.17393188, 1.43935314, 0.32606812, 1.04573155, 1.12419452), tolerance = 1e-2)
})
#
test_that("Coxph fixed intercept", {
  fname <- "dose.csv"
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(0, 0, 1, 0, 1)
  a_n <- c(-0.1, 0.1, -1, 0.1, -1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, 0) + lin_dose(dose, 0) + step_dose(dose, 0) + multiplicative(), df, a_n = a_n, keep_constant = keep_constant, control = control))
})

test_that("Coxph loglin_M Strata", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(0)
  a_n <- c(0.01)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox_Strata(t0, t1, lung, fac) ~ loglinear(dose, 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$beta_0, c(-0.106), tolerance = 1e-2)
  e <- CoxRun(Cox_Strata(t0, t1, lung, c(fac)) ~ loglinear(dose, 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$beta_0, c(-0.106), tolerance = 1e-2)
  e <- CoxRun(Cox_Strata(t0, t1, lung, c(fac, fac)) ~ loglinear(dose, 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$beta_0, c(-0.106), tolerance = 1e-2)
})
test_that("Coxph loglin_M Single", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + multiplicative(), df, a_n = a_n, control = control, single = TRUE)
  expect_equal(e$AIC, 1056.299, tolerance = 1e-2)
})
test_that("Coxph loglin_M Null", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(t0, t1, lung) ~ null(), df, control = control)
  expect_equal(e$AIC, 1052, tolerance = 1e-2)
})
#
test_that("Coxph loglin_M CENSOR", {
  fname <- "ll_cens_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$beta_0, c(-1.19, 0.08), tolerance = 1e-2)
})

test_that("Coxph censoring weight", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  df$weighting <- 1 - df$t1 / 20
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)
  e0 <- CoxRun(FineGray(t0, t1, lung, weighting) ~ loglinear(dose, fac, 0) + multiplicative(), df, a_n = a_n, control = control)
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)
  e1 <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + multiplicative(), df, a_n = a_n, control = control)
  #
  expect_equal(e0$LogLik - e1$LogLik, -9.474311, tolerance = 1e-2)
  #
  expect_error(CoxRun(FineGray(t0, t1, lung, bad_weighting) ~ loglinear(dose, fac, 0) + multiplicative(), df, a_n = a_n, control = control))
})

test_that("dose nondose combinations", {
  fname <- "dose.csv"
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  #
  df$dose2 <- df$dose * df$dose
  df$a <- df$dose + 0.001
  df$b <- df$dose2 + 0.001
  #
  keep_constant <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  a_n <- c(1.0, -0.1, -0.1, 1, 0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, 1, 1, 1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  options(warn = -1)
  expect_no_error(CoxRun(Cox(t0, t1, lung) ~ loglin_dose(dose, 0) + lin_dose(dose, 0) + quad(dose, 0) + step_dose(dose, 0) + lin_quad_dose(dose, 0) + lin_exp_dose(dose, 0) + loglinear(b, 1) + linear(b, 1) + plinear(b, 1) + additive(), df, a_n = a_n, keep_constant = keep_constant, control = control))
  options(warn = 0)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  expect_no_error(CoxRun(Cox(t0, t1, lung) ~ loglin_dose(a, 0) + lin_dose(a, 0) + quad(a, 0) + step_dose(a, 0) + lin_quad_dose(a, 0) + lin_exp_dose(a, 0) + loglinear(b, 0) + linear(b, 0) + plinear(b, 0) + loglinear(b, 1) + linear(b, 1) + plinear(b, 1) + a(), df, a_n = a_n, keep_constant = keep_constant, control = control))
  expect_no_error(CoxRun(Cox(t0, t1, lung) ~ loglin_dose(a, 0) + lin_dose(a, 0) + quad(a, 0) + step_dose(a, 0) + lin_quad_dose(a, 0) + lin_exp_dose(a, 0) + loglinear(b, 0) + linear(b, 0) + plinear(b, 0) + loglinear(b, 1) + linear(b, 1) + plinear(b, 1) + pae(), df, a_n = a_n, keep_constant = keep_constant, control = control))
})
