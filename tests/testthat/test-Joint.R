## ------------------------------------- ##
## Verify working results
## ------------------------------------- ##
test_that("Joint data generation, no error", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
  df$pyr <- df$t1 - df$t0
  events <- c("e0", "e1")
  keep_constant_shared <- c(0, 0)
  a_n_shared <- c(0.001, -0.02)
  #
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix()
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX()
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix()
  formula_list <- list(model_1, model_2, "shared" = model_s)
  #
  expect_no_error(PoisRunJoint(formula_list, df, ncores = 1))
  expect_no_error(PoisRunJoint(formula_list, df, ncores = 1, norm = "mean"))
  expect_no_error(PoisRunJoint(formula_list, df, ncores = 1, norm = "max"))
  expect_error(PoisRunJoint(formula_list, df, ncores = 1, norm = "bad"))
  expect_error(PoisRunJoint(formula_list, df, ncores = 1, bad = "wrong"))
})
test_that("Joint data generation fill defaults, no error", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
  df$pyr <- df$t1 - df$t0
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0)
  model_s <- Pois(pyr) ~ plinear(t0, 0)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  #
  expect_no_error(get_form_joint(formula_list, df))
  val <- get_form_joint(formula_list, df)
  model <- val$model
  expect_equal(model$term_n, rep(0, 3))
  expect_equal(model$tform, c("plin", "loglin", "loglin"))
  expect_equal(model$keep_constant, rep(0, 3))
})
test_that("Joint data generation, check results", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
  df$pyr <- df$t1 - df$t0
  events <- c("e0", "e1")
  #
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0)
  model_s <- Pois(pyr) ~ plinear(t0, 0)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  #
  val <- get_form_joint(formula_list, df)
  model <- val$model
  expect_equal(model$names, c("t0", "fac_e0", "fac_e1"))
  expect_equal(names(val$data), c("t0", "t1", "events", "e0", "e1", "fac", "pyr", "fac_e0", "fac_e1"))
})
test_that("Joint data regression, no error", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  f <- c(0, 1, 0, 0, 1, 1)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
  df$pyr <- df$t1 - df$t0
  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 2, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow"
  )
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0)
  model_s <- Pois(pyr) ~ plinear(t0, 0)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_no_error(PoisRunJoint(formula_list, df, control = control))
  formula_list <- list(model_1, model_2)
  expect_no_error(PoisRunJoint(formula_list, df, control = control))
})
test_that("Joint data regression, check results", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  f <- c(0, 1, 0, 0, 1, 1)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
  df$pyr <- df$t1 - df$t0
  events <- c("e0", "e1")
  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 10, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow"
  )
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0)
  model_s <- Pois(pyr) ~ plinear(t0, 0)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  e <- PoisRunJoint(formula_list, df, control = control)
  expect_equal(e$beta_0, c(-0.1845, 0.5742, -1.0347), tolerance = 1e-2)
  expect_equal(e$Converged, TRUE)
})
