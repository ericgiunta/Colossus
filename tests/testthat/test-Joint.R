## ------------------------------------- ##
## Verify working results
## ------------------------------------- ##
test_that("Joint data generation, GMIX", {
  a <- c(0, 0, 0, 1, 1, 1)
  b <- c(1, 1, 1, 2, 2, 2)
  c <- c(0, 1, 2, 2, 1, 0)
  d <- c(1, 1, 0, 0, 1, 1)
  e <- c(0, 1, 1, 1, 0, 0)
  f <- c(1, 0, 0, 1, 1, 0)
  df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e, "res" = f)
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
  # Errors with one term
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(r)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(r)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(r)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  # Correct with multiple terms
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + loglin(res, 1) + gmix(e, e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_no_error(PoisRunJoint(formula_list, df, ncores = 1))
  # one of the term gmix values is different
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + loglin(res, 1) + gmix(e, e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(r)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + loglin(res, 1) + gmix(e, e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + GMIX(e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(r)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  #
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(r)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + loglin(res, 1) + gmix(e, e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + loglin(res, 1) + gmix(e, e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + gmix(r)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  #
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(r)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + gmix(e)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + loglin(res, 1) + gmix(e, e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
  model_1 <- Pois(pyr, e0) ~ loglin(fac, 0) + gmix(e)
  model_2 <- Pois(pyr, e1) ~ loglin(fac, 0) + gmix(r)
  model_s <- Pois(pyr) ~ plinear(t0, 0) + loglin(res, 1) + gmix(e, e)
  formula_list <- list(model_1, model_2, "shared" = model_s)
  expect_error(PoisRunJoint(formula_list, df, ncores = 1))
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
  expect_identical(model$term_n, rep(0L, 3))
  expect_identical(model$tform, c("plin", "loglin", "loglin"))
  expect_identical(model$keep_constant, rep(0, 3))
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
  expect_identical(model$names, c("t0", "fac_e0", "fac_e1"))
  expect_named(val$data, c("t0", "t1", "events", "e0", "e1", "fac", "pyr", "fac_e0", "fac_e1"))
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
  expect_true(e$Converged)
})

test_that("Data production checks and errors", {
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
  name_list <- list("e0" = "fac", "e1" = "fac", "shared" = "t0")
  term_n_list <- list("e0" = 0, "e1" = 0, "shared" = 0)
  tform_list <- list("e0" = "loglin", "e1" = "loglin", "shared" = "loglin")
  keep_constant_list <- list("e0" = 0, "e1" = 0, "shared" = 0)
  a_n_list <- list("e0" = 0.1, "e1" = 0.1, "shared" = 0.1)
  # Test it works
  expect_no_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list, keep_constant_list = keep_constant_list, a_n_list = a_n_list))
  # Check for mismatched number of items in term
  term_n_list <- list("e0" = c(0, 1), "e1" = 0, "shared" = 0)
  expect_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list))
  # Add missing term value
  term_n_list <- list("e1" = 0, "shared" = 0)
  expect_no_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list))
  term_n_list <- list("e0" = 0, "e1" = 0, "shared" = 0)
  # Check for mismatched number of items in subterm
  tform_list <- list("e0" = c("loglinear", "plinear"), "e1" = "loglinear", "shared" = "loglinear")
  expect_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list))
  # Add missing subterm value
  tform_list <- list("e1" = "loglin", "shared" = "loglin")
  expect_no_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list))
  tform_list <- list("e0" = "loglin", "e1" = "loglin", "shared" = "loglin")
  # Check for mismatched number of items in constant
  keep_constant_list <- list("e0" = c(0, 1), "e1" = 0, "shared" = 0)
  expect_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list, keep_constant_list = keep_constant_list, a_n_list = a_n_list))
  # Add missing constant value
  keep_constant_list <- list("e1" = 0, "shared" = 0)
  expect_no_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list, keep_constant_list = keep_constant_list, a_n_list = a_n_list))
  keep_constant_list <- list("e0" = 0, "e1" = 0, "shared" = 0)
  # Check for mismatched number of items in a_n
  a_n_list <- list("e0" = c(0.1, 0.1), "e1" = 0.1, "shared" = 0.1)
  expect_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list, keep_constant_list = keep_constant_list, a_n_list = a_n_list))
  # Add missing a_n value
  a_n_list <- list("e1" = 0.1, "shared" = 0.1)
  expect_no_error(Joint_Multiple_Events(df, events, name_list = name_list, term_n_list = term_n_list, tform_list = tform_list, keep_constant_list = keep_constant_list, a_n_list = a_n_list))
  a_n_list <- list("e0" = 0.1, "e1" = 0.1, "shared" = 0.1)
})
