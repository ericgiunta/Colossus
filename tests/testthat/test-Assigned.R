## ------------------------------------- ##
## Verify working results
## ------------------------------------- ##
test_that("Poisson Assigned Events, check results", {
  df <- data.table::data.table(
    "UserID" = c(112, 114, 213, 214, 115, 116, 117),
    "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
    "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
    "Cancer_Status" = c(12, 10, 18, 6, 1, 11, 4),
    "a" = c(0, 1, 1, 0, 1, 0, 1),
    "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
    "c" = c(10, 11, 10, 11, 12, 9, 11),
    "d" = c(0, 0, 0, 1, 1, 1, 1)
  )

  df$pyr <- df$Ending_Age - df$Starting_Age
  a_n <- c(-0.75, 0.1, -0.05, -1.5)

  keep_constant <- c(0, 0, 0, 0)

  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0
  )
  #
  poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + linear(b, c, 1) + plinear(d, 2), df, a_n = a_n, control = control)
  e <- EventAssignment(poisres, df)

  e0 <- e$predict
  e1 <- e$caused
  expect_equal(e1[, 3], df$Cancer_Status, tolerance = 1e-1)
  expect_equal(sum(e0), 108.757, tolerance = 1e-2)
  expect_equal(sum(e1), 124, tolerance = 1e-2)

  expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
  expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  #
  expect_no_error(df$bk <- e0[, 1])
  #
  expect_error(EventAssignment(poisres, df, bad = "wrong"))
  expect_error(EventAssignment(poisres, df, assign_control = "wrong"))
})
test_that("Poisson Assigned Events, check results strata", {
  df <- data.table::data.table(
    "UserID" = c(112, 114, 213, 214, 115, 116, 117, 118),
    "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18, 56),
    "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55, 57),
    "Cancer_Status" = c(12, 10, 18, 6, 1, 11, 4, 0),
    "a" = c(0, 1, 1, 0, 1, 0, 1, 1),
    "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2, 1),
    "c" = c(10, 11, 10, 11, 12, 9, 11, 5),
    "d" = c(0, 0, 0, 1, 1, 1, 1, 2)
  )
  set.seed(3742)
  df$pyr <- df$Ending_Age - df$Starting_Age
  a_n <- c(-0.75, 0.1, -0.05)

  keep_constant <- c(0, 0, 0)

  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 0.2, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0
  )
  #
  term_n <- c(0, 1, 2)
  model_control <- list("strata" = FALSE)
  poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + loglinear(b, 1) + loglinear(c, 2), df, a_n = a_n, control = control)
  poisres_strata <- PoisRun(Pois_Strata(pyr, Cancer_Status, d) ~ loglinear(a, 0) + loglinear(b, 1) + loglinear(c, 2), df, a_n = a_n, control = control)
  expect_no_error(EventAssignment(poisres, df))
  expect_no_error(EventAssignment(poisres_strata, df))
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 0.2, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0
  )
  res_0 <- c(166.1197, 185.6144, 310.4459, 152.3147, 124.5314, 211.6529, 287.7323, 134.2735, 169.0388)
  res_1 <- c(29.293795245, 0.004552392, 16.705778339, 25.508603705, 0.018238072, 28.813274218, 21.602469009, 2.007582573, 0.131742868)
  res_2 <- c(29.346103747, 0.002581566, 15.799809687, 25.554827359, 0.013333945, 28.864420253, 21.396141586, 1.838636616, 0.100627725)
  res_3 <- c(29.293795245, 0.004552392, 16.705778339, 25.508603705, 0.018238072, 28.813274218, 21.602469009, 2.007582573, 0.131742868)
  for (i in 1:9) {
    a_n <- 2 * runif(3) - 1
    poisres$beta_0 <- a_n
    e <- EventAssignment(poisres, df, control = control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1]), res_0[i], tolerance = 1e-2)
    expect_equal(sum(e1[, 1]), res_1[i], tolerance = 1e-2)
    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)

    poisres_strata$beta_0 <- a_n
    e <- EventAssignment(poisres_strata, df, control = control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1]), res_2[i], tolerance = 1e-2)
    expect_equal(sum(e1[, 1]), res_3[i], tolerance = 1e-2)
    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  }
})

test_that("Poisson Assigned Events, combinations", {
  df <- data.table::data.table(
    "UserID" = c(112, 114, 213, 214, 115, 116, 117),
    "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
    "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
    "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
    "a" = c(0, 1, 1, 0, 1, 0, 1),
    "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
    "c" = c(10, 11, 10, 11, 12, 9, 11),
    "d" = c(0, 0, 0, 1, 1, 1, 1)
  )
  # For the interval case
  df$pyr <- df$Ending_Age - df$Starting_Age
  pyr <- "pyr"
  a_n <- c(0.1, 0.1, 0.1)

  keep_constant <- c(0, 0, 0, 0)
  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5,
    "epsilon" = 1e-3, "deriv_epsilon" = 1e-3,
    "step_max" = 0.2, "change_all" = TRUE, "thres_step_max" = 100.0,
    "verbose" = 0, "ties" = "breslow"
  )
  poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + loglinear(b, 1) + loglinear(c, 2), df, a_n = a_n, control = control)
  df$Cancer_Status <- rep(0, nrow(df))
  expect_error(EventAssignment(poisres, df))
})

test_that("Poisson Assigned Events bounds, check results", {
  df <- data.table::data.table(
    "UserID" = c(112, 114, 213, 214, 115, 116, 117),
    "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
    "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
    "Cancer_Status" = c(12, 10, 18, 6, 1, 11, 4),
    "a" = c(0, 1, 1, 0, 1, 0, 1),
    "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
    "c" = c(10, 11, 10, 11, 12, 9, 11),
    "d" = c(0, 0, 0, 1, 1, 1, 1)
  )

  df$pyr <- df$Ending_Age - df$Starting_Age
  a_n <- c(-0.75, 0.1, -0.05, -1.5)
  keep_constant <- c(0, 0, 0, 0)
  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 0.2, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0
  )
  #
  poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + linear(b, c, 1) + plinear(d, 2), df, a_n = a_n, control = control, norm = "max")
  assign_control <- list(check_num = 4)
  e <- EventAssignment(poisres, df, assign_control = assign_control, z = 2)

  elow <- e$lower_limit$predict
  emid <- e$midpoint$predict
  eupp <- e$upper_limit$predict
  #
  expect_equal(sum(elow), 123.4787, tolerance = 1e-2)
  expect_equal(sum(emid), 123.5965, tolerance = 1e-2)
  expect_equal(sum(eupp), 124.0949, tolerance = 1e-2)
  #
  poisbound <- LikelihoodBound(poisres, df, para_number = 4, maxstep = 50)
  assign_control <- list(check_num = 4)
  e <- EventAssignment(poisbound, df, assign_control = assign_control, z = 2)
  #
  elow <- e$lower_limit$predict
  emid <- e$midpoint$predict
  eupp <- e$upper_limit$predict
  #
  expect_equal(sum(elow), 123.4753, tolerance = 1e-2)
  expect_equal(sum(emid), 123.5965, tolerance = 1e-2)
  expect_equal(sum(eupp), 124.5127, tolerance = 1e-2)
})
test_that("Poisson Assigned Events bounds single entry, check results", {
  df <- data.table::data.table(
    "UserID" = c(112, 114, 213, 214, 115, 116, 117),
    "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
    "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
    "Cancer_Status" = c(12, 10, 18, 6, 1, 11, 4),
    "a" = c(0, 1, 1, 0, 1, 0, 1),
    "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
    "c" = c(10, 11, 10, 11, 12, 9, 11),
    "d" = c(0, 0, 0, 1, 1, 1, 1)
  )

  df$pyr <- df$Ending_Age - df$Starting_Age
  a_n <- c(-0.75)
  keep_constant <- c(0)
  control <- list(
    "ncores" = 1, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 0.2, "change_all" = TRUE,
    "thres_step_max" = 100.0, "verbose" = 0
  )
  #
  model <- Pois(pyr, Cancer_Status) ~ plinear(a, 0)
  poisres <- PoisRun(model, df, a_n = a_n, control = control, norm = "max")
  assign_control <- list(check_num = 1)
  e <- EventAssignment(poisres, df, assign_control = assign_control, z = 2)

  elow <- e$lower_limit$predict
  emid <- e$midpoint$predict
  eupp <- e$upper_limit$predict
  #
  expect_equal(sum(elow), 202.9804, tolerance = 1e-2)
  expect_equal(sum(emid), 225.95, tolerance = 1e-2)
  expect_equal(sum(eupp), 248.9195, tolerance = 1e-2)
  #
  poisbound <- LikelihoodBound(poisres, df, para_number = 1, maxstep = 50)
  assign_control <- list(check_num = 1)
  e <- EventAssignment(poisbound, df, assign_control = assign_control, z = 2)
  #
  elow <- e$lower_limit$predict
  emid <- e$midpoint$predict
  eupp <- e$upper_limit$predict
  #
  expect_equal(sum(elow), 205.9671, tolerance = 1e-2)
  expect_equal(sum(emid), 225.95, tolerance = 1e-2)
  expect_equal(sum(eupp), 251.1475, tolerance = 1e-2)
})
