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
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
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
  expect_error(EventAssignment(poisres, df, bad = "wrong"))
  expect_error(EventAssignment(poisres, df, assign_control = "wrong"))
})
test_that("Poisson Assigned Events, check results strata", {
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
  set.seed(3742)
  df$pyr <- df$Ending_Age - df$Starting_Age
  a_n <- c(-0.75, 0.1, -0.05)

  keep_constant <- c(0, 0, 0)

  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
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
  for (i in 1:9) {
    a_n <- 2 * runif(3) - 1
    poisres$beta_0 <- a_n
    e <- EventAssignment(poisres, df)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)

    poisres_strata$beta_0 <- a_n
    e <- EventAssignment(poisres_strata, df)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  }
})
#
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
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5,
    "epsilon" = 1e-3, "deriv_epsilon" = 1e-3,
    "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0,
    "verbose" = 0, "ties" = "breslow"
  )
  poisres <- PoisRun(Pois(pyr, Cancer_Status) ~ loglinear(a, 0) + loglinear(b, 1) + loglinear(c, 2), df, a_n = a_n, control = control)
  df$Cancer_Status <- rep(0, nrow(df))
  expect_error(EventAssignment(poisres, df))
})
#
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
    "ncores" = 2, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
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
    "ncores" = 2, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE,
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
#
#
# test_that("Poisson Assigned Events, Epicure Validated", {
#  df <- fread("Assigned_Check.csv")
#  col_list <- c("a", "b")
#  val <- factorize(df, col_list)
#  df <- val$df
#  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
#    skip("Cran Skip")
#  }
#  sum_values <- c(1.23715194448793, 0.060895797790186, 1.29804774227811, 1.67430308814, 0.0824135004185432, 1.75671658855854, 4.55122765982382, 0.224023120607428, 4.77525078043124, 12.3715194448793, 0.60895797790186, 12.9804774227811, 1.23715194448793, 0.478385445077363, 1.71553738956529, 1.67430308814, 0.478385445077363, 2.15268853321736, 4.55122765982382, 0.478385445077363, 5.02961310490118, 12.3715194448793, 0.478385445077363, 12.8499048899566, 1.23715194448793, 0.060895797790186, 1.29804774227811, 1.67430308814, 0.0824135004185432, 1.75671658855854, 4.55122765982382, 0.224023120607428, 4.77525078043124, 12.3715194448793, 0.60895797790186, 12.9804774227811, 1.23715194448793, -1.17625614669774, 0.060895797790186, 1.67430308814, -1.59188958772145, 0.0824135004185432, 4.55122765982382, -4.32720453921639, 0.224023120607428, 12.3715194448793, -11.7625614669774, 0.60895797790186)
#  min_values <- c(0.1, 0.000230110913103737, 0.103747930428476, 0.135335283236613, 0.000311421256007305, 0.140407555497502, 0.367879441171442, 0.000846530741200557, 0.381667306687213, 1, 0.00230110913103743, 1.03747930428476, 0.1, 0.00170470357204963, 0.13669058432965, 0.135335283236613, 0.00170470357204963, 0.17281458752137, 0.367879441171442, 0.0017047035720496, 0.405358745456199, 1, 0.00170470357204966, 1.03747930428476, 0.1, 0.000230110913103737, 0.103747930428476, 0.135335283236613, 0.000311421256007305, 0.140407555497502, 0.367879441171442, 0.000846530741200557, 0.381667306687213, 1, 0.00230110913103743, 1.03747930428476, 0.1, -0.14873408300218, 0.000230110913103737, 0.135335283236613, -0.201289692500379, 0.000311421256007305, 0.367879441171442, -0.547162113379888, 0.000846530741200557, 1, -1.4873408300218, 0.00230110913103743)
#  max_values <- c(0.164872127070013, 0.016138044067833, 0.181010171137846, 0.22313016014843, 0.0218404676480512, 0.244970627796481, 0.606530659712633, 0.0593685463327452, 0.665899206045379, 1.64872127070013, 0.16138044067833, 1.81010171137846, 0.164872127070013, 0.0978821851493432, 0.262754312219356, 0.22313016014843, 0.0978821851493432, 0.321012345297773, 0.606530659712633, 0.0978821851493432, 0.704412844861977, 1.64872127070013, 0.0978821851493432, 1.74660345584947, 0.164872127070013, 0.016138044067833, 0.181010171137846, 0.22313016014843, 0.0218404676480512, 0.244970627796481, 0.606530659712633, 0.0593685463327452, 0.665899206045379, 1.64872127070013, 0.16138044067833, 1.81010171137846, 0.164872127070013, -0.0960313868406136, 0.016138044067833, 0.22313016014843, -0.129964349376792, 0.0218404676480512, 0.606530659712633, -0.353279729258435, 0.0593685463327452, 1.64872127070013, -0.960313868406136, 0.16138044067833)
#  event <- "e"
#  t0 <- "t0"
#  t1 <- "t1"
#  pyr <- "pyr"
#  names <- c("a", "b", "c", "CONST")
#  j <- 1
#    for (para in c(-2.3025850929940455, -2, -1, 0.0)) {
#      a_n <- c(0.1, 0.15, 0.02, para)
#      keep_constant <- rep(0, length(names))
#      keep_constant[4] <- 1
#
#      control <- list("ncores" = 2, "maxiters" = c(1, 1), "verbose" = FALSE, "step_max" = 0.1)
#      poisres <- PoisRun(Pois(pyr, e) ~ loglinear(a,b,c, 0) + linear(CONST, 1) + multiplicative(), df, a_n = a_n, control = control, keep_constant = keep_constant)
#      e <- EventAssignment(poisres, df)
#      e <- e$predict
#      BK <- e[, 1]
#      EX <- e[, 2]
#      FV <- e[, 3]
#
#      # BK values
#      expect_equal(sum(BK), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(BK), min_values[j], tolerance = 1e-2)
#      expect_equal(max(BK), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # EX values
#      expect_equal(sum(EX), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(EX), min_values[j], tolerance = 1e-2)
#      expect_equal(max(EX), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # FV values
#      expect_equal(sum(FV), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(FV), min_values[j], tolerance = 1e-2)
#      expect_equal(max(FV), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#    }
#    for (para in c(-2.3025850929940455, -2, -1, 0.0)) {
#      a_n <- c(0.1, 0.15, 0.02, para)
#      keep_constant <- rep(0, length(names))
#      keep_constant[4] <- 1
#
#      control <- list("ncores" = 2, "maxiters" = c(1, 1), "verbose" = FALSE, "step_max" = 0.1)
#      poisres <- PoisRun(Pois(pyr, e) ~ loglinear(a,b,c, 0) + linear(CONST, 1) + additive(), df, a_n = a_n, control = control, keep_constant = keep_constant)
#      e <- EventAssignment(poisres, df)
#      e <- e$predict
#      BK <- e[, 1]
#      EX <- e[, 2]
#      FV <- e[, 3]
#
#      # BK values
#      expect_equal(sum(BK), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(BK), min_values[j], tolerance = 1e-2)
#      expect_equal(max(BK), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # EX values
#      expect_equal(sum(EX), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(EX), min_values[j], tolerance = 1e-2)
#      expect_equal(max(EX), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # FV values
#      expect_equal(sum(FV), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(FV), min_values[j], tolerance = 1e-2)
#      expect_equal(max(FV), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#    }
#    for (para in c(-2.3025850929940455, -2, -1, 0.0)) {
#      a_n <- c(0.1, 0.15, 0.02, para)
#      keep_constant <- rep(0, length(names))
#      keep_constant[4] <- 1
#
#      control <- list("ncores" = 2, "maxiters" = c(1, 1), "verbose" = FALSE, "step_max" = 0.1)
#      poisres <- PoisRun(Pois(pyr, e) ~ loglinear(a,b,c, 0) + linear(CONST, 1) + pae(), df, a_n = a_n, control = control, keep_constant = keep_constant)
#      e <- EventAssignment(poisres, df)
#      e <- e$predict
#      BK <- e[, 1]
#      EX <- e[, 2]
#      FV <- e[, 3]
#
#      # BK values
#      expect_equal(sum(BK), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(BK), min_values[j], tolerance = 1e-2)
#      expect_equal(max(BK), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # EX values
#      expect_equal(sum(EX), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(EX), min_values[j], tolerance = 1e-2)
#      expect_equal(max(EX), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # FV values
#      expect_equal(sum(FV), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(FV), min_values[j], tolerance = 1e-2)
#      expect_equal(max(FV), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#    }
#    for (para in c(-2.3025850929940455, -2, -1, 0.0)) {
#      a_n <- c(0.1, 0.15, 0.02, para)
#      keep_constant <- rep(0, length(names))
#      keep_constant[4] <- 1
#
#      control <- list("ncores" = 2, "maxiters" = c(1, 1), "verbose" = FALSE, "step_max" = 0.1)
#      poisres <- PoisRun(Pois(pyr, e) ~ loglinear(a,b,c, 0) + linear(CONST, 1) + pa(), df, a_n = a_n, control = control, keep_constant = keep_constant)
#      e <- EventAssignment(poisres, df)
#      e <- e$predict
#      BK <- e[, 1]
#      EX <- e[, 2]
#      FV <- e[, 3]
#
#      # BK values
#      expect_equal(sum(BK), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(BK), min_values[j], tolerance = 1e-2)
#      expect_equal(max(BK), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # EX values
#      expect_equal(sum(EX), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(EX), min_values[j], tolerance = 1e-2)
#      expect_equal(max(EX), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#
#      # FV values
#      expect_equal(sum(FV), sum_values[j], tolerance = 1e-2)
#      expect_equal(min(FV), min_values[j], tolerance = 1e-2)
#      expect_equal(max(FV), max_values[j], tolerance = 1e-2)
#      j <- j + 1
#    }
# })
