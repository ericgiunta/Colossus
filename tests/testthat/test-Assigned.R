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
  pyr <- "pyr"
  event <- "Cancer_Status"
  names <- c("a", "b", "c", "d")
  term_n <- c(0, 1, 1, 2)
  tform <- c("loglin", "lin", "lin", "plin")
  modelform <- "M"
  a_n <- c(-0.75, 0.1, -0.05, -1.5)

  keep_constant <- c(0, 0, 0, 0)

  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "abs_max" = 1.0, "change_all" = TRUE,
    "dose_abs_max" = 100.0, "verbose" = 0, "double_step" = 1
  )
  #
  e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control)

  e0 <- e$predict
  e1 <- e$caused

  expect_equal(sum(e0), 162.8914, tolerance = 1e-2)
  expect_equal(sum(e1), 124, tolerance = 1e-2)

  expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
  expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
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
  pyr <- "pyr"
  event <- "Cancer_Status"
  names <- c("a", "b", "c")
  term_n <- c(0, 1, 2)
  tform <- c("loglin", "loglin", "loglin")
  modelform <- "M"
  a_n <- c(-0.75, 0.1, -0.05)

  keep_constant <- c(0, 0, 0)

  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "abs_max" = 1.0, "change_all" = TRUE,
    "dose_abs_max" = 100.0, "verbose" = 0, "double_step" = 1
  )
  #
  term_n <- c(0, 1, 2)
  model_control <- list("strata" = FALSE)
  expect_no_error(RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control))
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (i in 1:3) {
    a_n <- 2 * runif(3) - 1
    model_control <- list("strata" = FALSE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
    model_control <- list("strata" = TRUE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "d", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  }
  term_n <- c(0, 2, 1)
  for (i in 1:3) {
    a_n <- 2 * runif(3) - 1
    model_control <- list("strata" = FALSE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
    model_control <- list("strata" = TRUE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "d", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  }
  term_n <- c(1, 0, 2)
  for (i in 1:3) {
    a_n <- 2 * runif(3) - 1
    model_control <- list("strata" = FALSE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
    model_control <- list("strata" = TRUE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "d", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
  }
  term_n <- c(2, 0, 1)
  for (i in 1:3) {
    a_n <- 2 * runif(3) - 1
    model_control <- list("strata" = FALSE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control)

    e0 <- e$predict
    e1 <- e$caused

    expect_equal(sum(e0[, 1:2]), sum(e0[, 3]), tolerance = 1e-2)
    expect_equal(sum(e1[, 1:2]), sum(e1[, 3]), tolerance = 1e-2)
    model_control <- list("strata" = TRUE)
    e <- RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "d", model_control)

    e0 <- e$predict
    e1 <- e$caused

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
  time1 <- "Starting_Age"
  time2 <- "Ending_Age"
  df$pyr <- df$Ending_Age - df$Starting_Age
  pyr <- "pyr"
  event <- "Cancer_Status"
  names <- c("a", "b", "c", "d")
  term_n <- c(0, 1, 1, 2)
  tform <- c("loglin", "lin", "lin", "plin")
  modelform <- "M"
  a_n <- c(0.1, 0.1, 0.1, 0.1)

  keep_constant <- c(0, 0, 0, 0)
  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 1, "halfmax" = 5,
    "epsilon" = 1e-3, "deriv_epsilon" = 1e-3,
    "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0,
    "verbose" = 0, "ties" = "breslow", "double_step" = 1
  )
  model_control <- list("strata" = TRUE)
  expect_error(RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, "null", model_control))
  keep_constant <- c(1, 1, 1, 1)
  expect_error(RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control))
  names <- c("a", "b", "CONST", "d")
  keep_constant <- c(0, 0, 0, 0)
  expect_no_error(RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control))

  df$Cancer_Status <- rep(0, nrow(df))
  expect_error(RunPoissonEventAssignment(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control))
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
  pyr <- "pyr"
  event <- "Cancer_Status"
  names <- c("a", "b", "c", "d")
  term_n <- c(0, 1, 1, 2)
  tform <- c("loglin", "lin", "lin", "plin")
  modelform <- "M"
  a_n <- c(-0.75, 0.1, -0.05, -1.5)
  keep_constant <- c(0, 0, 0, 0)
  control <- list(
    "ncores" = 2, "lr" = 0.75, "maxiter" = 100, "halfmax" = 5, "epsilon" = 1e-3,
    "deriv_epsilon" = 1e-3, "abs_max" = 1.0, "change_all" = TRUE,
    "dose_abs_max" = 100.0, "verbose" = 0, "double_step" = 1
  )
  #

  e0 <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control)

  e <- RunPoissonEventAssignment_bound(df, pyr, event, e0, keep_constant, modelform, 4, 2, control)

  elow <- e$lower_limit$predict
  emid <- e$midpoint$predict
  eupp <- e$upper_limit$predict
  #
  expect_equal(sum(elow), 96.07807, tolerance = 1e-2)
  expect_equal(sum(emid), 123.6017, tolerance = 1e-2)
  expect_equal(sum(eupp), 151.1252, tolerance = 1e-2)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  #
  for (i in 2:4) {
    for (j in c(1, 2, 10)) {
      e <- RunPoissonEventAssignment_bound(df, pyr, event, e0, keep_constant, modelform, i, j, control)
      elow <- e$lower_limit$predict
      emid <- e$midpoint$predict
      eupp <- e$upper_limit$predict
      expect_equal(elow[, 1], emid[, 1], tolerance = 1e-2)
      expect_equal(elow[, 1], eupp[, 1], tolerance = 1e-2)
      expect_equal(eupp[, 1], emid[, 1], tolerance = 1e-2)
    }
  }
  #
})


test_that("Poisson Assigned Events, Epicure Validated", {
  df <- fread("Assigned_Check.csv")
  col_list <- c("a", "b")
  val <- factorize(df, col_list)
  df <- val$df
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  sum_values <- c(1.23715194448793, 0.060895797790186, 1.29804774227811, 1.67430308814, 0.0824135004185432, 1.75671658855854, 4.55122765982382, 0.224023120607428, 4.77525078043124, 12.3715194448793, 0.60895797790186, 12.9804774227811, 1.23715194448793, 0.478385445077363, 1.71553738956529, 1.67430308814, 0.478385445077363, 2.15268853321736, 4.55122765982382, 0.478385445077363, 5.02961310490118, 12.3715194448793, 0.478385445077363, 12.8499048899566, 1.23715194448793, 0.060895797790186, 1.29804774227811, 1.67430308814, 0.0824135004185432, 1.75671658855854, 4.55122765982382, 0.224023120607428, 4.77525078043124, 12.3715194448793, 0.60895797790186, 12.9804774227811, 1.23715194448793, -1.17625614669774, 0.060895797790186, 1.67430308814, -1.59188958772145, 0.0824135004185432, 4.55122765982382, -4.32720453921639, 0.224023120607428, 12.3715194448793, -11.7625614669774, 0.60895797790186)
  min_values <- c(0.1, 0.000230110913103737, 0.103747930428476, 0.135335283236613, 0.000311421256007305, 0.140407555497502, 0.367879441171442, 0.000846530741200557, 0.381667306687213, 1, 0.00230110913103743, 1.03747930428476, 0.1, 0.00170470357204963, 0.13669058432965, 0.135335283236613, 0.00170470357204963, 0.17281458752137, 0.367879441171442, 0.0017047035720496, 0.405358745456199, 1, 0.00170470357204966, 1.03747930428476, 0.1, 0.000230110913103737, 0.103747930428476, 0.135335283236613, 0.000311421256007305, 0.140407555497502, 0.367879441171442, 0.000846530741200557, 0.381667306687213, 1, 0.00230110913103743, 1.03747930428476, 0.1, -0.14873408300218, 0.000230110913103737, 0.135335283236613, -0.201289692500379, 0.000311421256007305, 0.367879441171442, -0.547162113379888, 0.000846530741200557, 1, -1.4873408300218, 0.00230110913103743)
  max_values <- c(0.164872127070013, 0.016138044067833, 0.181010171137846, 0.22313016014843, 0.0218404676480512, 0.244970627796481, 0.606530659712633, 0.0593685463327452, 0.665899206045379, 1.64872127070013, 0.16138044067833, 1.81010171137846, 0.164872127070013, 0.0978821851493432, 0.262754312219356, 0.22313016014843, 0.0978821851493432, 0.321012345297773, 0.606530659712633, 0.0978821851493432, 0.704412844861977, 1.64872127070013, 0.0978821851493432, 1.74660345584947, 0.164872127070013, 0.016138044067833, 0.181010171137846, 0.22313016014843, 0.0218404676480512, 0.244970627796481, 0.606530659712633, 0.0593685463327452, 0.665899206045379, 1.64872127070013, 0.16138044067833, 1.81010171137846, 0.164872127070013, -0.0960313868406136, 0.016138044067833, 0.22313016014843, -0.129964349376792, 0.0218404676480512, 0.606530659712633, -0.353279729258435, 0.0593685463327452, 1.64872127070013, -0.960313868406136, 0.16138044067833)
  event <- "e"
  t0 <- "t0"
  t1 <- "t1"
  pyr <- "pyr"
  names <- c("a", "b", "c", "CONST")
  j <- 1
  for (modelform in c("M", "A", "PAE", "PA")) {
    for (para in c(-2.3025850929940455, -2, -1, 0.0)) {
      names <- c("a", "b", "c", "CONST")
      Term_n <- rep(0, length(names))
      Term_n[3] <- 1
      tform <- rep("loglin", length(names))
      tform[3] <- "lin"
      a_n <- c(0.1, 0.15, 0.02, para)
      keep_constant <- rep(0, length(names))
      keep_constant[4] <- 1

      control <- list("Ncores" = 2, "maxiters" = c(1, 1), "verbose" = FALSE, "abs_max" = 0.1)
      model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "CR" = FALSE)
      e <- RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, control)
      e <- e$predict
      BK <- e[, 1]
      EX <- e[, 2]
      FV <- e[, 3]

      # BK values
      expect_equal(sum(BK), sum_values[j], tolerance = 1e-2)
      expect_equal(min(BK), min_values[j], tolerance = 1e-2)
      expect_equal(max(BK), max_values[j], tolerance = 1e-2)
      j <- j + 1

      # EX values
      expect_equal(sum(EX), sum_values[j], tolerance = 1e-2)
      expect_equal(min(EX), min_values[j], tolerance = 1e-2)
      expect_equal(max(EX), max_values[j], tolerance = 1e-2)
      j <- j + 1

      # FV values
      expect_equal(sum(FV), sum_values[j], tolerance = 1e-2)
      expect_equal(min(FV), min_values[j], tolerance = 1e-2)
      expect_equal(max(FV), max_values[j], tolerance = 1e-2)
      j <- j + 1
    }
  }
})
