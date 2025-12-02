test_that("Coxph strata_gradient_CR", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  df$weighting <- df$t1 / 20
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.001, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  verbose <- FALSE
  j_iterate <- 1
  a_n <- c(-0.1, -0.1)
  e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  e <- CoxRun(Cox_Strata(t0, t1, lung, rand) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  e <- CoxRun(FineGray(t0, t1, lung, weighting) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  e <- CoxRun(FineGray_Strata(t0, t1, lung, rand, weighting) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
})

test_that("Coxph gradient methods", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  event <- "lung"
  keep_constant <- c(0, 0)
  model_control <- list("gradient" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  for (method in c("momentum", "adadelta", "adam")) {
    model_control <- list("gradient" = TRUE)
    model_control[[method]] <- TRUE
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    gradient_control <- list()
    gradient_control[[method]] <- TRUE
    e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = gradient_control)
    expect_equal(e$Status, "PASSED")
  }
  #
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  a_n <- c(1.0, 1.0)
  gradient_control <- list("adadelta" = TRUE)
  e1 <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, 0) + loglinear(fac, 1) + A(), df, a_n = a_n, control = control, gradient_control = gradient_control)
  e2 <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, 0) + loglinear(fac, 1) + ME(), df, a_n = a_n, control = control, gradient_control = gradient_control)
  e3 <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, 0) + loglinear(fac, 1) + GMIX(), df, a_n = a_n, control = control, gradient_control = gradient_control)
  expect_equal(e1$beta_0, c(0.5944258, 0.5397812), tolerance = 1e-3)
  expect_equal(e2$beta_0, c(0.6594536, 0.7276018), tolerance = 1e-3)
  expect_equal(e3$beta_0, c(0.6594536, 0.7276018), tolerance = 1e-3)
})

test_that("Pois strata_gradient", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$pyr <- df$t1 - df$t0
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  control <- list("ncores" = 2, "lr" = 0.001, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"
  verbose <- FALSE
  j_iterate <- 1
  e <- PoisRun(Poisson(pyr, lung) ~ loglinear(dose, rand, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  e <- PoisRun(Poisson_Strata(pyr, lung, fac) ~ loglinear(dose, rand, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
})

test_that("Logit Gradient Test", {
  if (system.file(package = "survival") != "") {
    data(cancer, package = "survival")
    veteran %>% setDT()
    df <- copy(veteran)

    # Make the same adjustments as Epicure example 6.5
    karno <- df$karno
    karno[93] <- 20
    df$karno <- karno
    df$trt <- df$trt - 1
    df$trt <- as.integer(df$trt == 0)
    cell_lvl <- c("large", "squamous", "smallcell", "adeno")
    df$cell <- as.integer(factor(df$celltype, level = cell_lvl)) - 1
    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)
    control <- list(verbose = 0, step_max = 0.1, maxiter = 100, ncores = 2)
    #
    def_rate <- log(sum(df$status) / length(df$status))
    a_n <- c(0.001, -0.95, def_rate)
    model <- logit(status) ~ plinear(karno50) + loglinear(trt, CONST)
    e <- LogisticRun(model, df, control = control, a_n = a_n, verbose = 0, gradient_control = list("adadelta" = TRUE))
    expect_equal(e$LogLik, -34.7568, tolerance = 1e-3)
  }
})

test_that("Constraint Check", {
  if (system.file(package = "survival") != "") {
    data(cancer, package = "survival")
    veteran %>% setDT()
    df <- copy(veteran)
    # Make the same adjustments as Epicure example 6.5
    karno <- df$karno
    karno[93] <- 20
    df$karno <- karno
    df$trt <- df$trt - 1
    df$trt <- as.integer(df$trt == 0)
    cell_lvl <- c("large", "squamous", "smallcell", "adeno")
    df$cell <- as.integer(factor(df$celltype, level = cell_lvl)) - 1
    df$karno50 <- df$karno - 50
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2)
    #
    e2 <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control, cons_mat = cons_mat0, cons_vec = c(0.0), gradient_control = list("adam" = TRUE))
    expect_equal(e2$beta_0, c(0.158278, -0.158278), tolerance = 1e-3)
    #
    e3 <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control, cons_mat = cons_mat0, cons_vec = c(0.0), gradient_control = list("momentum" = TRUE))
    e4 <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control, cons_mat = cons_mat0, cons_vec = c(0.0), gradient_control = list("adadelta" = TRUE))
    expect_equal(e3$beta_0, c(0.9, -0.9), tolerance = 1e-3)
    expect_equal(e4$beta_0, c(0.04037, -0.04037), tolerance = 1e-3)
    #
    cons_mat0 <- matrix(c(0, 0, 0, 0, 1), nrow = 1)
    e5 <- PoisRun(Pois(time, status) ~ loglinear(factor(cell, level = c(-1, 0, 1, 2, 3)), trt), df, control = control, cons_mat = cons_mat0, cons_vec = c(0.0), gradient_control = list("adam" = TRUE))
    expect_equal(e5$LogLik, -365.5018, tolerance = 1e-1)
  }
})
