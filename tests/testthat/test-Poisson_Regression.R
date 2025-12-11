test_that("Poisson time column missing", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 1, "lr" = 0.95, "maxiter" = -1, "halfmax" = 1, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 1.0, "verbose" = 0)
  expect_error(PoisRun(Poisson(a_bad, c) ~ loglinear(d, 0), df, control = control, a_n = a_n))
})
test_that("Poisson no events", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  fac <- c(0, 0, 0, 0, 1, 0, 0)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  keep_constant <- c(0)
  a_n <- c(-0.1)
  control <- list("ncores" = 1, "lr" = 0.95, "maxiter" = -1, "halfmax" = 1, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 1.0, "verbose" = 0)
  expect_error(PoisRun(Poisson(a, c) ~ loglinear(d, 0), df, control = control, a_n = a_n))
  #
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 1, 1, 0, 0, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  fac <- c(0, 0, 0, 0, 1, 0, 0)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "fac" = fac)
  expect_no_error(PoisRun(Poisson_Strata(a, c, fac) ~ loglinear(d, 0), df, control = control, a_n = a_n))
})

test_that("Pois loglin_M Single", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$pyr <- df$t1 - df$t0
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- PoisRun(Poisson(pyr, lung) ~ loglinear(dose, fac, 0), df, control = control, a_n = a_n, single = TRUE)
  expect_equal(e$AIC, 2354.055, tolerance = 1e-2)
})
test_that("Pois loglin_M Strata", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$pyr <- df$t1 - df$t0
  keep_constant <- c(0)
  a_n <- c(0.01)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- PoisRun(Poisson_Strata(pyr, lung, fac) ~ loglinear(dose, 0), df, control = control, a_n = a_n)
  expect_equal(e$beta_0, c(0.05476188), tolerance = 1e-1)
})

test_that("Checking pois strata default values", {
  if (system.file(package = "survival") != "") {
    data(cancer, package = "survival")
    veteran %>% setDT()
    df <- copy(veteran)
    # Make the same adjustments as Epicure example 6.5
    karno <- df$karno
    karno[93] <- 20
    df$karno <- karno / 100
    df$trt <- df$trt - 1
    df$trt <- as.integer(df$trt == 0)
    cell_lvl <- c("large", "squamous", "smallcell", "adeno")
    df$cell <- as.integer(factor(df$celltype, level = cell_lvl)) - 1
    df$karno50 <- df$karno - 50
    control <- list(ncores = 1, maxiter = 20, halfmax = 1)
    #
    a_n <- c(-1, 0.1, 0.1)
    model <- Pois_Strata(time, status, cell) ~ loglinear(CONST, trt, 0) + linear(karno, 1) + A()
    poisres <- PoisRun(model, df, a_n = a_n, control = control)
    expect_equal(poisres$beta_0, c(-0.4828601, -1.4264494, 1.1726181), tolerance = 1e-3)
    #
    a_n <- c(-0.4, -1, 1.17, -0.01)
    model <- Pois_Strata(time, status, cell) ~ loglinear(CONST, trt, 0) + loglin - dose(karno, 1) + PA()
    poisres <- PoisRun(model, df, a_n = a_n, control = control)
    expect_equal(poisres$beta_0, c(1.6661200, -0.2537878, 1.1700000, -2.7180987), tolerance = 1e-3)
    #
    a_n <- c(-1, 0.1, 0.1, 0.5)
    model <- Pois_Strata(time, status, cell) ~ loglinear(CONST, trt, 0) + linear - dose(karno, 1) + PAE()
    poisres <- PoisRun(model, df, a_n = a_n, control = control)
    expect_equal(poisres$beta_0, c(0.1685862, -0.1423654, -0.1216567, 0.1000000), tolerance = 1e-3)
    #
    a_n <- c(-1, 0.1, 0.1)
    model <- Pois_Strata(time, status, cell) ~ loglinear(CONST, trt, 0) + plinear(karno, 1) + GMIX()
    poisres <- PoisRun(model, df, a_n = a_n, control = control)
    expect_equal(poisres$beta_0, c(1.1841113, -0.2996062, -1.9301192), tolerance = 1e-3)
  }
})
