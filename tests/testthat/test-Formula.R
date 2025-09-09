test_that("Basic factor application to formula", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 2, 1, 1, 2, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)

  model <- Cox(a, b, c) ~ loglinear(d, factor(e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(x = e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(e, levels = c(1, 2), labels = c("low", "high")))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_high"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, factor(e, levels = c(2, 1), labels = c("high", "low")))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_low"), e$names)
})

test_that("Basic ns and bs application to formula", {
  df <- data.table("a" = 1:100, "b" = 2:101, "c" = c(rep(0, 20), rep(1, 80)), "d" = c(rep(1, 20), rep(2, 50), rep(3, 30)), "e" = 0:99)
  model <- Cox(a, b, c) ~ loglinear(d, ns(e, df = 2))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_ns1", "e_ns2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, ns(e, intercept = TRUE))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_ns1", "e_ns2"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, bs(e))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_bs1", "e_bs2", "e_bs3"), e$names)
  model <- Cox(a, b, c) ~ loglinear(d, bs(e, Boundary.knots = c(0, 99)))
  e <- get_form(model, df)$model
  expect_equal(c("d", "e_bs1", "e_bs2", "e_bs3"), e$names)
})

test_that("Basic generic function application to formula", {
  df <- data.table("a" = 1:100, "b" = 2:101, "c" = c(rep(0, 20), rep(1, 80)), "d" = c(rep(1, 20), rep(2, 50), rep(3, 30)), "e" = 1:100)
  for (exp_string in c("log", "sqrt", "exp", "cos")) {
    model <- as.formula(paste("Cox(a, b, c) ~ loglinear(d, ", exp_string, "(e))", sep = ""))
    expect_no_error(e <- get_form(model, df))
  }
  for (exp_string in c("sqrt")) {
    model <- as.formula(paste("Cox(a, b, c) ~ loglinear(d, ", exp_string, "(e))", sep = ""))
    expect_no_error(e <- CoxRun(model, df, control = list("ncores" = 2)))
    expect_no_error(f <- RelativeRisk(e, df))
  }
})

test_that("Checking formula works with result modification", {
  fname <- "dose.csv"
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  #
  df$dose2 <- df$dose * df$dose
  df$a <- df$dose + 0.001
  df$b <- df$dose2 + 0.001
  #
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"

  model <- Cox(t0, t1, lung) ~ loglinear(ns(a), b)
  expect_no_error(e <- CoxRun(model, df, control = list("ncores" = 2)))
  f <- RelativeRisk(e, df)
  expect_equal(sum(f$Risk), 548.6874, tolerance = 1e-2)
  expect_no_error(e <- PoisRun(Pois(t1, lung) ~ loglinear(CONST, ns(b)), df, control = list("ncores" = 2)))
  f <- Residual(e, df)
  expect_equal(f$Residual_Sum, 0.5121325, tolerance = 1e-2)
})

test_that("Checking interaction works in formula and call results", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  #
  df$dose2 <- df$dose * df$dose
  df$a <- df$dose + 0.001
  df$b <- df$dose2 + 0.001
  df$rand0 <- floor(runif(nrow(df)) * 3)
  df$rand1 <- floor(runif(nrow(df)) * 2)
  #
  #
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"

  model <- Cox(t0, t1, lung) ~ loglinear(a, a * rand0, a * rand0 * rand1)
  e <- get_form(model, df)
  expect_equal(e$model$names, c("a", "rand0", "a:rand0", "rand1", "a:rand1", "rand0:rand1", "a:rand0:rand1"))
  expect_no_error(e <- CoxRun(model, df, control = list(ncores = 2)))
  f <- RelativeRisk(e, df)
  expect_equal(sum(f$Risk), 1607.914, tolerance = 1e-2)
})

test_that("Checking power works in formula and call results", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  #
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"

  model <- Cox(t0, t1, lung) ~ loglinear(dose, I(dose^2))
  e <- get_form(model, df)
  expect_equal(e$model$names, c("dose", "dose^2"))
  e <- CoxRun(model, df, control = list(ncores = 2))
  f <- RelativeRisk(e, df)
  expect_equal(sum(f$Risk), 862.3834, tolerance = 1e-2)
})
