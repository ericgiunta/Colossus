test_that("Basic factor application to formula", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 2, 1, 1, 2, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)

  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + PA()
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

test_that("Basic gmix application to formula", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 2, 1, 1, 2, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)

  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix(0.5, "e")
  expect_no_error(get_form(model, df))
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix - r(0.5)
  expect_no_error(get_form(model, df))
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix - e(0.5)
  expect_no_error(get_form(model, df))
  #
  model <- Cox(a, b, c) ~ loglinear(d, 0) + loglinear(factor(e), 1) + gmix()
  expect_no_error(get_form(model, df))
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix(0.5)
  expect_no_error(get_form(model, df))
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix("e")
  expect_no_error(get_form(model, df))
  #
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix - r()
  expect_no_error(get_form(model, df))
  model <- Cox(a, b, c) ~ loglinear(d, factor(e)) + gmix - e()
  expect_no_error(get_form(model, df))
})

test_that("Run basic errors and checks", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 2, 1, 1, 2, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)

  model <- Cox(a, b, c) ~ loglinear(d) ~ M()
  expect_error(CoxRun(model, df, ncores = 2))
  model <- Cox(a, b, c) ~ loglinear(d)
  expect_no_error(CoxRun(model, df, ncores = 2))
  expect_error(CoxRun(model, df, ncores = 2, bad = "wrong")) # arguement not in control list
  expect_error(CoxRun(model, df, control = c(2))) # control wasn't a list
  e <- get_form(Cox(tstart = a, event = c) ~ loglinear(d), df)
  model <- e$model
  expect_no_error(CoxRun(model, df, ncores = 2))
  #
  model <- Pois(b, c) ~ loglinear(d)
  expect_no_error(PoisRun(model, df, ncores = 2))
  expect_error(PoisRun(model, df, control = c(2))) # control wasn't a list
  expect_error(PoisRun(model, df, ncores = 2, bad = "wrong")) # arguement not in control list
  #
  model <- CaseCon(c) ~ loglinear(d)
  e <- get_form(model, df)
  model <- e$model
  expect_no_error(CaseControlRun(model, df, ncores = 2, keep_constant = c(0)))
  expect_error(CaseControlRun("bad", df, ncores = 2)) # wasn't a formula or model object
  expect_error(CaseControlRun(model, df, control = 2)) # control wasn't a list
  #
  model <- Pois(b, c) ~ loglinear(d)
  e <- get_form(model, df)
  model <- e$model
  expect_no_error(PoisRunJoint(model, df, ncores = 2))
  expect_error(PoisRunJoint("bad", df, ncores = 2)) # wasn't a formula or model object
  expect_error(PoisRunJoint(model, df, control = 2)) # control wasn't a list
  expect_no_error(PoisRunJoint(model, df, ncores = 2, a_n = c(0.1), keep_constant = c(0)))
  #
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
    model <- as.formula(paste("Cox(tend = b, event = c) ~ loglinear(d, ", exp_string, "(e))", sep = ""))
    expect_no_error(e <- CoxRun(model, df, control = list("ncores" = 2)))
    expect_no_error(f <- RelativeRisk(e, df))
  }
  model <- Cox(a, b, c) ~ loglinear(factor(d))
  expect_no_error(e <- CoxRun(model, df, ncores = 2))
  expect_no_error(f <- RelativeRisk(e, df))
})

test_that("Basic factor application to formula with formula column", {
  df <- data.table("a" = 1:100, "b" = 2:101, "c" = c(rep(0, 20), rep(1, 80)), "d" = c(rep(1, 20), rep(2, 50), rep(3, 30)), "e" = 1:100)
  model <- Cox(a, b, c) ~ loglinear(factor(d))
  expect_no_error(e <- CoxRun(model, df, ncores = 2))
  expect_no_error(f <- RelativeRisk(e, df))
  df$d <- factor(df$d)
  model <- Cox(a, b, c) ~ loglinear(d)
  expect_no_error(e <- CoxRun(model, df, ncores = 2))
  expect_no_error(f <- RelativeRisk(e, df))
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
  model <- Cox(t0, t1, lung) ~ loglinear(ns(x = a), b)
  expect_no_error(e <- CoxRun(model, df, control = list("ncores" = 2)))
  f <- RelativeRisk(e, df)
  expect_equal(sum(f$Risk), 548.6874, tolerance = 1e-2)
  expect_no_error(e <- PoisRun(Pois(t1, lung) ~ loglinear(CONST, bs(x = b)), df, control = list("ncores" = 2)))
  f <- Residual(e, df)
  expect_equal(f$Residual_Sum, 0.497, tolerance = 1e-2)
  #
  model <- Cox(t0, t1, lung) ~ loglinear(bs(a), b)
  expect_no_error(e <- CoxRun(model, df, control = list("ncores" = 2)))
  expect_no_error(f <- RelativeRisk(e, df))
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
  df$rand1 <- factor(floor(runif(nrow(df)) * 2))
  #
  #
  model <- Cox(t0, t1, lung) ~ loglinear(a, a * rand0, a * rand0 * rand1)
  e <- get_form(model, df)
  expect_equal(e$model$names, c("a", "rand0", "a:rand0", "rand1_1", "a:rand1_1", "rand0:rand1_1", "a:rand0:rand1_1"))
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
  model <- Cox(t0, t1, lung) ~ loglinear(dose, I(dose^2))
  e <- get_form(model, df)
  expect_equal(e$model$names, c("dose", "dose^2"))
  e <- CoxRun(model, df, control = list(ncores = 2))
  f <- RelativeRisk(e, df)
  expect_equal(sum(f$Risk), 862.3834, tolerance = 1e-2)
})

test_that("Checking model and formula can be input", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  #
  model <- Cox(t0, t1, lung) ~ loglinear(dose, I(dose^2))
  e <- get_form(model, df)
  expect_no_error(e0 <- CoxRun(model, df, control = list(ncores = 2)))
  expect_no_error(e1 <- CoxRun(e$model, df, control = list(ncores = 2)))
  expect_equal(e0$beta_0, e1$beta_0)

  model <- Pois(t1, lung) ~ loglinear(CONST, dose, I(dose^2))
  e <- get_form(model, df)
  expect_no_error(e0 <- PoisRun(model, df, control = list(ncores = 2)))
  expect_no_error(e1 <- PoisRun(e$model, df, control = list(ncores = 2)))
  expect_equal(e0$beta_0, e1$beta_0)
})

test_that("Joint Form Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$rand1 <- floor(runif(nrow(df)) * 5)
  #
  model0 <- Pois(t1, lung) ~ loglinear(dose)
  model1 <- Pois(t1, lung) ~ loglinear(I(dose^2))
  models <- Pois(t1, lung) ~ loglinear(t0)
  expect_no_error(get_form_joint(list(model0, model1, "shared" = models), df))
  expect_error(get_form_joint("bad", df)) # didn't pass a list or formula
  model0 <- Cox(t1, lung) ~ loglinear(dose)
  model1 <- Pois(t1, lung) ~ loglinear(I(dose^2))
  models <- Pois(t1, lung) ~ loglinear(t0)
  expect_error(get_form_joint(list(model0, model1, "shared" = models), df)) # didn't use all pois models
  model0 <- Pois(t1, lung) ~ loglinear(dose)
  model1 <- Pois(t0, lung) ~ loglinear(I(dose^2))
  models <- Pois(t1, lung) ~ loglinear(t0)
  expect_error(get_form_joint(list(model0, model1, "shared" = models), df)) # didn't use same person-year time
  model0 <- Pois_Strata(t1, lung, rand0) ~ loglinear(dose)
  model1 <- Pois_Strata(t1, lung, rand0) ~ loglinear(I(dose^2))
  models <- Pois_Strata(t1, lung, rand0) ~ loglinear(t0)
  expect_no_error(get_form_joint(list(model0, model1, "shared" = models), df)) # same strata
  model0 <- Pois_Strata(t1, lung, rand0) ~ loglinear(dose)
  model1 <- Pois_Strata(t0, lung, rand1) ~ loglinear(I(dose^2))
  models <- Pois_Strata(t1, lung, rand0) ~ loglinear(t0)
  expect_error(get_form_joint(list(model0, model1, "shared" = models), df)) # different strata
  model0 <- Pois(t1, lung) ~ loglinear(dose) + M()
  model1 <- Pois(t1, lung) ~ loglinear(I(dose^2)) + A()
  models <- Pois(t1, lung) ~ loglinear(t0) + M()
  expect_error(get_form_joint(list(model0, model1, "shared" = models), df)) # different modelform
  model0 <- Pois(t1, lung) ~ loglinear(dose) + M()
  model1 <- Pois(t1, lung) ~ loglinear(I(dose^2)) + M()
  models <- Pois(t1, lung) ~ loglinear(t0) + A()
  expect_error(get_form_joint(list(model0, model1, "shared" = models), df)) # different modelform
})

test_that("General Form Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$rand1 <- floor(runif(nrow(df)) * 5)
  #
  expect_no_error(get_form(Cox(t0, t1, lung) ~ loglinear(dose), df))
  expect_error(get_form(bad ~ loglinear(dose), df)) # no ( on left side
  expect_error(get_form(also_bad(t0) ~ loglinear(dose), df)) # not defined left side
  #
  expect_error(get_form(Cox(t0, t1, lung) ~ dose, df)) # not defined right side term
  expect_error(get_form(Cox(t0, t1, lung) ~ loglinear(I(dose^a)), df)) # not a numeric power
  expect_error(get_form(Cox(t0, t1, lung) ~ loglinear(dose * a), df)) # missing interaction column
  #
  expect_error(get_form(Cox(t0, t1, lung) ~ loglinear(dose) + M() + A(), df)) # modelform defined twice
  expect_error(get_form(Cox(t0, t1, lung) ~ loglinear(dose) + weird(), df)) # unkown term with ()
})

test_that("Colossus Surv Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$rand1 <- floor(runif(nrow(df)) * 5)
  df$weight <- df$t1 / 100
  #
  expect_no_error(get_form(Cox(t0, t1, lung) ~ loglinear(dose), df))
  expect_no_error(get_form(Cox(t0, t1, event = lung) ~ loglinear(dose), df))
  expect_no_error(get_form(Cox(tstart = t0, tend = t1, event = lung) ~ loglinear(dose), df))
  expect_error(get_form(Cox(lung) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Cox(tstart = t0, tend = t1) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Cox(t0, t1, lung, lung) ~ loglinear(dose), df)) # too many
  expect_error(get_form(Cox(alpha = t0, t1, lung) ~ loglinear(dose), df)) # wrong named
  #
  expect_no_error(get_form(Cox_Strata(t0, t1, lung, rand0) ~ loglinear(dose), df))
  expect_no_error(get_form(Cox_Strata(t0, t1, lung, strata = rand0) ~ loglinear(dose), df))
  expect_error(get_form(Cox_Strata(lung, strata = rand0) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Cox_Strata(t0, t1, lung, lung, strata = rand0) ~ loglinear(dose), df)) # too many
  expect_error(get_form(Cox_Strata(t0, t1, lung, alpha = rand0) ~ loglinear(dose), df)) # wrong named
  #
  expect_no_error(get_form(FineGray(t0, t1, lung, weight) ~ loglinear(dose), df))
  expect_no_error(get_form(FineGray(t0, t1, lung, weight = weight) ~ loglinear(dose), df))
  expect_error(get_form(FineGray(lung, weight = weight) ~ loglinear(dose), df)) # too few
  expect_error(get_form(FineGray(t0, t1, lung, lung, weight = weight) ~ loglinear(dose), df)) # too many
  expect_error(get_form(FineGray(t0, t1, lung, alpha = weight) ~ loglinear(dose), df)) # wrong named
  #
  expect_no_error(get_form(FineGray_Strata(t0, t1, lung, rand0, weight) ~ loglinear(dose), df))
  expect_no_error(get_form(FineGray_Strata(t0, t1, lung, rand0, weight = weight) ~ loglinear(dose), df))
  expect_error(get_form(FineGray_Strata(lung, rand0, weight = weight) ~ loglinear(dose), df)) # too few
  expect_error(get_form(FineGray_Strata(t0, t1, lung, lung, rand0, weight = weight) ~ loglinear(dose), df)) # too many
  expect_error(get_form(FineGray_Strata(t0, t1, lung, rand0, alpha = weight) ~ loglinear(dose), df)) # wrong named
  #
  #
  expect_no_error(get_form(Pois(t1, lung) ~ loglinear(dose), df))
  expect_error(get_form(Pois(lung) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Pois(t1, lung, rand0) ~ loglinear(dose), df)) # too many
  expect_error(get_form(Pois(alpha = t1, lung) ~ loglinear(dose), df)) # wrong named
  #
  expect_no_error(get_form(Pois_Strata(pyr = t1, event = lung, rand0, rand1) ~ loglinear(dose), df))
  expect_error(get_form(Pois_Strata(pyr = t1, event = lung) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Pois_Strata(t1, lung) ~ loglinear(dose), df)) # too few
  expect_error(get_form(Pois_Strata(alpha = t1, lung, rand0) ~ loglinear(dose), df)) # wrong named
  expect_error(get_form(Pois_Strata(pyr = t1, event = lung, rand0, rand1) ~ linear(dose), df)) # can't correct with a default value of 0
  #
})

test_that("Pois multi_surv nonerror", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$rand1 <- floor(runif(nrow(df)) * 5)
  df$weight <- df$t1 / 100
  #
  expect_no_error(get_form(Pois_Strata(pyr = t1, event = lung, rand0, rand1) ~ loglinear(dose), df))
  expect_no_error(get_form(Pois_Strata(t1, lung, rand0, rand1) ~ loglinear(dose), df))
  #
  df$rand0 <- factor(df$rand0)
  df$rand1 <- factor(df$rand1)
  expect_no_error(get_form(Pois_Strata(t1, lung, rand0) ~ loglinear(dose), df))
  expect_no_error(get_form(Pois_Strata(t1, lung, rand0, rand1) ~ loglinear(dose), df))
})

test_that("CaseCon Surv Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$rand1 <- floor(runif(nrow(df)) * 5)
  df$weight <- df$t1 / 100
  #
  expect_no_error(get_form(CaseCon(lung) ~ loglinear(dose), df))
  expect_error(get_form(CaseCon() ~ loglinear(dose), df)) # too few
  expect_error(get_form(CaseCon(t1, lung) ~ loglinear(dose), df)) # too many
  #
  expect_no_error(get_form(CaseCon_time(t0, t1, lung) ~ loglinear(dose), df))
  expect_error(get_form(CaseCon_time(t1) ~ loglinear(dose), df)) # too few
  expect_error(get_form(CaseCon_time(t0, t1, t1, lung) ~ loglinear(dose), df)) # too many
  expect_error(get_form(CaseCon_time(alpha = t0, t1, lung) ~ loglinear(dose), df)) # wrong name
  # Getting the named entries, three case
  expect_no_error(get_form(CaseCon_time(tstart = t0, tend = t1, lung) ~ loglinear(dose), df))
  expect_no_error(get_form(CaseCon_time(tstart = t0, t1, event = lung) ~ loglinear(dose), df))
  expect_no_error(get_form(CaseCon_time(t0, tend = t1, event = lung) ~ loglinear(dose), df))
  # Getting the named entries, two case
  expect_no_error(get_form(CaseCon_time(tend = t1, lung) ~ loglinear(dose), df))
  expect_no_error(get_form(CaseCon_time(t1, event = lung) ~ loglinear(dose), df))
  #
  expect_no_error(get_form(CaseCon_strata(lung, rand0) ~ loglinear(dose), df))
  expect_error(get_form(CaseCon_strata(lung) ~ loglinear(dose), df)) # too few
  expect_error(get_form(CaseCon_strata(t1, lung, rand0) ~ loglinear(dose), df)) # too many
  expect_error(get_form(CaseCon_strata(alpha = lung, rand0) ~ loglinear(dose), df)) # wrong name
  #
  expect_no_error(get_form(CaseCon_strata_time(t0, t1, lung, rand0) ~ loglinear(dose), df))
  expect_error(get_form(CaseCon_strata_time(lung, rand0) ~ loglinear(dose), df)) # too few
  expect_error(get_form(CaseCon_strata_time(t0, t0, t1, lung, rand0) ~ loglinear(dose), df)) # too many
  expect_error(get_form(CaseCon_strata_time(t0, t1, lung, alpha = rand0) ~ loglinear(dose), df)) # wrong name
  #
})

test_that("Object Validation Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$col_bad <- "a"
  #
  control <- list(ncores = 2, maxiter = -1, maxiters = c(-1, -1))
  #
  true_cox <- get_form(Cox(t0, t1, lung) ~ loglinear(dose), df)
  cox_model <- copy(true_cox$model)
  expect_no_error(CoxRun(cox_model, df, control = control))
  true_pois <- get_form(Pois(t1, lung) ~ loglinear(dose), df)
  pois_model <- copy(true_pois$model)
  expect_no_error(CoxRun(cox_model, df, control = control))
  ##
  cox_model <- copy(true_cox$model)
  cox_model$names <- c("col_bad")
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  expect_error(CoxRun(cox_model, df, control = control, keep_constant = c(-1)))
  expect_error(CoxRun(cox_model, df, control = control, keep_constant = c(2)))
  expect_error(CoxRun(cox_model, df, control = control, keep_constant = c("a")))
  cox_model <- copy(true_cox$model)
  cox_model$tform <- c("bad_bad")
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  expect_error(CoxRun(cox_model, df, a_n = c(1, 1, 1, 1, 1), control = control))
  cox_model <- copy(true_cox$model)
  cox_model$names <- c("dose2dose")
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  cox_model$modelform <- "weird"
  expect_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_pois$model)
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  cox_model$event <- ""
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  cox_model$start_age <- "bad"
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  cox_model$end_age <- "bad"
  expect_error(CoxRun(cox_model, df, control = control))
  cox_model <- copy(true_cox$model)
  cox_model$event <- "bad"
  expect_error(CoxRun(cox_model, df, control = control))
  #
  pois_model <- copy(true_cox$model)
  expect_error(PoisRun(pois_model, df, control = control))
  pois_model <- copy(true_pois$model)
  pois_model$event <- ""
  expect_error(PoisRun(pois_model, df, control = control))
  pois_model <- copy(true_pois$model)
  pois_model$person_year <- "bad"
  expect_error(PoisRun(pois_model, df, control = control))
  pois_model <- copy(true_pois$model)
  pois_model$event <- "bad"
  expect_error(PoisRun(pois_model, df, control = control))
  #
})

test_that("Multiplicative model check", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$col_bad <- "a"
  control <- list(ncores = 2)
  #
  model <- Cox(t0, t1, lung) ~ loglinear(dose, rand0, 0) + M()
  res0 <- CoxRun(model, df, control = control)
  #
  model <- Cox(t0, t1, lung) ~ loglinear(dose, 0) + loglinear(rand0, 1) + M()
  res1 <- CoxRun(model, df, control = control)
  #
  expect_equal(res0$beta_0, res1$beta_0, tolerance = 1e-2)
})

test_that("Formula Validation Errors", {
  fname <- "dose.csv"
  set.seed(3742)
  colTypes <- c("double", "double", "double", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$rand0 <- floor(runif(nrow(df)) * 5)
  df$col_bad <- "a"
  #
  control <- list(ncores = 2, maxiter = -1, maxiters = c(-1, -1))
  #
  true_cox <- get_form(Cox(t0, t1, lung) ~ loglinear(dose, 0) + loglinear(rand0, 1), df)
  cox_model <- copy(true_cox$model)
  # Term_n errors
  cox_model <- copy(true_cox$model)
  cox_model$term_n <- c("a", "b")
  expect_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$term_n <- c(1, 3)
  expect_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$term_n <- c(0)
  expect_no_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$term_n <- c(0, 1, 2, 3)
  expect_no_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$keep_constant <- c(0.1, 0, 0, 0, 0, 0, 0)
  expect_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$tform <- c("loglin")
  expect_no_error(CoxRun(cox_model, df, control = control))
  #
  cox_model <- copy(true_cox$model)
  cox_model$tform <- c("loglin", "loglin", "loglin")
  expect_no_error(CoxRun(cox_model, df, control = control))
  #
  a_n <- list(c(0.1, 0.1))
  cox_model <- copy(true_cox$model)
  expect_no_error(CoxRun(cox_model, df, control = control, a_n = a_n))
  #
  a_n <- list(c(0.1, 0.1), c(0.1, 0.1, 0.2))
  expect_error(CoxRun(cox_model, df, control = control, a_n = a_n))
  a_n <- list(c(0.1), c(0.1))
  expect_no_error(CoxRun(cox_model, df, control = control, a_n = a_n))
  a_n <- list(c(0.1, 0.1, 0.2), c(0.1, 0.1, 0.2))
  expect_error(CoxRun(cox_model, df, control = control, a_n = a_n))
})

test_that("Basic formula passes and fails", {
  df <- data.table(
    "a" = c(0, 0, 0, 1, 0, 1),
    "b" = c(1, 1, 1, 1, 1, 1),
    "d" = c(1, 2, 3, 4, 3, 2)
  )
  expect_no_error(get_form(logit(b, a) ~ loglinear(d), df)) # all good
  expect_error(get_form(logit() ~ loglinear(d), df)) # no entries
  expect_error(get_form(logit(b, a, a) ~ loglinear(d), df)) # too many entries
  expect_error(get_form(logit(b, bad = a) ~ loglinear(d), df)) # bad named entry
  expect_error(get_form(logit(event = b, events = a) ~ loglinear(d), df)) # cannot do both
  # now matching correctly
  expect_no_error(get_form(logit(trials = b, a) ~ loglinear(d), df)) # all good
  expect_no_error(get_form(logit(b, event = a) ~ loglinear(d), df)) # all good
  expect_no_error(get_form(logit(b, events = a) ~ loglinear(d), df)) # all good
})
