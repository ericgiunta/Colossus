test_that("Basic formula and regression passes and fails", {
  df <- data.table(
    "a" = c(0, 0, 0, 1, 0, 1),
    "b" = c(1, 1, 1, 1, 1, 1),
    "d" = c(1, 2, 3, 4, 3, 2)
  )
  model <- logit(b, a) ~ loglinear(d)
  e <- get_form(model, df)$model
  expect_equal("b", e$trials)
  model <- logit(a) ~ loglinear(d)
  e <- get_form(model, df)$model
  expect_equal("CONST", e$trials)
  expect_no_error(LogisticRun(model, df, ncores = 1))
  expect_error(LogisticRun(model, df, control = "bad"))

  df <- data.table(
    "a" = c(0, 0, 0, 4, 0, 1),
    "b" = c(1, 1, 1, 1, 1, 1),
    "d" = c(1, 2, 3, 4, 3, 2)
  )
  model <- logit(b, a) ~ loglinear(d)
  expect_error(get_form(model, df))
  expect_error(LogisticRun(c(0), df, control = control, a_n = a_n))
})

test_that("basic regression with link non-fail", {
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


    control <- list(verbose = 0, step_max = 0.1, maxiter = 100, ncores = 1)
    #
    def_rate <- log(sum(df$status) / length(df$status))
    a_n <- c(0.001, -0.95, def_rate)
    model <- logit(status) ~ plinear(karno50) + loglinear(trt, CONST)
    e <- LogisticRun(model, df, control = control, a_n = a_n, verbose = 0)
    expect_equal(e$LogLik, -44.3606, tolerance = 1e-3)
    e <- LogisticRun(model, df, control = control, link = "odds", a_n = a_n)
    expect_equal(e$LogLik, -44.3606, tolerance = 1e-3)
    e <- LogisticRun(model, df, control = control, link = "ident", a_n = a_n)
    expect_equal(e$LogLik, -74.77318, tolerance = 1e-3)
    e <- LogisticRun(model, df, control = control, link = "loglink", a_n = a_n)
    expect_equal(e$LogLik, -91.35022, tolerance = 1e-3)
    #
    model <- logit(status) ~ plinear(karno50) + loglinear(trt, CONST)
    res <- get_form(model, df)
    expect_no_error(LogisticRun(res$model, df, control = control, a_n = a_n, norm = "max"))
    expect_no_error(LogisticRun(res$model, df, control = control, a_n = a_n, norm = "mean"))
    expect_error(LogisticRun(model, df, control = control, link = "bad", a_n = a_n))
    # observed_info
    expect_no_error(LogisticRun(res$model, df, control = control, a_n = a_n, observed_info = TRUE))
  }
})

test_that("epicure check", {
  df <- fread("sholom.csv", nThread = min(c(detectCores(), 2)), data.table = TRUE)

  control <- list(verbose = 0, ncores = 1)
  a_n <- c(0.4)
  model <- logit(n, x) ~ linear(CONST)
  e <- LogisticRun(model, df, control = control, a_n = a_n)
  expect_equal(e$beta_0, c(0.1222156), tolerance = 1e-4)

  a_n <- c(0.12, 0.1, 0.1)
  model <- logit(n, x) ~ linear(CONST, factor(alcohol))
  e <- LogisticRun(model, df, control = control, a_n = a_n)
  expect_equal(e$beta_0, c(0.09883952, 0.01229789, 0.11170645), tolerance = 1e-4)

  a_n <- c(0.12, 0.1, 0.1)
  model <- logit(n, x) ~ linear(CONST, 0) + linear(factor(alcohol), 1)
  e <- LogisticRun(model, df, control = control, a_n = a_n)
  expect_equal(e$beta_0, c(0.09881444, 0.12444199, 1.13051162), tolerance = 1e-4)

  a_n <- c(2)
  model <- logit(n, x) ~ linear(CONST)
  e <- LogisticRun(model, df, control = control, a_n = a_n, link = "loglink")
  expect_equal(e$beta_0, c(2.217405), tolerance = 1e-4)

  a_n <- c(2, 0.1, 0.7)
  model <- logit(n, x) ~ linear(CONST, factor(alcohol))
  e <- LogisticRun(model, df, control = control, a_n = a_n, link = "loglink")
  expect_equal(e$beta_0, c(2.4087301, -0.1061499, -0.6595241), tolerance = 1e-4)

  a_n <- c(0.01, 0.01, 0.01)
  model <- logit(n, x) ~ linear(CONST, factor(alcohol))
  e <- LogisticRun(model, df, control = control, a_n = a_n, link = "ident")
  expect_equal(e$beta_0, c(0.08989964, 0.01009402, 0.08397130), tolerance = 1e-4)
})
