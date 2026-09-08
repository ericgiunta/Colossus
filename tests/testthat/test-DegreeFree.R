test_that("Case-Control", {
  if (system.file(package = "survival") != "") {
    data(cancer, package = "survival")
    veteran |> setDT()
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


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(847.13843, 619.59108, 619.59108, 1125.09712, 918.33064, 918.33064, 61.45791, 52.75154, 49.85893, 63.53347, 63.53347, 63.53347)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)
    extra_bool <- "pass"
    for (time_bool in c(TRUE, FALSE)) {
      for (strat_bool in c(TRUE, FALSE)) {
        if (time_bool) {
          if (strat_bool) {
            model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon_Time(time, status) ~ loglinear(karno50, trt)
          }
        } else {
          if (strat_bool) {
            model <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon(status) ~ loglinear(karno50, trt)
          }
        }
        for (thres in c(0, 40, 100)) {
          e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, a_n = a_n)
          expect_equal(e$AIC, 2 * (e$FreeParameters + e$FreeSets) - 2 * e$LogLik, tolerance = 1e-3)
        }
      }
    }
    #
    for (time_bool in c(TRUE, FALSE)) {
      for (strat_bool in c(TRUE, FALSE)) {
        if (time_bool) {
          if (strat_bool) {
            model <- CaseCon_Strata_Time(time, status, cell) ~ null()
          } else {
            model <- CaseCon_Time(time, status) ~ null()
          }
        } else {
          if (strat_bool) {
            model <- CaseCon_Strata(status, cell) ~ null()
          } else {
            model <- CaseCon(status) ~ null()
          }
        }
        for (thres in c(0, 40, 100)) {
          e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, a_n = a_n)
          expect_equal(e$AIC, 2 * (e$FreeParameters + e$FreeSets) - 2 * e$LogLik, tolerance = 1e-3)
        }
      }
    }
  }
})

test_that("Cox Proportional Hazards", {
  fname <- "ll_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, verbose = FALSE, fill = TRUE)
  keep_constant <- 0
  a_n <- 0.01
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox_Strata(t0, t1, lung, fac) ~ loglinear(dose, 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
  e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, factor(fac), 0) + multiplicative(), df, a_n = a_n, control = control)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
  #
  e <- CoxRun(Cox_Strata(t0, t1, lung, fac) ~ null(), df, a_n = a_n, control = control)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
  e <- CoxRun(Cox(t0, t1, lung) ~ null(), df, a_n = a_n, control = control)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
})

test_that("Logistic", {
  df <- fread("sholom.csv", nThread = min(c(detectCores(), 2)), data.table = TRUE)
  control <- list(verbose = 0, ncores = 1)
  a_n <- c(0.12, 0.1, 0.1)
  model <- logit(n, x) ~ linear(CONST, factor(alcohol))
  e <- LogisticRun(model, df, control = control, a_n = a_n)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
  #
  model <- logit(n, x) ~ null()
  e <- LogisticRun(model, df, control = control)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) - 2 * e$LogLik, tolerance = 1e-3)
})

test_that("Poisson", {
  fname <- "ll_0.csv"
  col_types <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = col_types, verbose = FALSE, fill = TRUE)
  df$pyr <- df$t1 - df$t0
  keep_constant <- 0
  a_n <- 0.01
  control <- list("ncores" = 1, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- PoisRun(Poisson_Strata(pyr, lung, fac) ~ loglinear(dose, 0), df, control = control, a_n = a_n)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant) + e$strata_levels) + e$Deviance, tolerance = 1e-3)
  e <- PoisRun(Poisson(pyr, lung) ~ loglinear(dose, factor(fac), 0), df, control = control, a_n = a_n)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) + e$Deviance, tolerance = 1e-3)
  e <- PoisRun(Poisson_Strata(pyr, lung, fac) ~ null(), df, control = control, a_n = a_n)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant) + e$strata_levels - 1) + e$Deviance, tolerance = 1e-3)
  e <- PoisRun(Poisson(pyr, lung) ~ null(), df, control = control, a_n = a_n)
  expect_equal(e$AIC, 2 * (length(e$beta_0) - sum(e$model$keep_constant)) + e$Deviance, tolerance = 1e-3)
})
