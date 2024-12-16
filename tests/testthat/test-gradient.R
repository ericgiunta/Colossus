test_that("Coxph strata_gradient_CR", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  df$censor <- (df$lung == 0)
  event <- "censor"
  names <- c("dose", "fac")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  plot_options <- list("name" = paste(tempfile(), "run", sep = ""), "verbose" = FALSE, "studyid" = "studyid", "age_unit" = "years")
  dft <- GetCensWeight(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)
  #
  #
  t_ref <- dft$t
  surv_ref <- dft$surv
  t_c <- df$t1
  cens_weight <- approx(t_ref, surv_ref, t_c, rule = 2)$y
  df$weighting <- cens_weight
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.001, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  verbose <- FALSE
  j_iterate <- 1
  for (i in c(TRUE, FALSE)) {
    for (k in c(TRUE)) {
      for (l in c(TRUE, FALSE)) {
        model_control <- list("strata" = i, "gradient" = k, "cr" = l)
        a_n <- c(-0.1, -0.1)
        control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
        modelform <- "M"
        e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
        expect_equal(e$Status, "PASSED")
        modelform <- "A"
        e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
        expect_equal(e$Status, "PASSED")
      }
    }
  }
})

test_that("Coxph gradient methods", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  df$censor <- (df$lung == 0)
  event <- "censor"
  names <- c("dose", "fac")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  modelform <- "M"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  #
  event <- "lung"
  keep_constant <- c(0, 0)
  model_control <- list("gradient" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  modelform <- "M"
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control))
  for (method in c("momentum", "adadelta", "adam")) {
    model_control <- list("gradient" = TRUE)
    model_control[[method]] <- TRUE
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    modelform <- "M"
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, fir = fir, der_iden = der_iden, control = control, strat_col = "rand", model_control = model_control)
    expect_equal(e$Status, "PASSED")
  }
})

test_that("Pois strata_gradient", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  event <- "lung"
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  names <- c("dose", "rand", "rand")
  term_n <- c(2, 1, 0)
  tform <- c("loglin", "lin", "plin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(0.01, 0.1, 0.1)
  modelform <- "PAE"
  fir <- 0
  der_iden <- 0
  control <- list("ncores" = 2, "lr" = 0.001, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 0)
  strat_col <- "fac"
  verbose <- FALSE
  j_iterate <- 1
  for (i in c(TRUE, FALSE)) {
    for (j in c(TRUE)) {
      model_control <- list("strata" = i, "gradient" = j)
      if (verbose) {
        print(model_control)
      }
      a_n <- c(0.01, 0.1, 0.1)
      modelform <- "PAE"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, strat_col, model_control)
      expect_equal(e$Status, "PASSED")
    }
  }
})
