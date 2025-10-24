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
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  # plot_options <- list("name" = paste(tempfile(), "run", sep = ""), "verbose" = FALSE, "studyid" = "studyid", "age_unit" = "years")
  # dft <- GetCensWeight(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, plot_options)
  # #
  # #
  # t_ref <- dft$t
  # surv_ref <- dft$surv
  # t_c <- df$t1
  # cens_weight <- approx(t_ref, surv_ref, t_c, rule = 2)$y
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
  # for (i in c(TRUE, FALSE)) {
  #   for (k in c(TRUE)) {
  #     for (l in c(TRUE, FALSE)) {
  #       model_control <- list("strata" = i, "gradient" = k, "cr" = l)
  #       a_n <- c(-0.1, -0.1)
  #       control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #       modelform <- "M"
  #       e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
  #       expect_equal(e$Status, "PASSED")
  #       modelform <- "A"
  #       e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
  #       expect_equal(e$Status, "PASSED")
  #     }
  #   }
  # }
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
  control <- list("ncores" = 2, "lr" = 0.001, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  event <- "lung"
  keep_constant <- c(0, 0)
  model_control <- list("gradient" = TRUE)
  a_n <- c(-0.1, -0.1)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  modelform <- "M"
  # expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control))
  for (method in c("momentum", "adadelta", "adam")) {
    model_control <- list("gradient" = TRUE)
    model_control[[method]] <- TRUE
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    modelform <- "M"
    gradient_control <- list()
    gradient_control[[method]] <- TRUE
    e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = gradient_control)
    # e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control)
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
  # names <- c("dose", "rand", "rand")
  # term_n <- c(2, 1, 0)
  # tform <- c("loglin", "lin", "plin")
  keep_constant <- c(0, 0)
  a_n <- c(0.01, 0.1)
  modelform <- "PAE"
  control <- list("ncores" = 2, "lr" = 0.001, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"
  verbose <- FALSE
  j_iterate <- 1
  e <- PoisRun(Poisson(pyr, lung) ~ loglinear(dose, rand, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  e <- PoisRun(Poisson_Strata(pyr, lung, fac) ~ loglinear(dose, rand, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  expect_equal(e$Status, "PASSED")
  # for (i in c(TRUE, FALSE)) {
  #   for (j in c(TRUE)) {
  #     model_control <- list("strata" = i, "gradient" = j)
  #     if (verbose) {
  #       print(model_control)
  #     }
  #     a_n <- c(0.01, 0.1, 0.1)
  #     modelform <- "PAE"
  #     e <- PoisRun(Pois(pyr, lung) ~ loglinear(dose, fac, 0) + m(), df, a_n = a_n, keep_constant = keep_constant, control = control, gradient_control = list())
  #     expect_equal(e$Status, "PASSED")
  #   }
  # }
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
    cell_string <- df$celltype
    cell <- case_when(
      cell_string == "squamous" ~ 1,
      cell_string == "smallcell" ~ 2,
      cell_string == "adeno" ~ 3,
      cell_string == "large" ~ 0
    )
    df$cell <- cell
    df$karno50 <- df$karno - 50
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2)
    #
    e2 <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control, cons_mat = cons_mat0, cons_vec = c(0.0), gradient_control = list("adam" = TRUE))
    expect_equal(e2$beta_0, c(0.158278, -0.158278), tolerance = 1e-3)
  }
})
