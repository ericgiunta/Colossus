test_that("Coxph basic_single_null match", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("dose")
  term_n <- c(0)
  tform <- c("loglin")
  keep_constant <- c(0)
  a_n <- c(0.0)
  modelform <- "ME"

  verbose <- FALSE

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "null" = FALSE)
  e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "ME", control = control, strat_col = "fac", model_control = model_control)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (j in c(TRUE, FALSE)) {
    for (k in c(TRUE, FALSE)) {
      for (l in c(TRUE, FALSE)) {
        model_control <- list("strata" = FALSE, "basic" = j, "single" = k, "null" = l)
        a_n <- c(0.0)
        e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "ME", control = control, strat_col = "fac", model_control = model_control)
        expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-2)
      }
    }
  }
  model_control <- list("strata" = TRUE, "basic" = FALSE, "single" = FALSE, "null" = FALSE)
  e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "ME", control = control, strat_col = "fac", model_control = model_control)
  for (j in c(TRUE, FALSE)) {
    for (k in c(TRUE, FALSE)) {
      for (l in c(TRUE, FALSE)) {
        model_control <- list("strata" = TRUE, "basic" = j, "single" = k, "null" = l)
        a_n <- c(0.0)
        e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "ME", control = control, strat_col = "fac", model_control = model_control)
        expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-2)
      }
    }
  }
})
test_that("Pois strata_single", {
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  verbose <- FALSE
  j_iterate <- 1
  LL_comp <- c(-468.7465, -464.8984, -462.4579, -462.4461, -3033.332, -2734.64, -1104.25, -1368.039)
  for (i in c(TRUE, FALSE)) {
    for (j in c(TRUE, FALSE)) {
      model_control <- list("strata" = i, "single" = j)
      a_n <- c(0.01, 0.1, 0.1)
      modelform <- "PAE"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
      expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
      j_iterate <- j_iterate + 1
      a_n <- c(0.01, 0.1, 0.1)
      modelform <- "A"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
      expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
      j_iterate <- j_iterate + 1
    }
  }
})
#
test_that("Pois comb_forms", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  event <- "lung"
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 1, max = 5))
  names <- c("dose", "rand", "rand", "dose", "dose")
  term_n <- c(1, 0, 0, 0, 0)
  tform <- c("loglin", "lin", "plin", "loglin_slope", "loglin_top")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(0.01, 0.1, 0.1, 1.0, 0.1)
  modelform <- "PAE"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"

  verbose <- FALSE
  modelforms <- c("A", "PAE", "ME", "PA")
  j_iterate <- 1
  LL_comp <- c(-1644.494, -544.7434, -544.7434, -464.709, -1395.197, -1831.403, -1546.554, -464.709)
  for (modelform in modelforms) {
    model_control <- list("strata" = FALSE, "single" = FALSE)
    a_n <- c(0.01, 0.1, 0.1, 1.0, 0.1)
    e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
    expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
    j_iterate <- j_iterate + 1
  }
  term_n <- c(1, 1, 1, 0, 0)
  for (modelform in modelforms) {
    model_control <- list("strata" = FALSE, "single" = FALSE)
    a_n <- c(0.01, 0.1, 0.1, 1.0, 0.1)
    e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
    expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
    j_iterate <- j_iterate + 1
  }
})
#
test_that("Pois strata_single expanded", {
  fname <- "ll_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  event <- "lung"
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  names <- c("dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand")
  term_n <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  tform <- c("loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope")
  keep_constant <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  a_n <- c(-0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)

  modelform <- "PAE"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  verbose <- FALSE
  j_iterate <- 1
  LL_comp <- c(-496.7366, -475.4213, -496.7366, -475.4213, -4497.178, -3577.953, -4304.506, -2590.778)
  for (i in c(TRUE, FALSE)) {
    for (j in c(TRUE, FALSE)) {
      model_control <- list("strata" = i, "single" = j)
      a_n <- c(-0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)
      modelform <- "PAE"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
      expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
      j_iterate <- j_iterate + 1
      a_n <- c(-0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)
      modelform <- "A"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
      expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
      j_iterate <- j_iterate + 1
    }
  }
})

test_that("risk check omnibus plain", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  #
  event <- "lung"
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 1)
  tform <- c("loglin", "lin", "lin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  modelform <- "ME"


  cens_weight <- c(0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  coxres <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, 0) + plinear(fac, 0) + linear(dose, 1) + plinear(fac, 1) + loglinear(rand, 1), df, a_n = a_n, control = control)
  verbose <- FALSE

  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)

  model_list <- c("ME", "A", "PA", "PAE")
  means <- c(1.039780, 1.074573, 0.1434144, 1.039780)
  medians <- c(1.017646, 1.049693, 0.0904543, 1.017646)
  sums <- c(1039.779856, 1074.573350, 143.4143522, 1039.779856)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (model_i in 1:4) {
    modelform <- model_list[model_i]
    #
    coxres$model$modelform <- modelform
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RelativeRisk(coxres, df)$Risk

    expect_equal(mean(e), means[model_i], tolerance = 1e-2)
    expect_equal(median(e), medians[model_i], tolerance = 1e-2)
    expect_equal(sum(e), sums[model_i], tolerance = 1e-2)
  }
})

test_that("risk check omnibus gmix", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  #
  event <- "lung"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE

  model_list <- c("GMIX-R", "GMIX-E", "GMIX")
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 2)
  tform <- c("loglin", "loglin", "plin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)

  means <- c(1.068010, 4.181047, 1.068010, 2.112403, 2.113895, 4.181047)
  medians <- c(1.063724, 4.167848, 1.063724, 2.109973, 2.105042, 4.167848)
  sums <- c(1068.009705, 4181.046716, 1068.009705, 2112.403365, 2113.894860, 4181.046716)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  cens_weight <- c(0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  coxres <- CoxRun(Cox(t0, t1, lung) ~ loglinear(dose, fac, 0) + plinear(dose, fac, 1) + loglinear(rand, 2) + gmix(1.0, e, e, e), df, control = control)

  model_list <- c("GMIX-R", "GMIX-E", "GMIX")
  count <- 0
  for (model_i in 1:3) {
    modelform <- model_list[model_i]
    if (modelform == "GMIX") {
      for (term_i in 0:3) {
        gmix_term <- c(0, term_i %% 2, floor(term_i / 2))
        coxres$modelcontrol$gmix_term <- gmix_term
        #
        control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
        e <- e <- RelativeRisk(coxres, df)$Risk
        count <- count + 1
        expect_equal(mean(e), means[count], tolerance = 1e-2)
        expect_equal(median(e), medians[count], tolerance = 1e-2)
        expect_equal(sum(e), sums[count], tolerance = 1e-2)
      }
    } else if (modelform == "GMIX-R") {
      coxres$modelcontrol$gmix_term <- c(0, 0, 0)
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
      e <- RelativeRisk(coxres, df)$Risk
      count <- count + 1
      expect_equal(mean(e), means[count], tolerance = 1e-2)
      expect_equal(median(e), medians[count], tolerance = 1e-2)
      expect_equal(sum(e), sums[count], tolerance = 1e-2)
    } else {
      coxres$modelcontrol$gmix_term <- c(1, 1, 1)
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
      e <- RelativeRisk(coxres, df)$Risk
      count <- count + 1
      expect_equal(mean(e), means[count], tolerance = 1e-2)
      expect_equal(median(e), medians[count], tolerance = 1e-2)
      expect_equal(sum(e), sums[count], tolerance = 1e-2)
    }
  }
})
#
test_that("check deviation calc, expected cox", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("dose", "fac", "rand")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)
  modelform <- "ME"


  cens_weight <- c(0)

  verbose <- FALSE

  devs <- c()

  modelform <- "ME"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  for (i in 1:3) {
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
  devs <- c(devs, sum(e$Standard_Deviation))

  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("dose", "fac", "rand")
  term_n <- c(0, 0, 0)
  tform <- c("loglin", "loglin", "plin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)

  for (i in 1:3) {
    a_n <- c(0.6428582, 0.4240752, 0.1507817)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6428582, 0.4240752, 0.1507817)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
  devs <- c(devs, sum(e$Standard_Deviation))

  expect_equal(devs, c(0.61445, 0.54101, 0.73858, 0.95015, 0.63646, 0.56292, 0.73815, 0.97195), tolerance = 1e-4)
})
test_that("check deviation calc, Observed cox", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("dose", "fac", "rand")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)
  modelform <- "ME"


  cens_weight <- c(0)

  verbose <- FALSE

  devs <- c()

  modelform <- "ME"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE, "observed_info" = TRUE)
  for (i in 1:3) {
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
  devs <- c(devs, sum(e$Standard_Deviation))

  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("dose", "fac", "rand")
  term_n <- c(0, 0, 0)
  tform <- c("loglin", "loglin", "plin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)

  for (i in 1:3) {
    a_n <- c(0.6428582, 0.4240752, 0.1507817)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6428582, 0.4240752, 0.1507817)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
  devs <- c(devs, sum(e$Standard_Deviation))
  expect_equal(devs, c(0.6091269, 0.5356671, 0.7385757, 0.9448081, 0.7051473, 0.5838560, 0.7381538, 0.9897501), tolerance = 1e-4)
})

test_that("check Linear Constraints", {
  fname <- "l_pl_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  df$pyr <- df$t1 - df$t0
  time1 <- "t0"
  time2 <- "t1"
  pyr <- "pyr"
  event <- "lung"
  names <- c("dose", "fac")
  term_n <- c(0, 0)
  tform <- c("loglin", "plin")
  keep_constant <- c(0, 0)
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "null" = FALSE, "constraint" = TRUE)
  Constraint_Matrix <- matrix(c(1, -1), nrow = 1)
  Constraint_const <- c(0.0)
  set.seed(3742)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (i in 1:5) {
    a_n <- 2 * runif(2) - 1
    del <- abs(a_n[1] - a_n[2])
    a_n0 <- rep(sum(a_n) / 2, 2)
    a_n <- a_n0 - c(-del / 2, del / 2)
    modelform <- "ME"
    control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col = "fac", model_control = model_control, cons_mat = Constraint_Matrix, cons_vec = Constraint_const)
    expect_equal(e$beta_0, c(0.357333, 0.357333), tolerance = 1e-2)
  }
  for (i in 1:5) {
    a_n <- 2 * runif(2) - 1
    del <- abs(a_n[1] - a_n[2])
    a_n0 <- rep(sum(a_n) / 2, 2)
    a_n <- a_n0 + c(-del / 2, del / 2)
    modelform <- "ME"
    control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col = "fac", model_control = model_control, cons_mat = Constraint_Matrix, cons_vec = Constraint_const)
    expect_equal(e$beta_0, c(-0.472812, -0.472812), tolerance = 1e-2)
  }
})
test_that("check deviation calc, poisson", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  event <- "lung"
  names <- c("dose", "fac", "rand")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)
  modelform <- "ME"
  devs <- c()

  modelform <- "ME"
  for (inma_type in c(T, F)) {
    model_control <- list("observed_info" = inma_type)
    for (i in 1:3) {
      a_n <- c(0.6465390, 0.4260961, 0.1572781)
      keep_constant <- c(0, 0, 0)
      keep_constant[i] <- 1
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)
      devs <- c(devs, sum(e$Standard_Deviation))
    }
    event <- "lung"
    names <- c("dose", "fac", "rand")
    term_n <- c(0, 0, 0)
    tform <- c("loglin", "loglin", "plin")
    keep_constant <- c(0, 0, 0)
    a_n <- c(-0.1, 0.1, 0.2)

    for (i in 1:3) {
      a_n <- c(0.6428582, 0.4240752, 0.1507817)
      keep_constant <- c(0, 0, 0)
      keep_constant[i] <- 1
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)
      devs <- c(devs, sum(e$Standard_Deviation))
    }
  }
  expect_equal(devs, c(0.029317931, 0.014226835, 0.030171059, 0.026452308, 0.008968795, 0.040982844, 0.026119220, 0.008023552, 0.040535859, 0.026082652, 0.007801193, 0.040982844), tolerance = 1e-4)
})

test_that("Various CoxRegressionOmnibus options", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("rand", "fac", "dose")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)
  modelform <- "ME"

  cens_weight <- c(0)
  verbose <- FALSE
  devs <- c()
  options(warn = -1)
  modelform <- "ME"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  keep_constant <- c(1, 1, 1)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  lung_temp <- df$lung
  df$lung <- rep(0, length(lung_temp))
  keep_constant <- c(0, 0, 0)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  df$lung <- lung_temp
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1, 1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  a_n <- list(c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781))
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  #
  names <- c("rand", "fac", "dose")
  term_n <- c(0, 0, 1)
  tform <- c("lin", "lin", "lin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, -0.1, 0.2)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  options(warn = 0)
})

test_that("Various RunPoissonRegression_Omnibus options", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  time1 <- "t0"
  time2 <- "t1"
  df$pyr <- df$t1 - df$t0
  pyr <- "pyr"
  #
  event <- "lung"
  names <- c("rand", "fac", "dose")
  term_n <- c(0, 0, 1)
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2)
  modelform <- "ME"

  cens_weight <- c(0)
  verbose <- FALSE
  devs <- c()
  modelform <- "ME"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  keep_constant <- c(1, 1, 1)
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  lung_temp <- df$lung
  df$lung <- rep(0, length(lung_temp))
  keep_constant <- c(0, 0, 0)
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  df$lung <- lung_temp
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1, 1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  a_n <- list(c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781))
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  #
  names <- c("rand", "fac", "dose")
  term_n <- c(0, 0, 1)
  tform <- c("lin", "lin", "lin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, -0.1, 0.2)
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
})

test_that("Pois various_fixes", {
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
  tform <- c("loglin", "loglin", "loglin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(0.01, 0.1, 0.1)
  modelform <- "PAE"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  strat_col <- "fac"

  verbose <- FALSE
  model_control <- list("strata" = FALSE, "single" = FALSE)
  a_n <- c(0.01, 0.1, 0.1)
  modelform <- "PAE"
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  keep_constant <- c(1, 1, 1)
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  keep_constant <- c(0, 0, 0)
  ev0 <- df$lung
  df$lung <- rep(0, length(ev0))
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  names <- c("dose", "rand", "CONST")
  df$lung <- ev0
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
})
