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
  modelform <- "M"

  verbose <- FALSE

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "null" = FALSE)
  e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "M", control = control, strat_col = "fac", model_control = model_control)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (j in c(TRUE, FALSE)) {
    for (k in c(TRUE, FALSE)) {
      for (l in c(TRUE, FALSE)) {
        model_control <- list("strata" = FALSE, "basic" = j, "single" = k, "null" = l)
        a_n <- c(0.0)
        e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "M", control = control, strat_col = "fac", model_control = model_control)
        expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-2)
      }
    }
  }
  model_control <- list("strata" = TRUE, "basic" = FALSE, "single" = FALSE, "null" = FALSE)
  e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "M", control = control, strat_col = "fac", model_control = model_control)
  for (j in c(TRUE, FALSE)) {
    for (k in c(TRUE, FALSE)) {
      for (l in c(TRUE, FALSE)) {
        model_control <- list("strata" = TRUE, "basic" = j, "single" = k, "null" = l)
        a_n <- c(0.0)
        e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = c("loglin"), keep_constant = keep_constant, a_n = a_n, modelform = "M", control = control, strat_col = "fac", model_control = model_control)
        expect_equal(e0$LogLik, e1$LogLik, tolerance = 1e-2)
      }
    }
  }
})
test_that("Coxph strata_basic_single_CR", {
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  plot_options <- list("name" = paste(tempfile(), "run", sep = ""), "verbose" = FALSE, "studyid" = "studyid", "age_unit" = "years")
  dft <- GetCensWeight(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, plot_options)
  #
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE
  j_iterate <- 1
  LL_comp <- c(-81.12079, -81.12079, -77.97632, -77.97632, -68.60263, -68.67108, -75.34028, -75.3691, -81.12079, -81.12079, -77.97632, -77.97632, -68.60263, -68.67108, -75.34028, -75.3691, -122.8909, -122.8909, -119.9814, -119.9814, -109.6211, -109.6742, -117.0147, -117.0539, -122.8909, -122.8909, -119.9814, -119.9814, -109.6211, -109.6742, -117.0147, -117.0539)
  dev_comp <- c(0, 0, 0, 0, 0.6850986, 0.6804879, 0.7627068, 0.7600753, 0, 0, 0, 0, 0.6850986, 0.6804879, 0.7627068, 0.7600753, 0, 0, 0, 0, 0.6720932, 0.6688576, 0.7401029, 0.7389506, 0, 0, 0, 0, 0.6720932, 0.6688576, 0.7401029, 0.7389506)
  for (i in c(TRUE, FALSE)) {
    for (j in c(TRUE, FALSE)) {
      for (k in c(TRUE, FALSE)) {
        for (l in c(TRUE, FALSE)) {
          model_control <- list("strata" = i, "basic" = j, "single" = k, "cr" = l)
          a_n <- c(-0.1, -0.1)
          control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
          e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
          expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-3)
          if (k == FALSE) {
            expect_equal(sum(e$Standard_Deviation), dev_comp[j_iterate], tolerance = 1e-3)
          }
          j_iterate <- j_iterate + 1
          a_n <- c(-0.1, -0.1)
          control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "efron", "double_step" = 0)
          e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "weighting")
          expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-3)
          if (k == FALSE) {
            expect_equal(sum(e$Standard_Deviation), dev_comp[j_iterate], tolerance = 1e-3)
          }
          j_iterate <- j_iterate + 1
        }
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 0)
  strat_col <- "fac"
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  verbose <- FALSE
  j_iterate <- 1
  LL_comp <- c(-468.7465, -464.9279, -461.2769, -462.1182, -3033.332, -2734.64, -992.622, -1334.36)
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 0)
  strat_col <- "fac"

  verbose <- FALSE
  modelforms <- c("A", "PAE", "M", "PA")
  j_iterate <- 1
  LL_comp <- c(-790.3591, -471.0312, -471.0312, -463.1375, -707.56, -678.2228, -678.2228, -471.4805)
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
  names <- c("dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand", "rand")
  term_n <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  tform <- c("loglin_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope")
  keep_constant <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  a_n <- c(1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)

  modelform <- "PAE"

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 0)
  strat_col <- "fac"
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  verbose <- FALSE
  j_iterate <- 1
  LL_comp <- c(-496.736611833105, -475.421339181129, -466.96134309201, -463.218078300204, -4497.1782941901, -3577.95315136081, -2557.93228186619, -2350.19459732101)
  for (i in c(TRUE, FALSE)) {
    for (j in c(TRUE, FALSE)) {
      model_control <- list("strata" = i, "single" = j)
      a_n <- c(1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)
      modelform <- "PAE"
      e <- RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
      expect_equal(e$LogLik, LL_comp[j_iterate], tolerance = 1e-2)
      j_iterate <- j_iterate + 1
      a_n <- c(1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)
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
  df$censor <- (df$lung == 0)
  #
  event <- "lung"
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 1)
  tform <- c("loglin", "lin", "lin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  modelform <- "M"


  cens_weight <- c(0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE

  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)

  model_list <- c("M", "A", "PA", "PAE")
  means <- c(0.053270537, 0.053270537, 0.053270537, 0.053270537, 0.053270537, 0.223545566, 0.223545566, 0.223545566, 0.223545566, 0.223545566, 0.007932817, 0.007932817, 0.007932817, 0.007932817, 0.007932817, 0.053270537, 0.053270537, 0.053270537, 0.053270537, 0.053270537)
  medians <- c(0.08229779, 0.08229779, 0.08229779, 0.08229779, 0.08229779, 0.15570474, 0.15570474, 0.15570474, 0.15570474, 0.15570474, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.08229779, 0.08229779, 0.08229779, 0.08229779, 0.08229779)
  sums <- c(53.270537, 53.270537, 53.270537, 53.270537, 53.270537, 223.545566, 223.545566, 223.545566, 223.545566, 223.545566, 7.932817, 7.932817, 7.932817, 7.932817, 7.932817, 53.270537, 53.270537, 53.270537, 53.270537, 53.270537)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (model_i in 1:4) {
    modelform <- model_list[model_i]
    for (i in 1:5) {
      model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
      #
      df_order$order <- sample(df_order$order)
      setorderv(df_order, c("order"))
      term_n <- df_order$term_n
      tform <- df_order$tform
      keep_constant <- df_order$keep_constant
      a_n <- df_order$a_n
      names <- df_order$names
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
      e <- Cox_Relative_Risk(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)$Risk

      a_i <- (model_i - 1) * 5 + i
      expect_equal(mean(e), means[a_i], tolerance = 1e-2)
      expect_equal(median(e), medians[a_i], tolerance = 1e-2)
      expect_equal(sum(e), sums[a_i], tolerance = 1e-2)
    }
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
  df$censor <- (df$lung == 0)
  #
  event <- "lung"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE

  model_list <- c("GMIX-R", "GMIX-E", "GMIX")
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 2)
  tform <- c("loglin", "loglin", "plin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)

  means <- c(0.6918334, 2.9080098, 0.6918334, 1.2728663, 1.8214142, 2.9080098)
  medians <- c(0.5871897, 2.8144283, 0.5871897, 1.0226682, 1.7472830, 2.8144283)
  sums <- c(691.8334, 2908.0098, 691.8334, 1272.8663, 1821.4142, 2908.0098)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  count <- 0
  for (model_i in 1:3) {
    modelform <- model_list[model_i]
    if (modelform == "GMIX") {
      for (term_i in 0:3) {
        model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
        model_control$gmix_term <- c(0, term_i %% 2, floor(term_i / 2))
        #
        df_order$order <- sample(df_order$order)
        setorderv(df_order, c("order"))
        term_n <- df_order$term_n
        tform <- df_order$tform
        keep_constant <- df_order$keep_constant
        a_n <- df_order$a_n
        names <- df_order$names
        #
        control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
        e <- Cox_Relative_Risk(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)$Risk
        count <- count + 1
        expect_equal(mean(e), means[count], tolerance = 1e-2)
        expect_equal(median(e), medians[count], tolerance = 1e-2)
        expect_equal(sum(e), sums[count], tolerance = 1e-2)
      }
    } else {
      model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
      #
      df_order$order <- sample(df_order$order)
      setorderv(df_order, c("order"))
      term_n <- df_order$term_n
      tform <- df_order$tform
      keep_constant <- df_order$keep_constant
      a_n <- df_order$a_n
      names <- df_order$names
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
      e <- Cox_Relative_Risk(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)$Risk
      count <- count + 1
      expect_equal(mean(e), means[count], tolerance = 1e-2)
      expect_equal(median(e), medians[count], tolerance = 1e-2)
      expect_equal(sum(e), sums[count], tolerance = 1e-2)
    }
  }
})

test_that("gmix omnibus use", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  df$censor <- (df$lung == 0)
  #
  event <- "lung"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE

  model_list <- c("GMIX-R", "GMIX-E", "GMIX")
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 2)
  tform <- c("loglin", "loglin", "plin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  options(warn = -1)
  count <- 0
  cens_weight <- c(0)
  for (model_i in 1:3) {
    modelform <- model_list[model_i]
    if (modelform == "GMIX") {
      for (term_i in 0:3) {
        model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
        model_control$gmix_term <- c(0, term_i %% 2, floor(term_i / 2))
        #
        df_order$order <- sample(df_order$order)
        setorderv(df_order, c("order"))
        term_n <- df_order$term_n
        tform <- df_order$tform
        keep_constant <- df_order$keep_constant
        a_n <- df_order$a_n
        names <- df_order$names
        #
        control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
        expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "null"))
      }
    } else {
      model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
      #
      df_order$order <- sample(df_order$order)
      setorderv(df_order, c("order"))
      term_n <- df_order$term_n
      tform <- df_order$tform
      keep_constant <- df_order$keep_constant
      a_n <- df_order$a_n
      names <- df_order$names
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
      expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "rand", model_control = model_control, cens_weight = "null"))
    }
  }
})

test_that("risk check omnibus dose", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))

  time1 <- "t0"
  time2 <- "t1"
  df$censor <- (df$lung == 0)
  #
  event <- "lung"
  names <- c("dose", "fac", "dose", "fac", "rand")
  term_n <- c(0, 0, 1, 1, 1)
  tform <- c("loglin", "lin", "lin", "plin", "loglin")
  keep_constant <- c(0, 0, 0, 0, 0)
  a_n <- c(-0.1, 0.1, 0.2, 0.3, -0.5)
  modelform <- "M"


  cens_weight <- c(0)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)

  verbose <- FALSE

  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:5)

  model_list <- c("M", "A", "PA", "PAE")
  names <- c("dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "dose", "fac", "fac", "fac", "fac", "fac", "fac", "fac", "fac", "fac", "fac", "fac")
  term_n <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  tform <- c("loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope", "loglin_top", "lin_slope", "lin_int", "quad_slope", "step_slope", "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope")
  keep_constant <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  a_n <- c(-0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1, -0.1, -0.1, 1, -0.1, 1, 2, 0.3, 1.5, 0.2, 0.7, 1)

  df_order <- data.table("term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n, "names" = names, "order" = 1:22)
  means <- c(3.268830, 2.672975, 1.753389, 3.268830)
  medians <- c(2.898413, 2.558679, 1.636409, 2.898413)
  sums <- c(3268.830, 2672.975, 1753.389, 3268.830)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (model_i in 1:4) {
    modelform <- model_list[model_i]
    for (i in 1:5) {
      model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
      #
      df_order$order <- sample(df_order$order)
      setorderv(df_order, c("order"))
      term_n <- df_order$term_n
      tform <- df_order$tform
      keep_constant <- df_order$keep_constant
      a_n <- df_order$a_n
      names <- df_order$names
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
      e <- Cox_Relative_Risk(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)$Risk
      expect_equal(mean(e), means[model_i], tolerance = 1e-2)
      expect_equal(median(e), medians[model_i], tolerance = 1e-2)
      expect_equal(sum(e), sums[model_i], tolerance = 1e-2)
    }
  }
})

test_that("Coxph relative risk combinations", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  time1 <- "a"
  time2 <- "b"
  event <- "c"
  names <- c("d")
  term_n <- c(0)
  tform <- c("loglin")
  keep_constant <- c(0)
  a_n <- c(-0.1)
  modelform <- "M"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  keep_constant <- c(1)
  expect_error(Cox_Relative_Risk(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control))
  keep_constant <- c(0)
  names <- c("d", "CONST")
  term_n <- c(0, 0)
  tform <- c("loglin", "loglin")
  keep_constant <- c(0, 0)
  a_n <- c(-0.1, 0.1)
  expect_no_error(Cox_Relative_Risk(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control))
})

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
  modelform <- "M"


  cens_weight <- c(0)

  verbose <- FALSE

  devs <- c()

  modelform <- "M"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  for (i in 1:3) {
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6428582, 0.4240752, 0.1507817)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
  modelform <- "M"


  cens_weight <- c(0)

  verbose <- FALSE

  devs <- c()

  modelform <- "M"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE, "oberved_info" = TRUE)
  for (i in 1:3) {
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0, 0, 0)
    keep_constant[i] <- 1
    #
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
  }
  a_n <- c(0.6428582, 0.4240752, 0.1507817)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
    modelform <- "M"
    control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col = "fac", model_control = model_control, cons_mat = Constraint_Matrix, cons_vec = Constraint_const)
    expect_equal(e$beta_0, c(0.357333, 0.357333), tolerance = 1e-2)
  }
  for (i in 1:5) {
    a_n <- 2 * runif(2) - 1
    del <- abs(a_n[1] - a_n[2])
    a_n0 <- rep(sum(a_n) / 2, 2)
    a_n <- a_n0 + c(-del / 2, del / 2)
    modelform <- "M"
    control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
  modelform <- "M"
  devs <- c()

  modelform <- "M"
  for (inma_type in c(T, F)) {
    model_control <- list("oberved_info" = inma_type)
    for (i in 1:3) {
      a_n <- c(0.6465390, 0.4260961, 0.1572781)
      keep_constant <- c(0, 0, 0)
      keep_constant[i] <- 1
      #
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
      control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
  modelform <- "M"

  cens_weight <- c(0)
  verbose <- FALSE
  devs <- c()
  options(warn = -1)
  modelform <- "M"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  keep_constant <- c(1, 1, 1)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  lung_temp <- df$lung
  df$lung <- rep(0, length(lung_temp))
  keep_constant <- c(0, 0, 0)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  df$lung <- lung_temp
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1, 1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  a_n <- list(c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781))
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1, "guesses" = 1)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1, "guesses" = 10)
  expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
  #
  names <- c("rand", "fac", "dose")
  term_n <- c(0, 0, 1)
  tform <- c("lin", "lin", "lin")
  keep_constant <- c(0, 0, 0)
  a_n <- c(-0.1, -0.1, 0.2)
  expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "fac", model_control = model_control))
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
  modelform <- "M"

  cens_weight <- c(0)
  verbose <- FALSE
  devs <- c()
  modelform <- "M"
  model_control <- list("strata" = FALSE, "basic" = FALSE, "single" = FALSE, "cr" = FALSE)
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  keep_constant <- c(0, 0, 0)
  #
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1, 1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1)
  a_n <- list(c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781), c(0.6465390, 0.4260961, 0.1572781))
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  a_n <- c(0.6465390, 0.4260961, 0.1572781)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1, "guesses" = 50)
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  #
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1, "guesses" = 1)
  expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
  control <- list("ncores" = 2, "lr" = 0.75, "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 1, "guesses" = 10)
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

  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 0, "ties" = "breslow", "double_step" = 0)
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
  control$guesses <- 100
  expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control))
})
