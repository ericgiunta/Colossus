test_that("Coxph strata_basic_single_CR_null log_bound", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)
  df$weighting <- df$t1 / 20
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  expect_no_error(coxres <- CoxRun(Cox(tend = t1, event = lung) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  expect_no_error(coxres_s <- CoxRun(Cox_Strata(t0, t1, lung, rand) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  expect_no_error(coxres_c <- CoxRun(FineGray(t0, t1, lung, weighting) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  expect_no_error(coxres_sc <- CoxRun(FineGray_Strata(t0, t1, lung, rand, weighting) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  #
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(LikelihoodBound(coxres, df, curve_control, control = control, bad = FALSE))
  expect_error(LikelihoodBound(coxres, df, curve_control, control = control, norm = "bad"))
  expect_error(LikelihoodBound(coxres, df, curve_control = "bad"))
  #
  for (m in c(TRUE, FALSE)) {
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    curve_control <- list("manual" = m)
    expect_no_error(e <- LikelihoodBound(coxres, df, curve_control, control = control, bisect = FALSE))
    expect_no_error(LikelihoodBound(coxres_s, df, curve_control, control = control))
    expect_no_error(LikelihoodBound(coxres_c, df, curve_control, control = control))
    expect_no_error(LikelihoodBound(coxres_sc, df, curve_control, control = control))
    zz <- file(paste(tempfile(), ".txt", sep = ""), open = "wt")
    sink(zz)
    sink(zz, type = "message")
    print(e)
    sink(type = "message")
    sink(NULL)
    close(zz)
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "efron")
  }
  for (m in c(TRUE, FALSE)) {
    curve_control <- list("manual" = m)
    model_control <- list("null" = TRUE, "log_bound" = TRUE, "manual" = m)
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    coxres$modelcontrol$null <- TRUE
    expect_error(LikelihoodBound(coxres, df, curve_control, control = control))
    coxres$modelcontrol$null <- FALSE
    model_control <- list("single" = TRUE, "log_bound" = TRUE, "manual" = m)
    coxres$modelcontrol$single <- TRUE
    expect_error(LikelihoodBound(coxres, df, curve_control, control = control))
    coxres$modelcontrol$single <- FALSE
    model_control <- list("gradient" = TRUE, "log_bound" = TRUE, "manual" = m)
    coxres$modelcontrol$gradient <- TRUE
    expect_error(LikelihoodBound(coxres, df, curve_control, control = control))
    coxres$modelcontrol$gradient <- FALSE
  }
})
test_that("Poisson strata_single log_bound", {
  fname <- "ll_comp_0.csv"
  colTypes <- c("double", "double", "double", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  set.seed(3742)
  df$rand <- floor(runif(nrow(df), min = 0, max = 5))
  df$pyr <- df$t1 - df$t0
  keep_constant <- c(1, 0)
  a_n <- c(0, 0)

  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = 20, "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  #
  event <- "lung"
  a_n <- c(-0.1, -0.1)
  keep_constant <- c(0, 0)

  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  verbose <- FALSE
  expect_no_error(poisres <- PoisRun(Pois(pyr, lung) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  expect_no_error(poisres_s <- PoisRun(Pois_Strata(pyr, lung, rand) ~ loglinear(dose, fac, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  #
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_error(LikelihoodBound(poisres, df, curve_control, control = control, bad = FALSE))
  expect_error(LikelihoodBound(poisres, df, curve_control, control = control, norm = "bad"))
  expect_error(LikelihoodBound(poisres, df, curve_control = "bad"))
  #
  for (m in c(TRUE, FALSE)) {
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    curve_control <- list("manual" = m)
    expect_no_error(LikelihoodBound(poisres, df, curve_control, control = control, bisect = FALSE))
    expect_no_error(LikelihoodBound(poisres_s, df, curve_control, control = control))
  }
  for (m in c(TRUE, FALSE)) {
    model_control <- list("strata" = FALSE, "single" = TRUE, "log_bound" = TRUE, "manual" = m)
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    poisres$modelcontrol$single <- TRUE
    expect_error(LikelihoodBound(poisres, df, curve_control, control = control))
    poisres$modelcontrol$single <- FALSE
    model_control <- list("strata" = FALSE, "gradient" = TRUE, "log_bound" = TRUE, "manual" = m)
    a_n <- c(-0.1, -0.1)
    control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
    poisres$modelcontrol$gradient <- TRUE
    expect_error(LikelihoodBound(poisres, df, curve_control, control = control))
    poisres$modelcontrol$gradient <- FALSE
  }
})
test_that("Coxph EPICURE validated answers, loglin", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(0, 0)
  #
  a_n <- c(-0.6067, 5.019)
  model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(coxres <- CoxRun(Cox(entry, exit, event) ~ loglinear(dose0, dose1, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  coxres$beta_0 <- c(-0.6067, 5.019)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  v_lower <- c(-0.6305960, -0.6572672, -0.6817293, -0.6929630, -0.7300938, -0.7537744, -0.7749381, -0.8001031, -0.8175117)
  v_upper <- c(-0.5828725, -0.5562505, -0.5318645, -0.5206756, -0.4837373, -0.4602148, -0.4392159, -0.4142752, -0.3970399)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  alphas <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    a_n <- c(-0.6067, 5.019)
    model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = alphas[alpha_i], "para_number" = 1)
    curve_control <- list("alpha" = alphas[alpha_i], "para_number" = 1)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }

  v_lower <- c(4.981497, 4.939337, 4.900838, 4.883211, 4.825191, 4.788380, 4.755608, 4.716794, 4.690041)
  v_upper <- c(5.057414, 5.100032, 5.139239, 5.157283, 5.217094, 5.255376, 5.289680, 5.330581, 5.358945)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    a_n <- c(-0.6067, 5.019)
    model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = alphas[alpha_i], "para_number" = 2)
    curve_control <- list("alpha" = alphas[alpha_i], "para_number" = 2)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }
})

test_that("Coxph EPICURE validated answers, loglin manual", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(-0.6067, 5.019)
  model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(coxres <- CoxRun(Cox(entry, exit, event) ~ loglinear(dose0, dose1, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  coxres$beta_0 <- c(-0.6067, 5.019)
  v_lower <- c(-0.6305960, -0.6572672, -0.6817293, -0.6929630, -0.7300938, -0.7537744, -0.7749381, -0.8001031, -0.8175117)
  v_upper <- c(-0.5828725, -0.5562505, -0.5318645, -0.5206756, -0.4837373, -0.4602148, -0.4392159, -0.4142752, -0.3970399)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 2, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  alphas <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    a_n <- c(-0.6067, 5.019)
    model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = alphas[alpha_i], "para_number" = 1, "manual" = TRUE)
    curve_control <- list("alpha" = alphas[alpha_i], "para_number" = 1, "manual" = TRUE)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }

  v_lower <- c(4.981497, 4.939337, 4.900838, 4.883211, 4.825191, 4.788380, 4.755608, 4.716794, 4.690041)
  v_upper <- c(5.057414, 5.100032, 5.139239, 5.157283, 5.217094, 5.255376, 5.289680, 5.330581, 5.358945)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    a_n <- c(-0.6067, 5.019)
    model_control <- list("basic" = TRUE, "log_bound" = TRUE, "alpha" = alphas[alpha_i], "para_number" = 2, "manual" = TRUE)
    curve_control <- list("alpha" = alphas[alpha_i], "para_number" = 2, "manual" = TRUE)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }
})

test_that("Coxph, lin both", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0, 0)
  # a_n <- c(0.2462, 5.020, -0.5909)
  a_n <- c(-1.493177, 5.020007, 1.438377)
  model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = FALSE, "alpha" = 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(-1, -1), "halfmax" = -1, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  coxres <- CoxRun(Cox(entry, exit, event) ~ loglinear(dose0, dose1, 0) + linear(dose0, 1) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant)
  coxres$beta_0 <- c(-1.493177, 5.020007, 1.438377)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(2, 2), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  alpha <- 0.005
  a_n <- c(-1.493177, 5.020007, 1.438377)
  model_control <- list("basic" = FALSE, "maxstep" = 2, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = FALSE)
  curve_control <- list("alpha" = alpha, "para_number" = 2, "manual" = FALSE, "maxstep" = 2)
  expect_no_error(LikelihoodBound(coxres, df, curve_control, control = control))
  a_n <- c(-1.493177, 5.020007, 1.438377)
  model_control <- list("basic" = FALSE, "maxstep" = 2, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = TRUE)
  curve_control <- list("alpha" = alpha, "para_number" = 2, "manual" = TRUE, "maxstep" = 2)
  expect_no_error(LikelihoodBound(coxres, df, curve_control, control = control))
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  alpha_list <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(10, 10), "halfmax" = 5, "epsilon" = 1e-4, "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  v_lower <- c(4.98213, 4.939947, 4.901443, 4.883815, 4.825766, 4.788947, 4.756168, 4.717337, 4.690578)
  # c(4.97252283668956, 4.9349945105648, 4.89804715665926, 4.88084912208962, 4.82369762341988, 4.78721237571926, 4.7546530342797, 4.71603055250556, 4.68938287303871)
  v_upper <- c(5.058015, 5.100656, 5.139868, 5.157913, 5.217754, 5.256045, 5.290357, 5.331277, 5.359649)
  # c(5.06762896572498, 5.10561529697034, 5.14327069556976, 5.16088604918614, 5.21982880792394, 5.2577860471215, 5.29187760184654, 5.33258757872226, 5.36084782852899)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    alpha <- alpha_list[alpha_i]
    a_n <- c(-1.493177, 5.020007, 1.438377)
    model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = TRUE)
    curve_control <- list("alpha" = alpha, "para_number" = 2, "manual" = TRUE, "maxstep" = 100)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }

  v_lower <- c(1.174142, 0.8706736, 0.565792, 0.4011414, -0.2129912, 0.2830731, -0.07262359, -0.1306881, 0.3065214)
  # c(-0.643365949558998, -0.677336655540846, -0.706075250211414, -0.718165409196492, -0.753647332793819, -0.773208334303991, -0.789018704115451, -0.806061085000755, -0.816875114954096)
  v_upper <- c(1.702503, 2.0018595, 2.283068, 2.4149146, 2.8653084, 3.1659745, 3.44475771, 3.7898521, 4.0379529)
  # c(-0.521472203247917, -0.444964438813732, -0.327862977142017, -0.235044092073815, 2.91573713669059, 3.21014641617297, 3.48490803194128, 3.82648584413642, 4.07272009904963)
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    alpha <- alpha_list[alpha_i]
    a_n <- c(-1.493177, 5.020007, 1.438377)
    model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 3, "manual" = TRUE)
    curve_control <- list("alpha" = alpha, "para_number" = 3, "manual" = TRUE, "maxstep" = 100)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    # expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-4)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-4)
  }
})

test_that("Poisson, lin both", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(-2.917, 0.06526)
  #
  model_control <- list("basic" = FALSE, "maxstep" = 10, "log_bound" = FALSE, "alpha" = 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  poisres <- PoisRun(Pois(exit, event) ~ loglinear(dose0, 0) + linear(dose1, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant)
  poisres$beta_0 <- c(-2.917, 0.06526)
  alpha <- 0.005
  a_n <- c(-2.917, 0.06526)
  model_control <- list("basic" = FALSE, "maxstep" = 3, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = FALSE)
  curve_control <- list("alpha" = alpha, "para_number" = 2, "manual" = FALSE, "maxstep" = 3)
  expect_no_error(LikelihoodBound(poisres, df, curve_control, control = control))
  a_n <- c(-2.917, 0.06526)
  model_control <- list("basic" = FALSE, "maxstep" = 3, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = TRUE)
  curve_control <- list("alpha" = alpha, "para_number" = 2, "manual" = TRUE, "maxstep" = 3)
  expect_no_error(LikelihoodBound(poisres, df, curve_control, control = control))
})

test_that("Coxph, lin both, curve search", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0, 0)
  # a_n <- c(0.2462, 5.020, -0.5909)
  a_n <- c(-1.493177, 5.020007, 1.438377)
  model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = FALSE, "alpha" = 0.005)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 1, "epsilon" = 1e-3, "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  expect_no_error(coxres <- CoxRun(Cox(entry, exit, event) ~ loglinear(dose0, dose1, 0) + linear(dose0, 1) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  coxres$beta_0 <- c(-1.493177, 5.020007, 1.438377)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(10, 10), "halfmax" = 5, "epsilon" = 1e-3, "deriv_epsilon" = 1e-3, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow", "guesses" = 10)
  #
  alpha <- 0.005
  a_n <- c(-1.493177, 5.020007, 1.438377)

  alpha_list <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(10, 10), "halfmax" = 2, "epsilon" = 1e-4, "deriv_epsilon" = 1e-4, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")

  model_control <- list("basic" = FALSE, "maxstep" = 10, "log_bound" = FALSE, "alpha" = 0.005, "para_number" = 3, "step_size" = 0.5)
  v_lower <- c(1.174125, 0.8706585, 0.5657879, 0.401115, -0.749947, -0.7704548, -0.7868122, -0.8043293, -0.8154377)
  v_upper <- c(1.702507, 2.0018841, 2.283073, 2.414909, 2.865287, 3.1660076, 3.4447552, 3.7898479, 4.0379559)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    alpha <- alpha_list[alpha_i]
    a_n <- c(-1.493177, 5.020007, 1.438377)
    model_control <- list("basic" = FALSE, "maxstep" = 20, "log_bound" = FALSE, "alpha" = alpha, "para_number" = 3, "step_size" = 0.5)
    curve_control <- list("alpha" = alpha, "para_number" = 3, "bisect" = TRUE, "step_size" = 0.5, "maxstep" = 20, "manual" = FALSE)
    e <- LikelihoodBound(coxres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-3)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-3)
  }
})

test_that("Poisson, curve search", {
  fname <- "base_example.csv"
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE)
  keep_constant <- c(0, 0)
  a_n <- c(-2.917, 0.06526)
  model_control <- list("basic" = FALSE, "maxstep" = 100, "log_bound" = FALSE, "alpha" = 0.1)
  control <- list("ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1), "halfmax" = 1, "epsilon" = 1e-4, "deriv_epsilon" = 1e-4, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  expect_no_error(poisres <- PoisRun(Pois(exit, event) ~ loglinear(dose0, 0) + plinear(dose1, 0) + multiplicative - excess(), df, a_n = a_n, control = control, keep_constant = keep_constant))
  poisres$beta_0 <- c(-2.917, 0.06526)
  control <- list("ncores" = 2, "lr" = 0.75, "maxiters" = c(10, 10), "halfmax" = 5, "epsilon" = 1e-4, "deriv_epsilon" = 1e-4, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  alpha_list <- c(0.75, 0.5, 1 - 0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
  v_lower <- c(0.05418212, 0.04221923, 0.03135497, 0.02635009, 0.0100537, -0.0002002051, -0.009355479, -0.02009767, -0.02748292)
  v_upper <- c(0.07591063, 0.0881787, 0.09953124, 0.10471923, 0.1221753, 0.1333447168, 0.143415518, 0.15550048, 0.16386229)
  if (!isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))) {
    skip("Cran Skip")
  }
  for (alpha_i in c(1, 5)) { # seq_along(alphas)) {
    alpha <- alpha_list[alpha_i]
    a_n <- c(-2.917, 0.06526)
    model_control <- list("basic" = FALSE, "maxstep" = 20, "log_bound" = TRUE, "alpha" = alpha, "para_number" = 2, "manual" = FALSE)
    curve_control <- list("alpha" = alpha, "para_number" = 2, "bisect" = TRUE, "step_size" = 0.5, "maxstep" = 20, "manual" = FALSE)
    e <- LikelihoodBound(poisres, df, curve_control, control = control)
    a <- e$Parameter_Limits
    expect_equal(a[1], v_lower[alpha_i], tolerance = 1e-3)
    expect_equal(a[2], v_upper[alpha_i], tolerance = 1e-3)
  }
})
