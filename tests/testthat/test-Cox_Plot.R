######################
## PLOTTING CHECKS
######################

# test_that("Coxph plot time column missing", {
#   a <- c(0, 1, 2, 3, 4, 5, 6)
#   b <- c(1, 2, 3, 4, 5, 6, 7)
#   c <- c(0, 1, 0, 0, 0, 1, 0)
#   d <- c(3, 4, 5, 6, 7, 8, 9)
#   df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
#   time1 <- "a_bad"
#   time2 <- "b"
#   event <- "c"
#   names <- c("d")
#   term_n <- c(0)
#   tform <- c("loglin")
#   keep_constant <- c(0)
#   a_n <- c(-0.1)
#   modelform <- "M"
#   control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
#   plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = FALSE, "surv_curv" = FALSE, "studyid" = "a", "verbose" = 0)
#   if (system.file(package = "ggplot2") != "") {
#     expect_error(RunCoxPlots(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, plot_options))
#   }
# })
# test_that("Coxph plot no events", {
#   a <- c(0, 1, 2, 3, 4, 5, 6)
#   b <- c(1, 2, 3, 4, 5, 6, 7)
#   c <- c(0, 0, 0, 0, 0, 0, 0)
#   d <- c(3, 4, 5, 6, 7, 8, 9)
#   df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
#   time1 <- "a"
#   time2 <- "b"
#   event <- "c"
#   names <- c("d")
#   term_n <- c(0)
#   tform <- c("loglin")
#   keep_constant <- c(0)
#   a_n <- c(-0.1)
#   modelform <- "M"
#   control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
#   plot_options <- list("type" = c("surv", "run"), "martingale" = FALSE, "surv_curv" = FALSE, "studyid" = "a", "verbose" = 0)
#   if (system.file(package = "ggplot2") != "") {
#     expect_error(RunCoxPlots(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control, plot_options))
#   }
# })
test_that("Coxph plot no type", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  plot_options <- list("martingale" = FALSE, "surv_curv" = FALSE, "studyid" = "a", "verbose" = 0)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot strata col not in data", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  plot_options <- list("type" = c("surv", "run"), "martingale" = FALSE, "surv_curv" = FALSE, "strat_haz" = TRUE, "strat_col" = "e", "studyid" = "a", "verbose" = 0)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot strata col not given", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  plot_options <- list("type" = c("surv", "run"), "martingale" = FALSE, "surv_curv" = FALSE, "strat_haz" = TRUE, "studyid" = "a", "verbose" = 0)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot martingale dose col not in data", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", "run"), "martingale" = TRUE, "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot martingale dose col not given", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", "run"), "martingale" = TRUE, "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot Survival ID col not in data", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", "run"), "martingale" = FALSE, "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "f", "verbose" = 0, "km" = TRUE)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})
test_that("Coxph plot Survival ID col not given", {
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", "run"), "martingale" = FALSE, "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "verbose" = 0, "km" = TRUE)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
})

test_that("Coxph risk too few unique values", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 0, 1, 0, 1, 0, 0)
  d <- c(3, 1, 3, 1, 2, 3, 1)
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("RISK", paste(tempfile(), "run", sep = "")), "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  d <- c(3, 3, 3, 3, 3, 3, 3)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  options(warn = -1)
  if (system.file(package = "ggplot2") != "") {
    expect_error(plot(e, df, plot_options))
  }
  options(warn = 0)
})
test_that("Coxph risk plotting above discrete step number limit", {
  a <- rep(c(0, 1, 2, 3, 4, 5, 6), 20)
  b <- rep(c(1, 2, 3, 4, 5, 6, 7), 20)
  c <- rep(c(1, 0, 1, 0, 1, 0, 0), 20)
  set.seed(3742)
  d <- runif(length(a))
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("RISK", paste(tempfile(), "run", sep = "")), "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(plot(e, df, type = c("RISK", paste(tempfile(), "run", sep = "")), studyid = "a", verbose = 0))
    expect_error(plot(e, df, type = c("RISK", paste(tempfile(), "run", sep = "")), studyid = "a", bad = 0))
    expect_error(plot(e, df, plot_options = c(0))$d)
    expect_no_error(ep <- plot(e, df, plot_options)$d)
    expect_equal(ep$y[100], 0.9756069, tolerance = 1e-4)
  }
})

test_that("Coxph plot no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 0, 1, 0, 1, 0, 0)
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "surv_curv" = TRUE, "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(tend = b, event = c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options, verbose = 0)$standard)
    expect_equal(ep$h, c(0.000000, 4.344177, 18.450904), tolerance = 1e-4)
  }
})
test_that("Coxph plot stratafied no error", {
  a <- c(0, 0, 0, 0, 0, 0, 1)
  b <- c(2, 3, 4, 2, 3, 2, 4)
  c <- c(1, 0, 1, 0, 1, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 0, 0, 1, 0, 1)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  time1 <- "a"
  time2 <- "b"
  event <- "c"
  names <- c("d")
  term_n <- c(0)
  tform <- c("loglin")
  keep_constant <- c(0)
  a_n <- c(-0.1)
  modelform <- "M"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "surv_curv" = TRUE, "strat_haz" = TRUE, "strat_col" = "e", "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox_Strata(a, b, c, e) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options)$stratified_survival)
    expect_equal(ep$surv[2], 0.5880531, tolerance = 1e-4)
  }
})

test_that("Coxph risk no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 0, 1, 0, 1, 0, 0)
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("RISK", paste(tempfile(), "run", sep = "")), "studyid" = "a", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options)$d)
    expect_equal(ep$y[7], 0.5488116, tolerance = 1e-4)
  }
})

test_that("Coxph schoenfeld no error", {
  fname <- "MULTI_COV.csv"
  colTypes <- c("double", "double", "integer", "integer", "integer")
  df <- fread(fname, nThread = min(c(detectCores(), 2)), data.table = TRUE, header = TRUE, colClasses = colTypes, verbose = FALSE, fill = TRUE)
  time1 <- "t0"
  time2 <- "t1"
  event <- "lung"
  names <- c("a", "b")
  term_n <- c(0, 1)
  tform <- c("loglin", "loglin")
  keep_constant <- c(0, 0)
  a_n <- c(0.01, -15)
  modelform <- "M"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("SCHOENFELD", paste(tempfile(), "run", sep = "")), "studyid" = "t0", "verbose" = 0)
  e <- CoxRun(Cox(t0, t1, lung) ~ loglinear(a, b, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options)$b)
    expect_equal(ep$y[16], 9.999993e-01, tolerance = 1e-4)
  }
})
test_that("Coxph martingale no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 2, 2, 3, 3, 3)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e)
  time1 <- "a"
  time2 <- "b"
  event <- "c"
  names <- c("d")
  term_n <- c(0)
  tform <- c("loglin")
  keep_constant <- c(0)
  a_n <- c(-0.1)
  modelform <- "M"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = TRUE, "cov_cols" = "d", "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "e", "verbose" = 0)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options)$d)
    expect_equal(ep$res_sum[2], -0.2262094, tolerance = 1e-4)
  }
})
test_that("Coxph martingale combinations", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 1, 0, 0, 0, 1, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  e <- c(1, 1, 2, 2, 3, 3, 3)
  f <- c(1, 2, 3, 4, 5, 6, 7)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = e, "f" = f)
  time1 <- "a"
  time2 <- "b"
  event <- "c"
  names <- c("d")
  term_n <- c(0)
  tform <- c("loglin")
  keep_constant <- c(0)
  a_n <- c(-0.1)
  modelform <- "M"
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = TRUE, "cov_cols" = "d", "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "Not_In", "verbose" = 0)
  keep_constant <- c(1)
  options(warn = -1)
  if (system.file(package = "ggplot2") != "") {
    e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
    expect_error(plot(e, df, plot_options))
    names <- c("d", "CONST")
    term_n <- c(0, 0)
    tform <- c("loglin", "loglin")
    keep_constant <- c(0, 0)
    a_n <- c(-0.1, 0.1)
    plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = TRUE, "cov_cols" = "d", "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "f", "verbose" = 0)
    e <- CoxRun(Cox(a, b, c) ~ loglinear(d, CONST, 0), df, control = control, a_n = a_n)
    expect_no_error(plot(e, df, plot_options))
    plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = TRUE, "cov_cols" = "Not_In", "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "f", "verbose" = 0)
    e <- CoxRun(Cox(a, b, c) ~ loglinear(d, CONST, 0), df, control = control, a_n = a_n)
    expect_error(plot(e, df, plot_options))
    plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "martingale" = TRUE, "cov_cols" = "d", "surv_curv" = FALSE, "strat_haz" = FALSE, "smooth_haz" = FALSE, "studyid" = "f", "verbose" = 0)
    e <- CoxRun(Cox(a, b, c) ~ loglinear(d, CONST, 0), df, control = control, a_n = a_n)
    df$c <- rep(0, nrow(df))
    expect_error(plot(e, df, plot_options))
  }
  options(warn = 0)
})
test_that("Coxph km no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 0, 1, 0, 1, 0, 0)
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
  control <- list("ncores" = 2, "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  plot_options <- list("type" = c("surv", paste(tempfile(), "run", sep = "")), "surv_curv" = TRUE, "studyid" = "a", "verbose" = 0, "km" = TRUE)
  e <- CoxRun(Cox(a, b, c) ~ loglinear(d, 0), df, control = control, a_n = a_n)
  if (system.file(package = "ggplot2") != "") {
    expect_no_error(ep <- plot(e, df, plot_options)$`kaplin-meier`)
    expect_equal(ep$n_t[3], 0.7142857, tolerance = 1e-4)
  }
})
