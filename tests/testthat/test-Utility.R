## ------------------------------------- ##
## Verify the system check code
## ------------------------------------- ##
test_that("System version", {
  expect_no_error(System_Version())
  expect_no_error(Rcpp_version())
})

## ------------------------------------- ##
## Verify verbosity code
## ------------------------------------- ##
test_that("Check verbose, T/F", {
  expect_no_error(Check_Verbose(TRUE))
  expect_no_error(Check_Verbose(FALSE))
  expect_no_error(Check_Verbose("TRUE"))
})
test_that("Check verbose, 0/1/2/3/4", {
  expect_no_error(Check_Verbose(0))
  expect_no_error(Check_Verbose(1))
  expect_no_error(Check_Verbose(2))
  expect_no_error(Check_Verbose(3))
  expect_no_error(Check_Verbose(4))
  expect_no_error(Check_Verbose("1"))
})
test_that("Check verbose, Fails", {
  expect_error(Check_Verbose(-1))
  expect_error(Check_Verbose(5))
  expect_error(Check_Verbose("true"))
})

## ------------------------------------- ##
## Duplicate Columns
## ------------------------------------- ##

test_that("No dupe columns", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "d"), c(0, 0, 0, 0), TRUE), c("a", "b", "c", "d"))
})
test_that("No columns", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_equal(Check_Dupe_Columns(df, c(), c(), TRUE), c())
})
test_that("One column with varying", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_equal(Check_Dupe_Columns(df, c("a"), c(0), TRUE), c("a"))
})
test_that("One column with constant", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_equal(Check_Dupe_Columns(df, c("c"), c(0), TRUE), c("c"))
})
test_that("One column with constant 0", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_equal(Check_Dupe_Columns(df, c("c"), c(0), TRUE), c())
})
test_that("One duplicate column", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = a)
  options(warn = -1)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "d", "e"), c(0, 0, 0, 0, 0), TRUE), c("a", "b", "c", "d"))
  options(warn = 0)
})
test_that("One duplicate column, different term", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = a)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "d", "e"), c(0, 0, 0, 1, 1), TRUE), c("a", "b", "c", "d", "e"))
})
test_that("Multiple duplicate columns", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d, "e" = a, "f" = b)
  options(warn = -1)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "e", "f"), c(0, 0, 0, 0, 0), TRUE), c("a", "b", "c"))
  options(warn = 0)
})
test_that("All duplicate columns, different terms", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = a, "c" = a, "d" = a, "e" = a, "f" = a)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "e", "f"), c(0, 1, 2, 3, 4), TRUE), c("a", "b", "c", "e", "f"))
})
test_that("Repeated duplicate columns", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = a, "e" = a, "f" = a)
  options(warn = -1)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c", "d", "f"), c(0, 0, 0, 0, 0), TRUE), c("a", "b", "c"))
  options(warn = 0)
})
test_that("All but one duplicate column with varying", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = a, "c" = a)
  options(warn = -1)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c"), c(0, 0, 0), TRUE), c("a"))
  options(warn = 0)
})
test_that("All but one duplicate column with constant", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = c, "b" = c, "c" = c)
  options(warn = -1)
  expect_equal(Check_Dupe_Columns(df, c("a", "b", "c"), c(0, 0, 0), TRUE), c())
  options(warn = 0)
})
test_that("Duplicate with column not in df error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = c, "b" = c, "c" = c)
  options(warn = -1)
  expect_error(Check_Dupe_Columns(df, c("a", "b", "c", "e"), c(0, 0, 0, 0), TRUE))
  expect_error(Check_Dupe_Columns(df, c("a", "e", "c", "c"), c(0, 0, 0, 0), TRUE))
  options(warn = 0)
})

## ------------------------------------- ##
## LRT
## ------------------------------------- ##

test_that("Improve Ratio test", {
  a <- list("LogLik" = -400)
  b <- list("LogLik" = -350)
  expect_equal(Likelihood_Ratio_Test(b, a)$value, 100)
})
test_that("Worse Ratio test", {
  a <- list("LogLik" = -300)
  b <- list("LogLik" = -350)
  expect_equal(Likelihood_Ratio_Test(b, a)$value, -100)
})
test_that("Same Ratio test", {
  a <- list("LogLik" = -300)
  b <- list("LogLik" = -300)
  expect_equal(Likelihood_Ratio_Test(a, b)$value, 0)
})
test_that("No Data Ratio test", {
  a <- list("baditem" = -300)
  b <- list("LogLik" = -300)
  expect_error(Likelihood_Ratio_Test(a, b))
})

######################################
# FACTORING
######################################

test_that("Factorize factor", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)
  col_list <- c("c")
  expect_equal(factorize(df, col_list, TRUE)$cols, c("c_1"))
})
test_that("Factorize discrete", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  df <- data.table("a" = a, "b" = b, "c" = c)
  col_list <- c("a")
  expect_equal(factorize(df, col_list, TRUE)$cols, c("a_0", "a_1", "a_2", "a_3", "a_4", "a_5", "a_6"))
})
test_that("Factorize missing", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(0, 0, 0, 0, 0, 0, 0)
  df <- data.table("a" = a, "b" = b, "c" = c)
  col_list <- c("d")
  expect_error(factorize(df, col_list, TRUE))
})
test_that("Factorize survival lung, test", {
  if (system.file(package = "survival") != "") {
    data(cancer, package = "survival")
    cancer %>% setDT()
    df <- copy(cancer)
    col_list <- c("inst")
    expect_no_error(factorize(df, col_list, TRUE))
  }
})


######################################
# Time Dependent Cov gens
######################################

test_that("Gen_time_dep time error", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a_bad"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin")


  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 2))
})
test_that("Gen_time_dep event error", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c_bad"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin")


  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 2))
})
test_that("Gen_time_dep function error", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c_bad"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    stop()
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin")


  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 2))
})
test_that("Gen_time_dep functional form error", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("badbad")


  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
})

test_that("Gen_time_dep no error lin cox", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin")


  expect_no_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new", sep = ""), func_form, 2))
})
test_that("Gen_time_dep, error length names, tform, func_form", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin", "lin", "lin", "lin")


  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
  func_form <- c("lin")
  expect_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f, grt_f, grt_f, grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
})
test_that("Gen_time_dep no error step cox", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("step?0g?7l?12a?18b?")


  expect_no_error(gen_time_dep(df, time1, time2, event, TRUE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
})

test_that("Gen_time_dep no error lin not cox", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("lin")


  expect_no_error(gen_time_dep(df, time1, time2, event, FALSE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
})
test_that("Gen_time_dep no error step not cox", {
  a <- c(20, 20, 5, 10, 15)
  b <- c(1, 2, 1, 1, 2)
  c <- c(0, 0, 1, 1, 1)
  df <- data.table("a" = a, "b" = b, "c" = c)

  time1 <- "%trunc%"
  time2 <- "a"
  event <- "c"
  control <- list("lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9, "deriv_epsilon" = 1e-9, "step_max" = 1.0, "change_all" = TRUE, "thres_step_max" = 100.0, "verbose" = 0, "ties" = "breslow")
  grt_f <- function(df, time_col) {
    return((df[, "b"] * df[, get(time_col)])[[1]])
  }
  func_form <- c("step?0g?7l?10u?12a?18b?")


  expect_no_error(gen_time_dep(df, time1, time2, event, FALSE, 0.01, c("grt"), c(), c(grt_f), paste(tempfile(), "test", "_new.csv", sep = ""), func_form, 2))
})

test_that("linked quad negative slope error", {
  tforms <- list("first" = "quad")
  paras <- list("first" = c(-0.1, 10))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked quad string slope error", {
  tforms <- list("first" = "quad")
  paras <- list("first" = c("a", 10))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked quad string threshold error", {
  tforms <- list("first" = "quad")
  paras <- list("first" = c(0.1, "a"))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked quad no error", {
  tforms <- list("first" = "quad")
  paras <- list("first" = c(0.1, 10))
  expect_no_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked exp negative slope error", {
  tforms <- list("first" = "exp")
  paras <- list("first" = c(-0.1, 10, 5))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked exp string slope error", {
  tforms <- list("first" = "exp")
  paras <- list("first" = c("a", 10, 5))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked exp string threshold error", {
  tforms <- list("first" = "exp")
  paras <- list("first" = c(0.1, "a", 5))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked exp string exp slope error", {
  tforms <- list("first" = "exp")
  paras <- list("first" = c(0.1, 10, "a"))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked exp no error", {
  tforms <- list("first" = "exp")
  paras <- list("first" = c(0.1, 10, 5))
  expect_no_error(Linked_Dose_Formula(tforms, paras, TRUE))
})
test_that("linked formula combinations", {
  tforms <- list("first" = "quad")
  paras <- list("first" = c(0.1, 10))
  expect_error(Linked_Dose_Formula(tforms, paras, verbose = "p"))
  paras <- list("first" = c(0.1, "10"))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
  #
  tforms <- list("first" = "exp")
  paras <- list("first" = c(0.1, "10", 5))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
  paras <- list("first" = c(0.1, 10, "5"))
  expect_error(Linked_Dose_Formula(tforms, paras, TRUE))
})

test_that("linked exp parameter low goal error", {
  y <- 10
  a0 <- 1
  a_goal <- 5
  expect_error(Linked_Lin_Exp_Para(y, a0, a_goal, TRUE))
})
test_that("linked exp parameter negative slope error", {
  y <- 10
  a0 <- -0.1
  a_goal <- 5
  expect_error(Linked_Lin_Exp_Para(y, a0, a_goal, TRUE))
})
test_that("linked exp parameter no error", {
  y <- 10
  a0 <- 0.1
  a_goal <- 5
  expect_no_error(Linked_Lin_Exp_Para(y, a0, a_goal, TRUE))
})

test_that("Missing Value missing column error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_error(Replace_Missing(df, c("a", "e"), 0.0, TRUE))
})
test_that("Missing Value NA replacement error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_error(Replace_Missing(df, c("a", "b", "c", "d"), NA, TRUE))
})
test_that("Missing Value no error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_no_error(Replace_Missing(df, c("a", "b", "c", "d"), 0.0, TRUE))
})
test_that("Missing Value checked replaced 0", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(NA, 0, 0, 1, 0, 0, 1)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)

  df0 <- Replace_Missing(df, c("a", "b"), 0.0, TRUE)
  expect_equal(c(sum(df0$a), sum(df0$b)), c(sum(df$a), 2))
})
test_that("Missing Value checked replaced 1", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(NA, 0, 0, 1, 0, 0, 1)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)

  df0 <- Replace_Missing(df, c("a", "b"), 1.0, TRUE)
  expect_equal(c(sum(df0$a), sum(df0$b)), c(sum(df$a), 3))
})
test_that("Missing Value verbose error", {
  a <- c(0, 1, 2, 3, 4, 5, 6)
  b <- c(1, 2, 3, 4, 5, 6, 7)
  c <- c(1, 1, 1, 1, 1, 1, 1)
  d <- c(3, 4, 5, 6, 7, 8, 9)
  df <- data.table("a" = a, "b" = b, "c" = c, "d" = d)
  expect_no_error(Replace_Missing(df, c("a", "b", "c", "d"), 0.0, TRUE))
  expect_error(Replace_Missing(df, c("a", "b", "c", "d"), 0.0, -1))
})

# test_that( "Check Date Shift", {
#     m0 <- c(1,1,2,2)
#     m1 <- c(2,2,3,3)
#     d0 <- c(1,2,3,4)
#     d1 <- c(6,7,8,9)
#     y0 <- c(1990,1991,1997,1998)
#     y1 <- c(2001,2003,2005,2006)
#     df <- data.table( "m0"=m0, "m1"=m1, "d0"=d0, "d1"=d1, "y0"=y0, "y1"=y1)
#     expect_no_error(Date_Shift(df,c( "m0", "d0", "y0" ),c( "m1", "d1", "y1" ), "date_since" ))
# })
test_that("Check Date Shift, exact value", {
  m0 <- c(1, 1, 2, 2)
  m1 <- c(2, 2, 3, 3)
  d0 <- c(1, 2, 3, 4)
  d1 <- c(6, 7, 8, 9)
  y0 <- c(1990, 1991, 1997, 1998)
  y1 <- c(2001, 2003, 2005, 2006)
  df <- data.table("m0" = m0, "m1" = m1, "d0" = d0, "d1" = d1, "y0" = y0, "y1" = y1)
  e <- Date_Shift(df, c("m0", "d0", "y0"), c("m1", "d1", "y1"), "date_since")
  expect_equal(as.numeric(e$date_since), c(4054, 4419, 2955, 2955))
})

# test_that( "Check Date Since", {
#     m0 <- c(1,1,2,2)
#     m1 <- c(2,2,3,3)
#     d0 <- c(1,2,3,4)
#     d1 <- c(6,7,8,9)
#     y0 <- c(1990,1991,1997,1998)
#     y1 <- c(2001,2003,2005,2006)
#     df <- data.table( "m0"=m0, "m1"=m1, "d0"=d0, "d1"=d1, "y0"=y0, "y1"=y1)
#     tref <- strptime( "3-22-1997", format = "%m-%d-%Y",tz = 'UTC' )
#     expect_no_error(Time_Since(df,c( "m1", "d1", "y1" ),tref, "date_since" ))
# })
test_that("Check Date Since", {
  m0 <- c(1, 1, 2, 2)
  m1 <- c(2, 2, 3, 3)
  d0 <- c(1, 2, 3, 4)
  d1 <- c(6, 7, 8, 9)
  y0 <- c(1990, 1991, 1997, 1998)
  y1 <- c(2001, 2003, 2005, 2006)
  df <- data.table("m0" = m0, "m1" = m1, "d0" = d0, "d1" = d1, "y0" = y0, "y1" = y1)
  tref <- "3-22-1997"
  expect_error(Time_Since(df, c("m1", "d1", "y1"), tref, "date_since"))
})
test_that("Check Date Since, exact value", {
  m0 <- c(1, 1, 2, 2)
  m1 <- c(2, 2, 3, 3)
  d0 <- c(1, 2, 3, 4)
  d1 <- c(6, 7, 8, 9)
  y0 <- c(1990, 1991, 1997, 1998)
  y1 <- c(2001, 2003, 2005, 2006)
  df <- data.table("m0" = m0, "m1" = m1, "d0" = d0, "d1" = d1, "y0" = y0, "y1" = y1)
  tref <- strptime("3-22-1997", format = "%m-%d-%Y", tz = "UTC")
  e <- Time_Since(df, c("m1", "d1", "y1"), tref, "date_since")
  expect_equal(as.numeric(e$date_since), c(1417, 2148, 2908, 3274))
})
