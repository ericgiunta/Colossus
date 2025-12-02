test_that("Checking basic function", {
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
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2)
    #
    expect_no_error(res <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control))
    expect_no_error(res_max <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control))
    expect_no_error(res_mean <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "mean", control = control))
    #
    expect_false(isTRUE(all.equal(res$First_Der, res_max$First_Der)))
    expect_false(isTRUE(all.equal(res$First_Der, res_mean$First_Der)))
    #
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), control = control, cons_mat = cons_mat0, cons_vec = c(0.0)))
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control, cons_mat = cons_mat0, cons_vec = c(0.0)))
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "mean", control = control, cons_mat = cons_mat0, cons_vec = c(0.0)))
    #
    expect_no_error(res_max <- PoisRun(Pois(time, status) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control))
    expect_no_error(res_max <- CaseControlRun(CaseCon(status) ~ loglinear(karno50, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control))
    #
    expect_no_error(res_max <- PoisRun(Pois(time, status) ~ loglinear(karno, trt), df, a_n = list(c(-0.1, 0.1), c(0.1, 0.1)), norm = "max", control = control))
  }
})

test_that("Checking values converted back", {
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
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2, maxiter = -1, halfmax = -1)
    for (i in c(0.1, -0.1, -0.035, 0.035)) {
      a_n <- c(i, -1 * i)
      res <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, control = control)
      res_max <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, norm = "max", control = control)
      res_mean <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, norm = "mean", control = control)
      #
      expect_equal(res$First_Der, res_max$First_Der, tolerance = 1e-2)
      expect_equal(res$First_Der, res_mean$First_Der, tolerance = 1e-2)
      #
      expect_equal(res$Second_Der, res_max$Second_Der, tolerance = 1e-2)
      expect_equal(res$Second_Der, res_mean$Second_Der, tolerance = 1e-2)
      #
      expect_equal(res$beta_0, res_max$beta_0, tolerance = 1e-2)
      expect_equal(res$beta_0, res_mean$beta_0, tolerance = 1e-2)
      #
      expect_equal(res$Standard_Deviation, res_max$Standard_Deviation, tolerance = 1e-2)
      expect_equal(res$Standard_Deviation, res_mean$Standard_Deviation, tolerance = 1e-2)
      #
      expect_equal(res$Covariance, res_max$Covariance, tolerance = 1e-2)
      expect_equal(res$Covariance, res_mean$Covariance, tolerance = 1e-2)
      #
      res <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, control = control, cons_mat = cons_mat0, cons_vec = c(0.0))
      res_max <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, norm = "max", control = control, cons_mat = cons_mat0, cons_vec = c(0.0))
      res_mean <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = a_n, norm = "mean", control = control, cons_mat = cons_mat0, cons_vec = c(0.0))
      #
      expect_equal(res$constraint_matrix, res_max$constraint_matrix, tolerance = 1e-2)
      expect_equal(res$constraint_matrix, res_mean$constraint_matrix, tolerance = 1e-2)
      #
    }
    res <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt) + plinear(karno), df, a_n = a_n, control = control)
    res_max <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt) + plinear(karno), df, a_n = a_n, norm = "max", control = control)
    expect_equal(res$beta_0, res_max$beta_0, tolerance = 1e-2)
  }
  #
  model <- Cox_Strata(time, status, cell) ~ loglinear(trt) + step - dose(karno, 1) + ME()
  a_n <- c(0.1, 0.2, 50)
  res <- CoxRun(model, df, a_n = a_n, control = control)
  res_max <- CoxRun(model, df, a_n = a_n, norm = "max", control = control)
  res_mean <- CoxRun(model, df, a_n = a_n, norm = "mean", control = control)
  #
  expect_equal(res$First_Der, res_max$First_Der, tolerance = 1e-2)
  expect_equal(res$First_Der, res_mean$First_Der, tolerance = 1e-2)
})

test_that("Checking combination with gradient/single/null", {
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
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2, maxiter = -1, halfmax = -1)
    #
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ null(), df, control = control, norm = "max"))
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control, gradient_control = list(momentum = TRUE)))
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, a_n = c(-0.1, 0.1), norm = "max", control = control, single = TRUE))
  }
})

test_that("Checking errors and warnings", {
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
    cons_mat0 <- matrix(c(1, 1), nrow = 1)
    control <- list(ncores = 2, maxiter = -1, halfmax = -1)

    df$trt_mean <- df$trt - mean(df$trt)
    #
    options(warn = -1)
    expect_no_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt_mean), df, a_n = c(-0.1, 0.1), norm = "max", control = control))
    expect_error(CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt_mean), df, a_n = c(-0.1, 0.1), norm = "bad", control = control))
    expect_error(PoisRun(Pois(time, status) ~ loglinear(karno, trt_mean), df, a_n = c(-0.1, 0.1), norm = "bad", control = control))
    options(warn = 0)
  }
})

test_that("Checking likelihood bound", {
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
    control <- list(ncores = 2)
    #
    res <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, control = control)
    res_max <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, norm = "max", control = control)
    res_mean <- CoxRun(Cox_Strata(time, status, cell) ~ loglinear(karno, trt), df, norm = "mean", control = control)
    #
    expect_no_error(bound <- LikelihoodBound(res, df))
    expect_no_error(bound_max <- LikelihoodBound(res_max, df))
    expect_no_error(bound_mean <- LikelihoodBound(res_mean, df))
    #
    expect_equal(bound$Parameter_Limits, bound_max$Parameter_Limits, tolerance = 1e-3)
    expect_equal(bound$Parameter_Limits, bound_mean$Parameter_Limits, tolerance = 1e-3)
    #
    res <- PoisRun(Pois_Strata(time, status, cell) ~ loglinear(karno, trt), df, control = control)
    res_max <- PoisRun(Pois_Strata(time, status, cell) ~ loglinear(karno, trt), df, norm = "max", control = control)
    res_mean <- PoisRun(Pois_Strata(time, status, cell) ~ loglinear(karno, trt), df, norm = "mean", control = control)
    #
    expect_no_error(bound <- LikelihoodBound(res, df))
    expect_no_error(bound_max <- LikelihoodBound(res_max, df))
    expect_no_error(bound_mean <- LikelihoodBound(res_mean, df))
    #
    expect_equal(bound$Parameter_Limits, bound_max$Parameter_Limits, tolerance = 1e-3)
    expect_equal(bound$Parameter_Limits, bound_mean$Parameter_Limits, tolerance = 1e-3)
  }
})
