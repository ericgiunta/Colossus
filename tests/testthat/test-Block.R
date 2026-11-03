test_that("Basic check and expected value", {
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
    df$status <- as.integer(df$status > 0)

    df$karno50 <- df$karno - 50

    control <- list(verbose = 0, step_max = 1.0, ncores = 1, maxiters = c(10, 10))
    cx_model <- Cox(time, status) ~          loglinear(factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    ps_model <- Pois(time, status) ~         loglinear(CONST, factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    lo_model <- Logit(CONST, status) ~       loglinear(CONST, factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    cc_model <- CaseCon_Time(time, status) ~ loglinear(factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    # full term subterm individual bad
    cx_ll <- c(-480.4063, -505.4994, -478.2765, -480.466)
    ps_ll <- c(-217.5047, -300.6482, -299.5807, -299.5324)
    lo_ll <- c(-30.4672, -33.83429, -29.58049, -32.52137)
    cc_ll <- c(-453.031, -480.3437, -453.0454, -454.0614)
    i_val <- 1
    for (des in c("full", "term", "subterm", "individual")) {
      expect_no_error(e0 <- CoxRun(cx_model, df, control = control, descent = des))
      expect_no_error(e1 <- PoisRun(ps_model, df, control = control, descent = des))
      expect_no_error(e2 <- LogisticRun(lo_model, df, control = control, descent = des))
      expect_no_error(e3 <- CaseControlRun(cc_model, df, control = control, descent = des))
      expect_equal(e0$LogLik, cx_ll[i_val], tolerance = 1e-3)
      expect_equal(e1$LogLik, ps_ll[i_val], tolerance = 1e-3)
      expect_equal(e2$LogLik, lo_ll[i_val], tolerance = 1e-3)
      expect_equal(e3$LogLik, cc_ll[i_val], tolerance = 1e-3)
      i_val <- i_val + 1
    }
  }
})

test_that("dose subterm checks", {
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
    df$status <- as.integer(df$status > 0)

    df$karno50 <- df$karno - 50

    cx_model <- Cox(time, status) ~ loglinear(karno50, 0) + loglinear - dose(trt, 1) + linear - dose(trt, 1) + quad - dose(trt, 1) + step - dose(trt, 1) + linear - quadratic - dose(trt, 1) + linear - exponential - dose(trt, 1) + ME()
    control <- list(verbose = 2, step_max = 1.0, ncores = 1, maxiters = c(1, 1))
    expect_warning(e0 <- CoxRun(cx_model, df, control = control, descent = "individual"))
    # full term subterm individual bad
    cx_ll <- c(-484.4242, -484.4243, -487.9913, -485.201)
    control <- list(verbose = 2, step_max = 1.0, ncores = 1, maxiters = c(-1, 14))
    i_val <- 1
    for (des in c("full", "term", "subterm", "individual")) {
      expect_no_error(e0 <- CoxRun(cx_model, df, control = control, descent = des))
      expect_equal(e0$LogLik, cx_ll[i_val], tolerance = 1e-3)
      i_val <- i_val + 1
    }
  }
})
