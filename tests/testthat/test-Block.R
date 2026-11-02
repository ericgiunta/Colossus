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

    control <- list(verbose = 0, step_max = 0.1, ncores = 1, maxiters = c(20, 20))
    cx_model <- Cox(time, status) ~          loglinear(factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    ps_model <- Pois(time, status) ~         loglinear(CONST, factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    lo_model <- Logit(CONST, status) ~       loglinear(CONST, factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    cc_model <- CaseCon_Time(time, status) ~ loglinear(factor(cell), 0) + linear(karno50, 1) + loglinear(trt, 1) + ME()
    # full term subterm individual bad
    cx_ll <- c(-519.0727, -483.9966, -484.2420, -488.2954)
    ps_ll <- c(-402.7944, -7136.1817, -3961.1570, -2924.7688)
    lo_ll <- c(-31.13919, -78.71926, -85.59210, -85.07171)
    cc_ll <- c(-494.0147, -490.0290, -459.6962, -460.2688)
    i_val <- 1
    for (des in c("full", "term", "subterm", "individual")) {
      expect_no_error(e0 <- CoxRun(cx_model, df, control = control, descent = des))
      expect_no_error(e1 <- PoisRun(ps_model, df, control = control, descent = des))
      expect_no_error(e2 <- LogisticRun(lo_model, df, control = control, descent = des))
      expect_no_error(e3 <- CaseControlRun(cc_model, df, control = control, descent = des))
      expect_equal(cx_ll[i_val], e0$LogLik, tolerance = 1e-3)
      expect_equal(ps_ll[i_val], e1$LogLik, tolerance = 1e-3)
      expect_equal(lo_ll[i_val], e2$LogLik, tolerance = 1e-3)
      expect_equal(cc_ll[i_val], e3$LogLik, tolerance = 1e-3)
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

    control <- list(verbose = 0, step_max = 0.1, ncores = 1, maxiters = c(20, 20))
    cx_model <- Cox(time, status) ~ loglinear(karno50, 0) + loglinear - dose(trt, 1) + linear - dose(trt, 1) + quad - dose(trt, 1) + step - dose(trt, 1) + linear - quadratic - dose(trt, 1) + linear - exponential - dose(trt, 1) + ME()

    # full term subterm individual bad
    cx_ll <- c(-490.6792, -484.4244, -484.9368, -485.8349)
    i_val <- 1
    for (des in c("full", "term", "subterm", "individual")) {
      expect_no_error(e0 <- CoxRun(cx_model, df, control = control, descent = des))
      expect_equal(cx_ll[i_val], e0$LogLik, tolerance = 1e-3)
      i_val <- i_val + 1
    }
  }
})
