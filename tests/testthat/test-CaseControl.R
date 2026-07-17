test_that("check errors", {
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
    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
    expect_no_error(CaseControlRun(model, df, control = control, norm = "mean"))
    expect_error(CaseControlRun(model, df, control = control, norm = "bad"))
    expect_error(CaseControlRun(model, df, control = control, bad = "wrong"))
  }
})

test_that("threshold nonfail", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(847.13843, 619.59108, 619.59108, 1125.09712, 918.33064, 918.33064, 61.45791, 52.75154, 49.85893, 63.53347, 63.53347, 63.53347)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)
    extra_bool <- "pass"
    for (time_bool in c(TRUE, FALSE)) {
      for (strat_bool in c(TRUE, FALSE)) {
        if (time_bool) {
          if (strat_bool) {
            model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon_Time(time, status) ~ loglinear(karno50, trt)
          }
        } else {
          if (strat_bool) {
            model <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon(status) ~ loglinear(karno50, trt)
          }
        }
        for (thres in c(0, 40, 100)) {
          e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, a_n = a_n)
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
})

test_that("threshold nonfail, single", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(2357.7843, 1097.6994, 1097.6994, 3384.4885, 1445.9809, 1445.9809, 120.5090, 113.2893, 107.7894, 123.0279, 123.0279, 123.0279)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)

    extra_bool <- "single"
    for (time_bool in c(TRUE, FALSE)) {
      for (strat_bool in c(TRUE, FALSE)) {
        if (time_bool) {
          if (strat_bool) {
            model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon_Time(time, status) ~ loglinear(karno50, trt)
          }
        } else {
          if (strat_bool) {
            model <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
          } else {
            model <- CaseCon(status) ~ loglinear(karno50, trt)
          }
        }
        for (thres in c(0, 40, 100)) {
          e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, single = TRUE, a_n = a_n)
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
})

test_that("threshold nonfail, null", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(887.89770, 662.22318, 662.22318, 1167.77892, 961.67111, 961.67111, 64.42914, 57.35914, 54.43160, 66.40498, 66.40498, 66.40498)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)

    extra_bool <- "null"
    for (time_bool in c(TRUE, FALSE)) {
      for (strat_bool in c(TRUE, FALSE)) {
        if (time_bool) {
          if (strat_bool) {
            model <- CaseCon_Strata_Time(time, status, cell) ~ null()
          } else {
            model <- CaseCon_Time(time, status) ~ null()
          }
        } else {
          if (strat_bool) {
            model <- CaseCon_Strata(status, cell) ~ null()
          } else {
            model <- CaseCon(status) ~ null()
          }
        }
        for (thres in c(0, 40, 100)) {
          e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, single = TRUE, a_n = a_n)
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
})

test_that("threshold nonfail, gradient", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(2365.2638, 861.7430, 895.6137, 662.8983, 619.5911, 645.7867, 662.8983, 619.5911, 645.7867)
    free_strat <- c(113, 113, 113, 0, 0, 0, 0, 0, 0)
    i_index <- 1
    extra_bool <- "gradient"
    time_bool <- TRUE
    strat_bool <- TRUE
    model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
    for (thres in c(0, 40, 100)) {
      for (method in c("momentum", "adadelta", "adam")) {
        gradient_control <- list()
        gradient_control[[method]] <- TRUE
        e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
        #
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
    devs <- c(3112.6391, 1139.7693, 1603.9994, 961.6711, 918.3317, 928.9564, 961.6711, 918.3317, 928.9564)
    free_strat <- c(96, 96, 96, 0, 0, 0, 0, 0, 0)
    i_index <- 1
    time_bool <- TRUE
    strat_bool <- FALSE
    model <- CaseCon_Time(time, status) ~ loglinear(karno50, trt)
    for (thres in c(0, 40, 100)) {
      for (method in c("momentum", "adadelta", "adam")) {
        gradient_control <- list()
        gradient_control[[method]] <- TRUE
        e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
        #
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
    devs <- c(68.68225, 60.01931, 63.52628, 57.36262, 52.81594, 53.76418, 99.15252, 49.86315, 50.35396)
    free_strat <- c(4, 4, 4, 1, 1, 1, 0, 0, 0)
    i_index <- 1
    time_bool <- FALSE
    strat_bool <- TRUE
    model <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
    for (thres in c(0, 40, 100)) {
      for (method in c("momentum", "adadelta", "adam")) {
        gradient_control <- list()
        gradient_control[[method]] <- TRUE
        e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
        #
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
    devs <- c(69.95769, 62.18752, 63.28939, 113.10995, 69.95769, 62.18752, 63.28939, 113.10995, 69.95769, 62.18752, 63.28939, 113.10995)
    free_strat <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    i_index <- 1
    time_bool <- FALSE
    strat_bool <- FALSE
    model <- CaseCon(status) ~ loglinear(karno50, trt)
    for (thres in c(0, 40, 100)) {
      for (method in c("momentum", "adadelta", "adam", "none")) {
        gradient_control <- list()
        gradient_control[[method]] <- TRUE
        e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
        #
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
  }
})

test_that("information matrix calculations", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    model0 <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
    model1 <- CaseCon_Strata(status, cell) ~ loglinear(karno50) + plinear(trt)
    model2 <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
    model3 <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50) + plinear(trt)

    model0_sdo <- c(0.01403495, 0.49013322, 0.02048067, 0.57330581, 0.02315670, 0.73698750)
    model0_sde <- c(0.01403495, 0.49013322, 0.02048067, 0.57330581, 0.02315670, 0.73698750)
    model1_sdo <- c(0.01376595, 0.43243431, 0.02048033, 0.39577180, 0.02315695, 0.51160976)
    model1_sde <- c(0.01406310, 0.48225112, 0.02048070, 0.39590589, 0.02315678, 0.51141699)
    model2_sdo <- c(0.004699090, 0.137902133, 0.005553661, 0.201869822, 0.005553661, 0.201869822)
    model2_sde <- c(0.004699090, 0.137902133, 0.005553661, 0.201869822, 0.005553661, 0.201869822)
    model3_sdo <- c(0.00442034, 0.17622447, 0.00555366, 0.16083601, 0.00555366, 0.16083601)
    model3_sde <- c(0.004250846, 0.148569581, 0.005553660, 0.160832577, 0.005553660, 0.160832577)


    i_index <- 1
    for (thres in c(0, 40, 100)) {
      eo <- CaseControlRun(model0, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = TRUE)
      ee <- CaseControlRun(model0, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = FALSE)
      expect_equal(eo$Standard_Error[1], model0_sdo[i_index], tolerance = 1e-3)
      expect_equal(eo$Standard_Error[2], model0_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[1], model0_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[2], model0_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model1, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = TRUE)
      ee <- CaseControlRun(model1, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = FALSE)
      expect_equal(eo$Standard_Error[1], model1_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Error[2], model1_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[1], model1_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[2], model1_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model2, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = TRUE)
      ee <- CaseControlRun(model2, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = FALSE)
      expect_equal(eo$Standard_Error[1], model2_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Error[2], model2_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[1], model2_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[2], model2_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model3, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = TRUE)
      ee <- CaseControlRun(model3, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = FALSE)
      expect_equal(eo$Standard_Error[1], model3_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Error[2], model3_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[1], model3_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Error[2], model3_sde[i_index + 1], tolerance = 1e-4)
      i_index <- i_index + 2
    }
  }
})

test_that("constraint non-fail, gradient and hessian", {
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

    df$karno50 <- df$karno - 50
    a_n <- c(0.1, 0.1)
    cons_mat0 <- matrix(c(1, 1), nrow = 1)

    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    vals <- c(-0.03625298, 1.5, -0.1414994, 1.5, -0.03565239, -0.0737104, 0.0406554, -0.009756468, -0.03565239, -0.0737104, 0.0406554, -0.009756468)
    i_index <- 1
    extra_bool <- "gradient"
    time_bool <- TRUE
    strat_bool <- TRUE
    model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
    for (thres in c(0, 40, 100)) {
      e <- CaseControlRun(model, df, control = control, conditional_threshold = thres, a_n = a_n, cons_mat = cons_mat0, cons_vec = 0.0)
      expect_equal(sum(e$beta_0), 0.0, tolerance = 1e-4)
      expect_equal(e$beta_0[1], vals[i_index], tolerance = 1e-4)
      i_index <- i_index + 1
      for (method in c("momentum", "adadelta", "adam")) {
        gradient_control <- list()
        gradient_control[[method]] <- TRUE
        e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n, cons_mat = cons_mat0, cons_vec = 0.0)
        expect_equal(sum(e$beta_0), 0.0, tolerance = 1e-4)
        expect_equal(e$beta_0[1], vals[i_index], tolerance = 1e-4)
        i_index <- i_index + 1
      }
    }
  }
})
