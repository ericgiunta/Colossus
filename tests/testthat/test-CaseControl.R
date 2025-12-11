test_that("check errors", {
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
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(847.13843, 619.59108, 619.59108, 1125.09711, 918.33064, 918.33064, 59.78693, 52.75154, 49.85893, 62.08826, 62.08826, 62.08826)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)
    extra_bool <- "pass"
    for (time_bool in c(T, F)) {
      for (strat_bool in c(T, F)) {
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
          expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
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
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(2357.7843, 1097.6994, 1097.6994, 3384.4885, 1445.9809, 1445.9809, 120.5090, 113.2893, 107.7894, 123.0279, 123.0279, 123.0279)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)

    extra_bool <- "single"
    for (time_bool in c(T, F)) {
      for (strat_bool in c(T, F)) {
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
          expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
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
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(887.89770, 662.22318, 662.22318, 1167.77892, 961.67111, 961.67111, 64.42914, 57.35914, 54.43160, 66.40498, 66.40498, 66.40498)
    free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)

    extra_bool <- "null"
    for (time_bool in c(T, F)) {
      for (strat_bool in c(T, F)) {
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
          expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
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
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    devs <- c(3290.233, 865.7832, 990.8613, 732.6527, 619.5911, 659.3394, 732.6527, 619.5911, 659.3394)
    free_strat <- c(113, 113, 113, 0, 0, 0, 0, 0, 0)
    i_index <- 1
    extra_bool <- "gradient"
    for (time_bool in c(T)) {
      for (strat_bool in c(T)) {
        model <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
        for (thres in c(0, 40, 100)) {
          for (method in c("momentum", "adadelta", "adam")) {
            gradient_control <- list()
            gradient_control[[method]] <- TRUE
            e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
            #
            expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
            expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
            i_index <- i_index + 1
          }
        }
      }
    }
    devs <- c(8811.506, 1146.021, 1587.319, 961.6711, 918.3543, 963.5374, 961.6711, 918.3543, 963.5374)
    free_strat <- c(96, 96, 96, 0, 0, 0, 0, 0, 0)
    i_index <- 1
    for (time_bool in c(T)) {
      for (strat_bool in c(F)) {
        model <- CaseCon_Time(time, status) ~ loglinear(karno50, trt)
        for (thres in c(0, 40, 100)) {
          for (method in c("momentum", "adadelta", "adam")) {
            gradient_control <- list()
            gradient_control[[method]] <- TRUE
            e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
            #
            expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
            expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
            i_index <- i_index + 1
          }
        }
      }
    }
    devs <- c(67.49968, 60.20573, 67.69674, 57.45191, 52.95613, 57.91124, 54.89959, 49.87332, 57.06744)
    free_strat <- c(4, 4, 4, 1, 1, 1, 0, 0, 0)
    i_index <- 1
    for (time_bool in c(F)) {
      for (strat_bool in c(T)) {
        model <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
        for (thres in c(0, 40, 100)) {
          for (method in c("momentum", "adadelta", "adam")) {
            gradient_control <- list()
            gradient_control[[method]] <- TRUE
            e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
            #
            expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
            expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
            i_index <- i_index + 1
          }
        }
      }
    }
    devs <- c(75.325, 62.48172, 63.38723, 113.11, 75.325, 62.48172, 63.38723, 113.11, 75.325, 62.48172, 63.38723, 113.11)
    free_strat <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    i_index <- 1
    for (time_bool in c(F)) {
      for (strat_bool in c(F)) {
        model <- CaseCon(status) ~ loglinear(karno50, trt)
        for (thres in c(0, 40, 100)) {
          for (method in c("momentum", "adadelta", "adam", "none")) {
            gradient_control <- list()
            gradient_control[[method]] <- TRUE
            e <- CaseControlRun(model, df, gradient_control = gradient_control, control = control, conditional_threshold = thres, a_n = a_n)
            #
            expect_equal(e$Deviation, devs[i_index], tolerance = 1e-3)
            expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
            i_index <- i_index + 1
          }
        }
      }
    }
  }
})

test_that("information matrix calculations", {
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
    a_n <- c(0.1, 0.1)


    control <- list(verbose = 0, step_max = 0.1, ncores = 1)
    #
    i_index <- 1
    #
    model0 <- CaseCon_Strata(status, cell) ~ loglinear(karno50, trt)
    model1 <- CaseCon_Strata(status, cell) ~ loglinear(karno50) + plinear(trt)
    model2 <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50, trt)
    model3 <- CaseCon_Strata_Time(time, status, cell) ~ loglinear(karno50) + plinear(trt)


    model0_sdo <- c(0.01609599, 0.54817862, 0.02048499, 0.57345304, 0.02315814, 0.73703345)
    model1_sdo <- c(0.01374865, 0.42821895, 0.02048498, 0.39593058, 0.02315809, 0.51108894)
    model2_sdo <- c(0.004699080, 0.137902029, 0.005553695, 0.201870857, 0.005553695, 0.201870857)
    model3_sdo <- c(0.004420340, 0.176224467, 0.005553695, 0.160828177, 0.005553695, 0.160828177)
    model0_sde <- c(0.01609599, 0.54817862, 0.02048499, 0.57345304, 0.02315814, 0.73703345)
    model1_sde <- c(0.01405449, 0.47879518, 0.02048499, 0.39593268, 0.02315808, 0.51107695)
    model2_sde <- c(0.004699080, 0.137902029, 0.005553695, 0.201870857, 0.005553695, 0.201870857)
    model3_sde <- c(0.004250846, 0.148569581, 0.005553695, 0.160827320, 0.005553695, 0.160827320)


    i_index <- 1
    for (thres in c(0, 40, 100)) {
      eo <- CaseControlRun(model0, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = T)
      ee <- CaseControlRun(model0, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = F)
      expect_equal(eo$Standard_Deviation[1], model0_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Deviation[2], model0_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[1], model0_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[2], model0_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model1, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = T)
      ee <- CaseControlRun(model1, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = F)
      expect_equal(eo$Standard_Deviation[1], model1_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Deviation[2], model1_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[1], model1_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[2], model1_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model2, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = T)
      ee <- CaseControlRun(model2, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = F)
      expect_equal(eo$Standard_Deviation[1], model2_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Deviation[2], model2_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[1], model2_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[2], model2_sde[i_index + 1], tolerance = 1e-4)
      #
      eo <- CaseControlRun(model3, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = T)
      ee <- CaseControlRun(model3, df, control = control, conditional_threshold = thres, a_n = a_n, observed_info = F)
      expect_equal(eo$Standard_Deviation[1], model3_sdo[i_index], tolerance = 1e-4)
      expect_equal(eo$Standard_Deviation[2], model3_sdo[i_index + 1], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[1], model3_sde[i_index], tolerance = 1e-4)
      expect_equal(ee$Standard_Deviation[2], model3_sde[i_index + 1], tolerance = 1e-4)
      i_index <- i_index + 2
    }
  }
})
