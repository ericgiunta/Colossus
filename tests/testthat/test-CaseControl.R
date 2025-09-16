test_that("threshold nonfail", {
  data(cancer, package = "survival")
  veteran %>% setDT()
  df <- copy(veteran)

  # Make the same adjustments as Epicure example 6.5
  karno <- df$karno
  karno[93] <- 20
  df$karno <- karno
  df$trt <- df$trt - 1
  df$trt <- as.integer(df$trt == 0)
  cell_string <- df$celltype
  cell <- case_when(
    cell_string == "squamous" ~ 1,
    cell_string == "smallcell" ~ 2,
    cell_string == "adeno" ~ 3,
    cell_string == "large" ~ 0
  )
  df$cell <- cell

  df$karno50 <- df$karno - 50
  a_n <- c(0.1, 0.1)


  control <- list(verbose = 0, step_max = 0.1, ncores = 2)
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
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
  }
})

test_that("threshold nonfail, single", {
  data(cancer, package = "survival")
  veteran %>% setDT()
  df <- copy(veteran)

  # Make the same adjustments as Epicure example 6.5
  karno <- df$karno
  karno[93] <- 20
  df$karno <- karno
  df$trt <- df$trt - 1
  df$trt <- as.integer(df$trt == 0)
  cell_string <- df$celltype
  cell <- case_when(
    cell_string == "squamous" ~ 1,
    cell_string == "smallcell" ~ 2,
    cell_string == "adeno" ~ 3,
    cell_string == "large" ~ 0
  )
  df$cell <- cell

  df$karno50 <- df$karno - 50
  a_n <- c(0.1, 0.1)


  control <- list(verbose = 0, step_max = 0.1, ncores = 2)
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
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
  }
})

test_that("threshold nonfail, null", {
  data(cancer, package = "survival")
  veteran %>% setDT()
  df <- copy(veteran)

  # Make the same adjustments as Epicure example 6.5
  karno <- df$karno
  karno[93] <- 20
  df$karno <- karno
  df$trt <- df$trt - 1
  df$trt <- as.integer(df$trt == 0)
  cell_string <- df$celltype
  cell <- case_when(
    cell_string == "squamous" ~ 1,
    cell_string == "smallcell" ~ 2,
    cell_string == "adeno" ~ 3,
    cell_string == "large" ~ 0
  )
  df$cell <- cell

  df$karno50 <- df$karno - 50
  a_n <- c(0.1, 0.1)


  control <- list(verbose = 0, step_max = 0.1, ncores = 2)
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
        expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
        expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
        i_index <- i_index + 1
      }
    }
  }
})

test_that("threshold nonfail, gradient", {
  data(cancer, package = "survival")
  veteran %>% setDT()
  df <- copy(veteran)

  # Make the same adjustments as Epicure example 6.5
  karno <- df$karno
  karno[93] <- 20
  df$karno <- karno
  df$trt <- df$trt - 1
  df$trt <- as.integer(df$trt == 0)
  cell_string <- df$celltype
  cell <- case_when(
    cell_string == "squamous" ~ 1,
    cell_string == "smallcell" ~ 2,
    cell_string == "adeno" ~ 3,
    cell_string == "large" ~ 0
  )
  df$cell <- cell

  df$karno50 <- df$karno - 50
  a_n <- c(0.1, 0.1)


  control <- list(verbose = 0, step_max = 0.1, ncores = 2)
  #
  i_index <- 1
  #
  devs <- c(2367.0306, 862.1579, 850.3904, 1107.9470, 619.5911, 620.3378, 1107.9470, 619.5911, 620.3378)
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
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
  devs <- c(3117.1567, 1139.9832, 1366.5472, 1452.8518, 918.3362, 928.9564, 1452.8518, 918.3362, 928.9564)
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
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
  devs <- c(68.79205, 60.14696, 69.71407, 57.34329, 52.86674, 53.64409, 99.15252, 49.86964, 50.35396)
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
          expect_equal(e$Deviance, devs[i_index], tolerance = 1e-3)
          expect_equal(e$FreeSets, free_strat[i_index], tolerance = 1e-3)
          i_index <- i_index + 1
        }
      }
    }
  }
  devs <- c(128.79809, 62.25028, 66.40201, 128.79809, 62.25028, 66.40201, 128.79809, 62.25028, 66.40201)
  free_strat <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)
  i_index <- 1
  for (time_bool in c(F)) {
    for (strat_bool in c(F)) {
      model <- CaseCon(status) ~ loglinear(karno50, trt)
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
    }
  }
})
