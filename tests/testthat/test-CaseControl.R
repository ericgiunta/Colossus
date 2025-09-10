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
  devs <- c(3307.64803, 862.15791, 850.39045, 662.89835, 619.59109, 645.78672, 662.89835, 619.59109, 645.78672, 4708.10011, 1139.98323, 1209.77460, 961.67111, 918.33615, 970.65878, 961.67111, 918.33615, 970.65878, 75.15047, 60.14696, 61.22819, 98.03485, 52.86674, 56.29411, 73.57696, 49.86964, 56.35650, 207.36346, 62.25028, 64.19787, 207.36346, 62.25028, 64.19787, 207.36346, 62.25028, 64.19787)
  free_strat <- c(113, 113, 113, 0, 0, 0, 0, 0, 0, 96, 96, 96, 0, 0, 0, 0, 0, 0, 4, 4, 4, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)

  extra_bool <- "gradient"
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
