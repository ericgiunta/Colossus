test_that("threshold nonfail", {
  data(cancer, package = "survival")
  df <- veteran

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
  # Convert the cell column into multiple factor columns
  fcols <- c("cell")
  val <- factorize(df, fcols) # Colossus function
  df <- val$df

  t0 <- "%trunc%"
  t1 <- "time"
  event <- "status"

  names <- c(
    "karno50", "trt"
  )
  tform_1 <- c(
    "loglin", "loglin"
  )

  term_n <- c(0, 0)
  a_n <- c(0.1, 0.1)

  control <- list(verbose = 0, abs_max = 0.1)
  devs <- c(2357.78433, 1097.69940, 1097.69940, 3384.48849, 1445.98094, 1445.98094, 120.50896, 113.28929, 107.78937, 123.02794, 123.02794, 123.02794, 3307.64803, 862.15791, 850.39045, 662.89835, 619.59109, 645.78672, 662.89835, 619.59109, 645.78672, 4708.10011, 1139.98323, 1209.77460, 961.67111, 918.33615, 970.65878, 961.67111, 918.33615, 970.65878, 75.15047, 60.14696, 61.22819, 98.03485, 52.86674, 56.29411, 73.57696, 49.86964, 56.35650, 207.36346, 62.25028, 64.19787, 207.36346, 62.25028, 64.19787, 207.36346, 62.25028, 64.19787, 887.89770, 662.22318, 662.22318, 1167.77892, 961.67111, 961.67111, 64.42914, 57.35914, 54.43160, 66.40498, 66.40498, 66.40498, 847.13843, 619.59108, 619.59108, 1125.09711, 918.33064, 918.33064, 59.78693, 52.75154, 49.85893, 62.08826, 62.08826, 62.08826)
  free_strat <- c(113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1, 113, 113, 113, 0, 0, 0, 0, 0, 0, 96, 96, 96, 0, 0, 0, 0, 0, 0, 4, 4, 4, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1, 113, 0, 0, 96, 0, 0, 4, 1, 0, 1, 1, 1)

  i_index <- 1

  for (extra_bool in c("single", "gradient", "null", "pass")) {
    for (time_bool in c(T, F)) {
      for (strat_bool in c(T, F)) {
        for (thres in c(0, 40, 100)) {
          if (extra_bool == "gradient") {
            for (method in c("momentum", "adadelta", "adam")) {
              model_control <- list("time_risk" = time_bool, "strata" = strat_bool, "conditional_threshold" = thres)
              model_control[extra_bool] <- TRUE
              model_control[[method]] <- TRUE
              e <- RunCaseControlRegression_Omnibus(
                df, t0, t1, event,
                names = names, tform = tform_1,
                strat_col = "cell", model_control = model_control,
                control = control, term_n = term_n, a_n = a_n
              )
              expect_equal(devs[i_index], e$Deviance, tolerance = 1e-3)
              expect_equal(free_strat[i_index], e$FreeSets, tolerance = 1e-3)
              i_index <- i_index + 1
            }
          } else {
            model_control <- list("time_risk" = time_bool, "strata" = strat_bool, "conditional_threshold" = thres)
            model_control[extra_bool] <- TRUE
            e <- RunCaseControlRegression_Omnibus(
              df, t0, t1, event,
              names = names, tform = tform_1,
              strat_col = "cell", model_control = model_control,
              control = control, term_n = term_n, a_n = a_n
            )
            expect_equal(devs[i_index], e$Deviance, tolerance = 1e-3)
            expect_equal(free_strat[i_index], e$FreeSets, tolerance = 1e-3)
            i_index <- i_index + 1
          }
        }
      }
    }
  }
})
