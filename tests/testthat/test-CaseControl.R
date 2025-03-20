test_that("Basic nonfail", {
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

  control <- list(verbose = 0)

  for (time_bool in c(F)) {
    for (strat_bool in c(F)) {
      model_control <- list("time_risk" = time_bool, "strata" = strat_bool)
      expect_no_error(RunCaseControlRegression_Omnibus(
        df, t0, t1, event,
        names = names, tform = tform_1,
        strat_col = "cell", model_control = model_control,
        control = control, term_n = term_n, a_n = a_n
      ))
    }
  }
})
