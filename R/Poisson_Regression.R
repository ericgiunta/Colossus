#' Performs basic Poisson regression using the omnibus function
#'
#' \code{RunPoissonRegression_Omnibus} uses user provided data, time/event columns,
#'  vectors specifying the model, and options to control the convergence and starting positions.
#'  Has additional options to starting with several initial guesses
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Poisson Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
#'   "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
#'   "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
#'   "a" = c(0, 1, 1, 0, 1, 0, 1),
#'   "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
#'   "c" = c(10, 11, 10, 11, 12, 9, 11),
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' # For the interval case
#' pyr <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c("a", "b", "c", "d")
#' a_n <- c(1.1, -0.1, 0.2, 0.5) # used to test at a specific point
#' term_n <- c(0, 1, 1, 2)
#' tform <- c("loglin", "lin", "lin", "plin")
#' modelform <- "M"
#' keep_constant <- c(0, 0, 0, 0)
#' control <- list(
#'   "ncores" = 2, "lr" = 0.75, "maxiter" = 5,
#'   "halfmax" = 5, "epsilon" = 1e-3,
#'   "deriv_epsilon" = 1e-3, "step_max" = 1.0,
#'   "thres_step_max" = 100.0, "verbose" = FALSE, "ties" = "breslow",
#'   "double_step" = 1
#' )
#' strat_col <- "e"
#' e <- RunPoissonRegression_Omnibus(
#'   df, pyr, event, names, term_n,
#'   tform, keep_constant,
#'   a_n, modelform,
#'   control, strat_col
#' )
#' @importFrom rlang .data
RunPoissonRegression_Omnibus <- function(df, pyr0 = "pyr", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  df <- df[get(pyr0) > 0, ]
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$strata == TRUE) {
    val <- factorize(df, strat_col)
    df0 <- val$df
    df <- val$df
    val_cols <- c()
    for (col in val$cols) {
      dftemp <- df[get(col) == 1, ]
      temp <- sum(dftemp[, get(event0)])
      if (temp == 0) {
        if (control$verbose >= 2) {
          warning(paste("Warning: no events for strata group:", col,
            sep = " "
          ))
        }
        df <- df[get(col) != 1, ]
        df0 <- df0[get(col) != 1, ]
      } else {
        val_cols <- c(val_cols, col)
      }
      data.table::setkeyv(df0, c(pyr0, event0))
    }
  } else {
    df0 <- data.table::data.table("a" = c(0, 0))
    val <- list(cols = c("a"))
    val_cols <- c("a")
  }
  data.table::setkeyv(df, c(pyr0, event0))
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(pyr0, event0)
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
  if (model_control$log_bound) {
    if ("maxiters" %in% names(control)) {
      # good
    } else {
      control$maxiters <- c(control$maxiter)
    }
    if ("guesses" %in% names(control)) {
      # good
    } else {
      control$guesses <- 10
    }
    e <- pois_Omnibus_Bounds_transition(
      as.matrix(df[, ce, with = FALSE]),
      term_n, tform, a_ns, dfc, x_all, 0, modelform, control,
      keep_constant, term_tot, as.matrix(df0[, val_cols, with = FALSE]),
      model_control, cons_mat, cons_vec
    )
    if ("Status" %in% names(e)) {
      if (e$Status != "PASSED") {
        stop(e$Status)
      }
    }
  } else {
    res <- Check_Iters(control, a_n)
    control <- res$control
    a_n <- res$a_n
    #
    e <- pois_Omnibus_transition(
      as.matrix(df[, ce, with = FALSE]),
      term_n, tform, matrix(a_ns,
        nrow = length(control$maxiters) - 1,
        byrow = TRUE
      ), dfc, x_all, 0,
      modelform, control, keep_constant,
      term_tot, as.matrix(df0[, val_cols,
        with = FALSE
      ]),
      model_control, cons_mat, cons_vec
    )
    if (is.nan(e$LogLik)) {
      stop(e$Status)
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Poisson"
  e$modelcontrol <- model_control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}

#' Predicts how many events are due to baseline vs excess
#'
#' \code{RunPoissonEventAssignment} uses user provided data, person-year/event columns, vectors specifying the model, and options to calculate background and excess events
#'
#' @noRd
#' @inheritParams R_template
#' @family Poisson Wrapper Functions
#' @return returns a list of the final results
#'
RunPoissonEventAssignment <- function(df, pyr0 = "pyr", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", model_control = list()) {
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  control$maxiters <- c(1, control$maxiter)
  control$guesses <- 1
  control <- Def_Control(control)
  df <- df[get(pyr0) > 0, ]
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$strata == TRUE) {
    val <- factorize(df, strat_col)
    df0 <- val$df
    df <- val$df
    val_cols <- c()
    for (col in val$cols) {
      dftemp <- df[get(col) == 1, ]
      temp <- sum(dftemp[, get(event0)])
      if (temp == 0) {
        if (control$verbose >= 2) {
          warning(paste("Warning: no events for strata group:", col,
            sep = " "
          ))
        }
        #        df <- df[get(col) != 1, ]
        #        df0 <- df0[get(col) != 1, ]
      } else {
        val_cols <- c(val_cols, col)
      }
      #      data.table::setkeyv(df0, c(pyr0, event0))
    }
  } else {
    df0 <- data.table::data.table("a" = c(0, 0))
    val <- list(cols = c("a"))
    val_cols <- c("a")
  }
  #  data.table::setkeyv(df, c(pyr0, event0))
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(pyr0, event0)
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
  e <- Assigned_Event_Poisson_transition(
    as.matrix(df[, ce, with = FALSE]),
    as.matrix(df0[, val_cols,
      with = FALSE
    ]), term_n, tform,
    a_n, dfc, x_all, 0,
    modelform, control, keep_constant,
    term_tot, model_control
  )
  return(e)
}

#' Calculates poisson residuals
#'
#' \code{RunPoissonRegression_Residual} uses user provided data, time/event columns,
#'  vectors specifying the model, and options. Calculates residuals or sum of residuals
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Poisson Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
#'   "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
#'   "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
#'   "a" = c(0, 1, 1, 0, 1, 0, 1),
#'   "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
#'   "c" = c(10, 11, 10, 11, 12, 9, 11),
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' # For the interval case
#' pyr <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c("a", "b", "c", "d")
#' a_n <- c(1.1, -0.1, 0.2, 0.5) # used to test at a specific point
#' term_n <- c(0, 1, 1, 2)
#' tform <- c("loglin", "lin", "lin", "plin")
#' modelform <- "M"
#' keep_constant <- c(0, 0, 0, 0)
#' control <- list(
#'   "ncores" = 2, "lr" = 0.75, "maxiter" = 5,
#'   "halfmax" = 5, "epsilon" = 1e-3,
#'   "deriv_epsilon" = 1e-3, "step_max" = 1.0,
#'   "thres_step_max" = 100.0, "verbose" = FALSE, "ties" = "breslow",
#'   "double_step" = 1
#' )
#' strat_col <- "e"
#' e <- RunPoissonRegression_Residual(
#'   df, pyr, event, names, term_n,
#'   tform, keep_constant,
#'   a_n, modelform,
#'   control, strat_col
#' )
#' @importFrom rlang .data
RunPoissonRegression_Residual <- function(df, pyr0 = "pyr", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", model_control = list()) {
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  cons_mat <- as.matrix(c(0))
  cons_vec <- c(0)
  control <- Def_Control(control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  df <- df[get(pyr0) > 0, ]
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$strata == TRUE) {
    val <- factorize(df, strat_col)
    df0 <- val$df
    df <- val$df
    val_cols <- c()
    for (col in val$cols) {
      dftemp <- df[get(col) == 1, ]
      temp <- sum(dftemp[, get(event0)])
      if (temp == 0) {
        if (control$verbose >= 2) {
          warning(paste("Warning: no events for strata group:",
            col,
            sep = " "
          ))
        }
        #        df <- df[get(col) != 1, ]
        #        df0 <- df0[get(col) != 1, ]
      } else {
        val_cols <- c(val_cols, col)
      }
      #      data.table::setkeyv(df0, c(pyr0, event0))
    }
  } else {
    df0 <- data.table::data.table("a" = c(0, 0))
    val <- list(cols = c("a"))
    val_cols <- c("a")
  }
  #  data.table::setkeyv(df, c(pyr0, event0))
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(pyr0, event0)
  e <- pois_Residual_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_n[[1]],
    dfc, x_all, 0, modelform,
    control, keep_constant,
    term_tot, as.matrix(df0[, val_cols,
      with = FALSE
    ]),
    model_control
  )
  return(e)
}

#' Calculates the likelihood curve for a poisson model directly
#'
#' \code{PoissonCurveSolver} solves the confidence interval for a poisson model, starting at the optimum point and
#' iteratively optimizing each point to using the bisection method
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Poisson Wrapper Functions
#' @importFrom rlang .data
PoissonCurveSolver <- function(df, pyr0 = "pyr", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  df <- df[get(pyr0) > 0, ]
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$strata == TRUE) {
    val <- factorize(df, strat_col)
    df0 <- val$df
    df <- val$df
    val_cols <- c()
    for (col in val$cols) {
      dftemp <- df[get(col) == 1, ]
      temp <- sum(dftemp[, get(event0)])
      if (temp == 0) {
        if (control$verbose >= 2) {
          warning(paste("Warning: no events for strata group:", col,
            sep = " "
          ))
        }
        df <- df[get(col) != 1, ]
        df0 <- df0[get(col) != 1, ]
      } else {
        val_cols <- c(val_cols, col)
      }
      data.table::setkeyv(df0, c(pyr0, event0))
    }
  } else {
    df0 <- data.table::data.table("a" = c(0, 0))
    val <- list(cols = c("a"))
    val_cols <- c("a")
  }
  data.table::setkeyv(df, c(pyr0, event0))
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(pyr0, event0)
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
  res <- Check_Iters(control, a_n)
  control <- res$control
  a_n <- res$a_n
  if ("alpha" %in% names(model_control)) {
    model_control["qchi"] <- qchisq(1 - model_control[["alpha"]], df = 1) / 2
  } else {
    model_control["alpha"] <- 0.05
    model_control["qchi"] <- qchisq(1 - model_control[["alpha"]], df = 1) / 2
  }
  para_number <- model_control$para_number
  keep_constant[para_number] <- 1
  if (min(keep_constant) == 1) {
    model_control["single"] <- TRUE
  }
  e <- pois_Omnibus_CurveSearch_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_ns, dfc, x_all, 0, modelform, control,
    keep_constant, term_tot, as.matrix(df0[, val_cols, with = FALSE]),
    model_control, cons_mat, cons_vec
  )
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Poisson"
  e$modelcontrol <- model_control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}

#' Performs Poisson regression using the omnibus function with multiple column realizations
#'
#' \code{RunPoisRegression_Omnibus_Multidose} uses user provided data, time/event columns,
#'       vectors specifying the model, and options to control the convergence
#'       and starting positions. Used for 2DMC column uncertainty methods.
#'       Returns optimized parameters, log-likelihood, and standard deviation for each realization.
#'       Has additional options for using stratification
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return returns a list of the final results for each realization
#' @family Poisson Wrapper Functions
#' @importFrom rlang .data
RunPoisRegression_Omnibus_Multidose <- function(df, pyr0 = "pyr", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), strat_col = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  #
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  #
  to_remove <- c("CONST", "%trunc%")
  to_keep <- c(pyr0, event0, names, realization_index, as.vector(realization_columns))
  if (model_control$strata == TRUE) {
    to_keep <- c(to_keep, strat_col)
  }
  to_keep <- unique(to_keep)
  to_keep <- to_keep[!to_keep %in% to_remove]
  to_keep <- to_keep[to_keep %in% names(df)]
  df <- df[, to_keep, with = FALSE]
  #
  ce <- c(pyr0, event0)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if ("MCML" %in% names(model_control)) {
    # fine
  } else {
    model_control$MCML <- FALSE
  }
  if (model_control$strata == TRUE) {
    val <- factorize(df, strat_col)
    df0 <- val$df
    df <- val$df
    val_cols <- c()
    for (col in val$cols) {
      dftemp <- df[get(col) == 1, ]
      temp <- sum(dftemp[, get(event0)])
      if (temp == 0) {
        if (control$verbose >= 2) {
          warning(paste("Warning: no events for strata group:", col,
            sep = " "
          ))
        }
        df <- df[get(col) != 1, ]
        df0 <- df0[get(col) != 1, ]
      } else {
        val_cols <- c(val_cols, col)
      }
      data.table::setkeyv(df0, c(pyr0, event0))
    }
  } else {
    df0 <- data.table::data.table("a" = c(0, 0))
    val <- list(cols = c("a"))
    val_cols <- c("a")
  }
  data.table::setkeyv(df, c(pyr0, event0))
  #
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  # replace_missing equivalent for the realization columns
  if (length(realization_index) == length(realization_columns[, 1])) {
    # pass
  } else {
    # the number of columns per realization does not match the number of indexes provided
    stop(paste("Error:", length(realization_index),
      " column indexes provided, but ",
      length(realization_columns[, 1]),
      " rows of realizations columns provided",
      sep = " "
    ))
  }
  if (all(realization_index %in% all_names)) {
    # pass
  } else {
    stop(paste("Error: Atleast one realization column provided was not used in the model", sep = " "))
  }
  #  all_names <- unique(c(all_names, as.vector(realization_columns)))
  dose_names <- unique(as.vector(realization_columns))
  if (all(dose_names %in% names(df))) {
    # pass
  } else {
    stop(paste("Error: Atleast one realization column provided was not in the data.table", sep = " "))
  }
  dfc <- match(names, all_names)
  dose_cols <- matrix(match(realization_columns, dose_names), nrow = nrow(realization_columns))
  dose_index <- match(realization_index, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  dose_all <- as.matrix(df[, dose_names, with = FALSE])
  e <- pois_multidose_Omnibus_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_n,
    as.matrix(dose_cols, with = FALSE), dose_index, dfc, x_all, dose_all,
    0, modelform, control,
    keep_constant, term_tot, as.matrix(df0[, val_cols, with = FALSE]),
    model_control, cons_mat, cons_vec
  )
  if ("Status" %in% names(e)) {
    if (e$Status != "PASSED") {
      stop(e$Status)
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  if (model_control$MCML) {
    e$Survival_Type <- "Pois_Multidose"
  } else {
    e$Survival_Type <- "Pois_Multidose"
  }
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  # df <- copy(df)
  return(e)
}
