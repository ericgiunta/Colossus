#' Performs basic Logistic regression using the omnibus function
#'
#' \code{RunLogisticRegression_Omnibus} uses user provided data, time/event columns,
#'  vectors specifying the model, and options to control the convergence and starting positions.
#'  Has additional options to starting with several initial guesses
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Logistic Wrapper Functions
#' @importFrom rlang .data
RunLogisticRegression_Omnibus <- function(df, trial0 = "CONST", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  initial_size <- nrow(df)
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  #
  if ("CONST" %in% c(trial0, names)) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  #
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  df <- df[get(trial0) > 0, ]
  if (min(df[, event0, with = FALSE]) < 0) {
    stop("Error: negative events in atleast one row")
  }
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  df0 <- data.table::data.table(a = c(0, 0))
  val <- list(cols = c("a"))
  val_cols <- c("a")
  data.table::setkeyv(df, c(event0, trial0))
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(event0, trial0)
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
  run_size <- nrow(df)
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
    if (model_control[["bisect"]]) {
      para_number <- model_control$para_number
      keep_constant[para_number] <- 1
      if (min(keep_constant) == 1) {
        model_control["single"] <- TRUE
      }
      e <- logist_Omnibus_CurveSearch_transition(
        as.matrix(df[, ce, with = FALSE]),
        term_n, tform, a_ns,
        dfc, x_all, 0,
        modelform, control,
        keep_constant, term_tot, model_control,
        cons_mat, cons_vec
      )
    } else {
      e <- logist_Omnibus_Bounds_transition(
        as.matrix(df[, ce, with = FALSE]),
        term_n, tform, a_ns,
        dfc, x_all, 0,
        modelform, control,
        keep_constant, term_tot, model_control,
        cons_mat, cons_vec
      )
    }
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
    e <- logist_Omnibus_transition(
      as.matrix(df[, ce, with = FALSE]),
      term_n, tform, matrix(a_ns,
        nrow = length(control$maxiters) - 1,
        byrow = TRUE
      ), dfc, x_all, 0,
      modelform, control, keep_constant,
      term_tot,
      model_control, cons_mat, cons_vec
    )
    if (is.nan(e$LogLik)) {
      stop(e$Status)
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Logistic"
  e$modelcontrol <- model_control
  e$control <- control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  e$UsedRecords <- run_size
  e$RejectedRecords <- initial_size - run_size
  e
}

#' Performs Logistic regression using the omnibus function with multiple event realizations
#'
#' \code{RunLogisticRegression_Omnibus_Multioutcome} uses user provided data, time/event
#'       columns, vectors specifying the model, and options to control the convergence and
#'       starting position.
#'       Returns optimized parameters, log-likelihood, and standard deviation for each realization.
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return returns a list of the final results
#' @family Logistic Wrapper Functions
#' @importFrom rlang .data
RunLogisticRegression_Omnibus_Multioutcome <- function(df, trial0 = "CONST", event0 = "event0", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", realization_columns = c("event0", "event1"), control = list(), model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  initial_size <- nrow(df)
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  #
  if ("CONST" %in% c(trial0, names)) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  #
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  event_names <- unique(as.vector(realization_columns))
  if (!all(event_names %in% names(df))) {
    stop("Error: Atleast one realization column provided was not in the data.table")
  }
  df <- df[get(trial0) > 0, ]
  if (!is.numeric(df[[trial0]])) {
    stop("Error: Trial column was not numeric: ", trial0)
  }
  for (i in seq(1, length(realization_columns))) {
    if (!is.numeric(df[[realization_columns[i]]])) {
      stop("Error: Event column was not numeric: ", realization_columns[i])
    }
    if (min(df[[realization_columns[i]]]) < 0) {
      stop("Error: Negative events in atleast one row of column ", realization_columns[i])
    }
    if (sum(df[[realization_columns[i]]]) == 0) {
      stop("Error: No events in column ", realization_columns[i])
    }
  }
  df0 <- data.table::data.table(a = c(0, 0))
  val <- list(cols = c("a"))
  val_cols <- c("a")
  data.table::setkeyv(df, trial0)
  #
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  event_cols <- match(realization_columns, event_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  event_all <- as.matrix(df[, event_names, with = FALSE])
  ce <- c(event0, trial0)
  run_size <- nrow(df)
  e <- logist_multioutcome_Omnibus_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_n,
    event_cols, dfc, x_all, event_all,
    0, modelform, control,
    keep_constant, term_tot,
    model_control, cons_mat, cons_vec
  )
  if ("Status" %in% names(e)) {
    if (all(e$Status != "PASSED")) {
      stop("Error: Every realization failed.")
    } else if (any(e$Status != "PASSED")) {
      if (control$verbose >= 2) {
        warning("Warning: Atleast one realization failed.")
      }
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Logistic_Multidose"
  e$modelcontrol <- model_control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  e$UsedRecords <- run_size
  e$RejectedRecords <- initial_size - run_size
  e
}

#' Performs Logistic regression using the omnibus function with multiple column realizations
#'
#' \code{RunLogisticRegression_Omnibus_Multidose} uses user provided data, time/event columns,
#'       vectors specifying the model, and options to control the convergence
#'       and starting positions. Used for 2DMC column uncertainty methods.
#'       Returns optimized parameters, log-likelihood, and standard error for each realization.
#'       Has additional options for using stratification
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return returns a list of the final results for each realization
#' @family Logistic Wrapper Functions
#' @importFrom rlang .data
RunLogisticRegression_Omnibus_Multidose <- function(df, trial0 = "CONST", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  initial_size <- nrow(df)
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  #
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  #
  to_remove <- c("CONST", "%trunc%")
  to_keep <- c(trial0, event0, names, realization_index, as.vector(realization_columns))
  to_keep <- unique(to_keep)
  to_keep <- to_keep[!to_keep %in% to_remove]
  to_keep <- to_keep[to_keep %in% names(df)]
  df <- df[, to_keep, with = FALSE]
  #
  ce <- c(event0, trial0)
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
  df <- df[get(trial0) > 0, ]
  data.table::setkeyv(df, c(event0, trial0))
  #
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  # replace_missing equivalent for the realization columns
  if (length(realization_index) == length(realization_columns[, 1])) {
    # pass
  } else {
    # the number of columns per realization does not match the number of indexes provided
    stop(
      "Error: ", length(realization_index),
      " column indexes provided, but ",
      length(realization_columns[, 1]),
      " rows of realizations columns provided"
    )
  }
  if (!all(realization_index %in% all_names)) {
    stop("Error: Atleast one realization column provided was not used in the model")
  }
  dose_names <- unique(as.vector(realization_columns))
  if (!all(dose_names %in% names(df))) {
    stop("Error: Atleast one realization column provided was not in the data.table")
  }
  dfc <- match(names, all_names)
  dose_cols <- matrix(match(realization_columns, dose_names), nrow = nrow(realization_columns))
  dose_index <- match(realization_index, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  dose_all <- as.matrix(df[, dose_names, with = FALSE])
  run_size <- nrow(df)
  e <- logist_multidose_Omnibus_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_n,
    as.matrix(dose_cols, with = FALSE), dose_index, dfc, x_all, dose_all,
    0, modelform, control,
    keep_constant, term_tot,
    model_control, cons_mat, cons_vec
  )
  if ("Status" %in% names(e)) {
    if (all(e$Status != "PASSED")) {
      stop("Error: Every realization failed.")
    } else if (any(e$Status != "PASSED")) {
      if (control$verbose >= 2) {
        warning("Warning: Atleast one realization failed.") # nocov
      }
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$modelcontrol <- model_control
  e$control <- control
  e$Survival_Type <- "Logistic_Multidose"
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  e$UsedRecords <- run_size
  e$RejectedRecords <- initial_size - run_size
  e
}

#' Calculates logistic residuals
#'
#' \code{RunLogisticRegression_Residual} uses user provided data, time/event columns,
#'  vectors specifying the model, and options. Calculates residuals or sum of residuals
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Poisson Wrapper Functions
#' @importFrom rlang .data
RunLogisticRegression_Residual <- function(df, trial0 = "trial", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), model_control = list()) {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  cons_mat <- as.matrix(c(0))
  cons_vec <- c(0)
  control <- Def_Control(control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  if (min(df[, trial0, with = FALSE]) < 0) {
    stop("Error: negative trials in atleast one row")
  }
  if (min(df[, event0, with = FALSE]) < 0) {
    stop("Error: negative events in atleast one row")
  }
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  if (sum(df[, trial0, with = FALSE]) == 0) {
    stop("Error: no duration")
  }

  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  #
  #  df$og_order <- seq_len(nrow(df))
  #  data.table::setkeyv(df, val_cols)
  #  print(df)
  #
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(event0, trial0)
  e <- logist_Residual_transition(
    as.matrix(df[, ce, with = FALSE]),
    term_n, tform, a_n[[1]],
    dfc, x_all, 0, modelform,
    control, keep_constant,
    term_tot,
    model_control
  )
  #  extra_names <- names(e)
  #  df_risk <- data.table("index" = df$og_order)
  #  for (ex_name in extra_names) {
  #    df_risk[[ex_name]] <- e[[ex_name]]
  #  }
  #  data.table::setkeyv(df_risk, "index")
  #  for (ex_name in extra_names) {
  #    e[[ex_name]] <- df_risk[[ex_name]]
  #  }
  e
}
