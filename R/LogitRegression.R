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
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  e$UsedRecords <- run_size
  e$RejectedRecords <- initial_size - run_size
  e
}
