#' Performs Matched Case-Control Conditional Logistic Regression
#'
#' \code{RunCaseControlRegression_Omnibus} uses user provided data, time/event columns,
#'   vectors specifying the model, and options to control the convergence
#'   and starting positions. Has additional options for starting with several
#'   initial guesses, using stratification and/or matching by time at risk,
#'   and calculation without derivatives
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Case Control Wrapper Functions
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
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c("a", "b", "c", "d")
#' a_n <- list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4))
#' # used to test at a specific point
#' term_n <- c(0, 1, 1, 2)
#' tform <- c("loglin", "lin", "lin", "plin")
#' modelform <- "M"
#' keep_constant <- c(0, 0, 0, 0)
#' control <- list(
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(5, 5, 5),
#'   "halfmax" = 5, "epsilon" = 1e-3, "deriv_epsilon" = 1e-3,
#'   "step_max" = 1.0, "thres_step_max" = 100.0,
#'   "verbose" = FALSE,
#'   "ties" = "breslow", "double_step" = 1
#' )
#' e <- RunCaseControlRegression_Omnibus(df, time1, time2, event,
#'   names, term_n, tform, keep_constant,
#'   a_n, modelform, control,
#'   model_control = list(
#'     "stata" = FALSE,
#'     "time_risk" = FALSE
#'   )
#' )
#' @importFrom rlang .data
RunCaseControlRegression_Omnibus <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", cens_weight = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
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
  if (model_control$time_risk == TRUE) {
    ce <- c(time1, time2, event0)
    t_check <- Check_Trunc(df, ce)
    df <- t_check$df
    ce <- t_check$ce
    ## Cox regression only uses intervals which contain an event time
    time1 <- ce[1]
    time2 <- ce[2]
    dfend <- df[get(event0) == 1, ]
    tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
    if (length(tu) == 0) {
      stop("Error: no events")
    }
    # remove rows that end before first event
    df <- df[get(time2) >= tu[1], ]
    # remove rows that start after the last event
    df <- df[get(time1) <= tu[length(tu)], ]
  }
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  #  if (model_control$linear_err == TRUE) {
  #    if (all(sort(unique(tform)) != c("loglin", "plin"))) {
  #      stop("Error: Linear ERR model used, but term formula wasn't only loglin and plin")
  #    }
  #    if (sum(tform == "plin") > 1) {
  #      stop("Error: Linear ERR model used, but more than one plin element was used")
  #    }
  #    if (length(unique(term_n)) > 1) {
  #      if (control$verbose >= 2) {
  #        warning("Warning: Linear ERR model used, but more than one term number used. Term numbers all set to 0")
  #      }
  #      term_n <- rep(0, length(term_n))
  #    }
  #    if (modelform != "M") {
  #      if (control$verbose >= 2) {
  #        warning("Warning: Linear ERR model used, but multiplicative model not used. Modelform corrected")
  #      }
  #      modelform <- "M"
  #    }
  #  }
  if (model_control$time_risk == TRUE) {
    if (model_control$strata == TRUE) {
      if (!is.null(levels(df[[strat_col]]))) {
        # The column is a factor, so we can convert to numbers
        factor_lvl <- levels(df[[strat_col]])
        df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
      } else if (is(typeof(df[[strat_col]]), "character")) {
        df[[strat_col]] <- factor(df[[strat_col]])
        factor_lvl <- levels(df[[strat_col]])
        df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
      }
      dfend <- df[get(event0) == 1, ]
      uniq_end <- unlist(unique(dfend[, strat_col, with = FALSE]),
        use.names = FALSE
      )
      df <- df[get(strat_col) %in% uniq_end, ]
      uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
        use.names = FALSE
      ))
      if (control$verbose >= 3) {
        message(paste("Note:", length(uniq), " strata used", sep = " ")) # nocov
      }
      data.table::setkeyv(df, c(strat_col, event0, time2, time1))
      ce <- c(time1, time2, strat_col, event0)
    } else {
      data.table::setkeyv(df, c(event0, time2, time1))
      uniq <- c(0)
      ce <- c(time1, time2, event0)
    }
  } else {
    if (model_control$strata == TRUE) {
      if (!is.null(levels(df[[strat_col]]))) {
        # The column is a factor, so we can convert to numbers
        factor_lvl <- levels(df[[strat_col]])
        df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
      } else if (is(typeof(df[[strat_col]]), "character")) {
        df[[strat_col]] <- factor(df[[strat_col]])
        factor_lvl <- levels(df[[strat_col]])
        df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
      }
      dfend <- df[get(event0) == 1, ]
      uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
        use.names = FALSE
      ))
      for (i in seq_along(uniq)) {
        df0 <- dfend[get(strat_col) == uniq[i], ]
        if (nrow(df0) == 0) {
          if (control$verbose >= 2) { # nocov
            warning(paste("Warning: no events for strata group:",
              uniq[i],
              sep = " "
            )) # nocov
          }
          df <- df[get(strat_col) != uniq[i], ] # nocov
        }
      }
      uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
        use.names = FALSE
      ))
      if (control$verbose >= 3) {
        message(paste("Note:", length(uniq), " strata used", sep = " ")) # nocov
      }
      data.table::setkeyv(df, c(strat_col, event0))
      ce <- c(strat_col, event0)
    } else {
      data.table::setkeyv(df, c(event0))
      ce <- c(event0)
      uniq <- c(0)
    }
  }
  if (model_control$time_risk == TRUE) {
    dfend <- df[get(event0) == 1, ]
    tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
    if (control$verbose >= 3) {
      message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
    }
  } else {
    tu <- c(0)
  }
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  # make sure any constant 0 columns are constant
  for (i in seq_along(keep_constant)) {
    if ((keep_constant[i] == 0) && (names[i] %in% names(df))) {
      if (names[i] != "CONST") {
        if (min(df[[names[i]]]) == max(df[[names[i]]])) {
          keep_constant[i] <- 1 # nocov
          if (control$verbose >= 2) { # nocov
            warning(paste("Warning: element ", i,
              " with column name ", names[i],
              " was set constant",
              sep = ""
            )) # nocov
          }
        }
      }
    }
  }
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free") # nocov
  }
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
  res <- Check_Iters(control, a_n)
  control <- res$control
  a_n <- res$a_n
  if (model_control$null) {
    a_ns <- matrix(a_ns)
  } else {
    a_ns <- matrix(a_ns, nrow = length(control$maxiters) - 1, byrow = TRUE)
  }
  #
  e <- caco_Omnibus_transition(
    term_n, tform, a_ns, dfc, x_all, 0,
    modelform, control, as.matrix(df[, ce, with = FALSE]), tu,
    keep_constant, term_tot, uniq, model_control,
    cons_mat, cons_vec
  )
  if ("Status" %in% names(e)) {
    if (is.nan(e$LogLik)) {
      stop(e$Status)
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "CaseControl"
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}
