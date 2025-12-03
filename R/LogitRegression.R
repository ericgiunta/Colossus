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
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' df <- data.table::data.table(
#'   "Trials" = c(30, 45, 57, 47, 36, 60, 55),
#'   "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
#'   "a" = c(0, 1, 1, 0, 1, 0, 1),
#'   "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
#'   "c" = c(10, 11, 10, 11, 12, 9, 11),
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' # For the interval case
#' trial <- "Trials"
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
#' e <- RunLogisticRegression_Omnibus(
#'   df, trial, event, names, term_n,
#'   tform, keep_constant,
#'   a_n, modelform,
#'   control
#' )
#' @importFrom rlang .data
RunLogisticRegression_Omnibus <- function(df, trial0 = "CONST", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
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
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (sum(df[, event0, with = FALSE]) == 0) {
    stop("Error: no events")
  }
  df0 <- data.table::data.table("a" = c(0, 0))
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
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Logistic"
  if (is.nan(e$LogLik)) {
    stop(e$Status)
  }
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}
