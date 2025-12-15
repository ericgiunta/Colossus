#' Fully runs a cox or fine-gray regression model, returning the model and results
#'
#' \code{CoxRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus cox or fine-gray regression function
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
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
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df,
#'   a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)),
#'   control = control
#' )
CoxRun <- function(model, df, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, observed_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), norm = "null", ...) {
  func_t_start <- Sys.time()
  if (is(model, "coxmodel")) {
    # using already prepped formula and data
    coxmodel <- copy(model)
    calls <- coxmodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    #
  } else if (is(model, "coxres")) {
    coxmodel <- model$model
    calls <- coxmodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    coxmodel$a_n <- model$beta_0
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    coxmodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or coxmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if (!missing(a_n)) {
    coxmodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    coxmodel$keep_constant <- keep_constant
  }
  #
  if ("CONST" %in% coxmodel$names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  ce <- c(coxmodel$start_age, coxmodel$end_age, coxmodel$event)
  val <- Check_Trunc(df, ce)
  if (any(val$ce != ce)) {
    df <- val$df
    ce <- val$ce
    coxmodel$start_age <- ce[1]
    coxmodel$end_age <- ce[1]
  }
  # Checks that the current coxmodel is valid
  validate_coxsurv(coxmodel, df)
  if (!coxmodel$null) {
    coxmodel <- validate_formula(coxmodel, df, control$verbose)
  }
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  a_n <- coxmodel$a_n
  modelform <- coxmodel$modelform
  strat_col <- coxmodel$strata
  cens_weight <- coxmodel$weight
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (coxmodel$null) {
    model_control["null"] <- TRUE
    #
    names <- c("CONST")
    term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0)
  } else {
    # check for basic and linear_err
    if (length(unique(term_n)) == 1) {
      modelform <- "M"
      if (all(unique(tform) == c("loglin"))) {
        model_control[["basic"]] <- TRUE
      } else if (identical(sort(unique(tform)), c("loglin", "plin")) && sum(tform == "plin") == 1) {
        model_control[["linear_err"]] <- TRUE
      }
    } else if (modelform == "GMIX") {
      model_control[["gmix_term"]] <- coxmodel$gmix_term
      model_control[["gmix_theta"]] <- coxmodel$gmix_theta
    }
  }
  if (all(coxmodel$strata != "NONE")) {
    model_control[["strata"]] <- TRUE
    #
    df$"_strata_col" <- format(df[, strat_col[1], with = FALSE]) # defining a strata column
    for (i in seq_len(length(strat_col) - 1)) {
      df$"_strata_col" <- paste(df$"_strata_col", format(df[, strat_col[i + 1], with = FALSE]), sep = "_") # interacting with any other strata columns
    }
    df$"_strata_col" <- factor(df$"_strata_col") # converting to a factor
  }
  if (coxmodel$weight != "NONE") {
    model_control[["cr"]] <- TRUE
  }
  if (!missing(cons_mat)) {
    model_control[["constraint"]] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  int_count <- 0.0
  if (!coxmodel$null) {
    norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
    a_n <- norm_res$a_n
    cons_mat <- norm_res$cons_mat
    norm_weight <- norm_res$norm_weight
    df <- norm_res$df
    if (any(norm_weight != 1.0)) {
      int_avg_weight <- 0.0
      for (i in seq_along(names)) {
        if (grepl("_int", tform[i])) {
          int_avg_weight <- int_avg_weight + norm_weight[i]
          int_count <- int_count + 1
        }
      }
      if (int_count > 0) {
        if (control$verbose >= 3) {
          message("Note: Threshold max step adjusted to match new weighting")
        }
        control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
      }
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, "_strata_col", cens_weight, model_control, cons_mat, cons_vec)
  if (int_count > 0) {
    control$thres_step_max <- control$thres_step_max * (int_avg_weight / int_count)
  }
  res$model <- coxmodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (!model_control$null) {
    if (model_control[["constraint"]]) {
      res$constraint_matrix <- cons_mat
      res$constraint_vector <- cons_vec
    }
    res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight, "tform" = tform), model_control)
  }
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  coxres <- new_coxres(res)
  coxres
}

#' Fully runs a poisson regression model, returning the model and results
#'
#' \code{PoisRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus poisson regression function
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Poisson Wrapper Functions
#' @examples
#' library(data.table)
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
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- Pois(Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- PoisRun(formula, df, a_n = c(1.1, -0.1, 0.2, 0.5), control = control)
PoisRun <- function(model, df, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, observed_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), norm = "null", ...) {
  func_t_start <- Sys.time()
  if (is(model, "poismodel")) {
    # using already prepped formula and data
    poismodel <- copy(model)
    calls <- poismodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
  } else if (is(model, "poisres")) {
    poismodel <- model$model
    calls <- poismodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    poismodel$a_n <- model$beta_0
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    poismodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or poismodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if (!missing(a_n)) {
    poismodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    poismodel$keep_constant <- keep_constant
  }
  #
  if ("CONST" %in% poismodel$names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  # Checks that the current coxmodel is valid
  validate_poissurv(poismodel, df)
  poismodel <- validate_formula(poismodel, df, control$verbose)
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  a_n <- poismodel$a_n
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (length(unique(term_n)) == 1) {
    modelform <- "M"
  } else if (modelform == "GMIX") {
    model_control[["gmix_term"]] <- poismodel$gmix_term
    model_control[["gmix_theta"]] <- poismodel$gmix_theta
  }
  if (all(poismodel$strata != "NONE")) {
    model_control["strata"] <- TRUE
  }
  if (ncol(cons_mat) > 1) {
    model_control["constraint"] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  int_count <- 0.0
  if (any(norm_weight != 1.0)) {
    int_avg_weight <- 0.0
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        int_avg_weight <- int_avg_weight + norm_weight[i]
        int_count <- int_count + 1
      }
    }
    if (int_count > 0) {
      if (control$verbose >= 3) {
        message("Note: Threshold max step adjusted to match new weighting")
      }
      control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  if (int_count > 0) {
    control$thres_step_max <- control$thres_step_max * (int_avg_weight / int_count)
  }
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight, "tform" = tform), model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  poisres <- new_poisres(res)
  poisres
}

#' Fully runs a logistic regression model, returning the model and results
#'
#' \code{LogisticRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus logistic regression function
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Logistic Wrapper Functions
#' @examples
#' library(data.table)
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
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- logit(Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- LogisticRun(formula, df, a_n = c(1.1, -0.1, 0.2, 0.5), control = control)
LogisticRun <- function(model, df, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), link = "odds", single = FALSE, observed_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), norm = "null", ...) {
  func_t_start <- Sys.time()
  if (is(model, "logitmodel")) {
    # using already prepped formula and data
    logitmodel <- copy(model)
    calls <- logitmodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    logitmodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or logitmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if (!missing(a_n)) {
    logitmodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    logitmodel$keep_constant <- keep_constant
  }
  #
  if ("CONST" %in% c(logitmodel$names, logitmodel$trials)) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  # Checks that the current coxmodel is valid
  validate_logitsurv(logitmodel, df)
  logitmodel <- validate_formula(logitmodel, df, control$verbose)
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  trial0 <- logitmodel$trials
  event0 <- logitmodel$event
  names <- logitmodel$names
  term_n <- logitmodel$term_n
  tform <- logitmodel$tform
  keep_constant <- logitmodel$keep_constant
  a_n <- logitmodel$a_n
  modelform <- logitmodel$modelform
  strat_col <- logitmodel$strata
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (length(unique(term_n)) == 1) {
    modelform <- "M"
  } else if (modelform == "GMIX") {
    model_control[["gmix_term"]] <- logitmodel$gmix_term
    model_control[["gmix_theta"]] <- logitmodel$gmix_theta
  }
  if (missing(link)) {
    model_control["logit_odds"] <- TRUE
  } else {
    # "logit_odds", "logit_ident", "logit_loglink"
    link <- tolower(link)
    acceptable <- c("logit_odds", "logit_ident", "logit_loglink", "odds", "ident", "loglink", "id", "odd", "log")
    if (link %in% acceptable) {
      if (link %in% c("logit_odds", "odds", "Odd")) {
        model_control["logit_odds"] <- TRUE
      } else if (link %in% c("logit_ident", "ident", "id")) {
        model_control["logit_ident"] <- TRUE
      } else if (link %in% c("logit_loglink", "loglink", "log")) {
        model_control["logit_loglink"] <- TRUE
      } else {
        stop(gettextf(
          "Error: Argument '%s' not matched to set link options",
          link
        ), domain = NA)
      }
    } else {
      stop(gettextf(
        "Error: Argument '%s' not matched to allowable link options",
        link
      ), domain = NA)
    }
  }
  if (ncol(cons_mat) > 1) {
    model_control["constraint"] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  int_count <- 0.0
  if (any(norm_weight != 1.0)) {
    int_avg_weight <- 0.0
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        int_avg_weight <- int_avg_weight + norm_weight[i]
        int_count <- int_count + 1
      }
    }
    if (int_count > 0) {
      if (control$verbose >= 3) {
        message("Note: Threshold max step adjusted to match new weighting")
      }
      control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunLogisticRegression_Omnibus(df, trial0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, model_control, cons_mat, cons_vec)
  if (int_count > 0) {
    control$thres_step_max <- control$thres_step_max * (int_avg_weight / int_count)
  }
  res$model <- logitmodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight, "tform" = tform), model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  logitres <- new_logitres(res)
  logitres
}

#' Fully runs a case-control regression model, returning the model and results
#'
#' \code{CaseControlRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus matched case-control regression function
#'
#' @param conditional_threshold threshold above which unconditional logistic regression is used to calculate likelihoods in a matched group
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Case Control Wrapper Functions
#' @examples
#' library(data.table)
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
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- CaseCon_Strata(Cancer_Status, e) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CaseControlRun(formula, df,
#'   a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)),
#'   control = control
#' )
CaseControlRun <- function(model, df, a_n = list(c(0)), keep_constant = c(0), control = list(), conditional_threshold = 50, gradient_control = list(), single = FALSE, observed_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), norm = "null", ...) {
  func_t_start <- Sys.time()
  if (is(model, "caseconmodel")) {
    # using already prepped formula and data
    caseconmodel <- copy(model)
    calls <- caseconmodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    #
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    caseconmodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or caseconmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if (!missing(a_n)) {
    caseconmodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    caseconmodel$keep_constant <- keep_constant
  }
  # Checks that the current coxmodel is valid
  validate_caseconsurv(caseconmodel, df)
  if (!caseconmodel$null) {
    caseconmodel <- validate_formula(caseconmodel, df, control$verbose)
  }
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  time1 <- caseconmodel$start_age
  time2 <- caseconmodel$end_age
  event0 <- caseconmodel$event
  names <- caseconmodel$names
  term_n <- caseconmodel$term_n
  tform <- caseconmodel$tform
  keep_constant <- caseconmodel$keep_constant
  a_n <- caseconmodel$a_n
  modelform <- caseconmodel$modelform
  strat_col <- caseconmodel$strata
  cens_weight <- caseconmodel$weight
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (caseconmodel$null) {
    model_control["null"] <- TRUE
    #
    names <- c("CONST")
    term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0)
  } else {
    # check for basic and linear_err
    if (length(unique(term_n)) == 1) {
      modelform <- "M"
    } else if (modelform == "GMIX") {
      model_control[["gmix_term"]] <- caseconmodel$gmix_term
      model_control[["gmix_theta"]] <- caseconmodel$gmix_theta
    }
  }
  if (all(caseconmodel$strata != "NONE")) {
    model_control[["strata"]] <- TRUE
    #
    df$"_strata_col" <- format(df[, strat_col[1], with = FALSE]) # defining a strata column
    for (i in seq_len(length(strat_col) - 1)) {
      df$"_strata_col" <- paste(df$"_strata_col", format(df[, strat_col[i + 1], with = FALSE]), sep = "_") # interacting with any other strata columns
    }
    df$"_strata_col" <- factor(df$"_strata_col") # converting to a factor
  }
  if (time1 != time2) {
    model_control[["time_risk"]] <- TRUE
  }
  #  if (ncol(cons_mat) > 1) {
  #    model_control[["constraint"]] <- TRUE
  #  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  model_control["conditional_threshold"] <- conditional_threshold
  control_def_names <- c(
    "single", "basic", "null", "time_risk",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  int_count <- 0.0
  if (!caseconmodel$null) {
    norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
    a_n <- norm_res$a_n
    cons_mat <- norm_res$cons_mat
    norm_weight <- norm_res$norm_weight
    df <- norm_res$df
    if (any(norm_weight != 1.0)) {
      int_avg_weight <- 0.0
      for (i in seq_along(names)) {
        if (grepl("_int", tform[i])) {
          int_avg_weight <- int_avg_weight + norm_weight[i]
          int_count <- int_count + 1
        }
      }
      if (int_count > 0) {
        if (control$verbose >= 3) {
          message("Note: Threshold max step adjusted to match new weighting")
        }
        control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
      }
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCaseControlRegression_Omnibus(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, "_strata_col", cens_weight, model_control) # , cons_mat, cons_vec)
  if (int_count > 0) {
    control$thres_step_max <- control$thres_step_max * (int_avg_weight / int_count)
  }
  res$model <- caseconmodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (!model_control$null) {
    if (model_control[["constraint"]]) {
      res$constraint_matrix <- cons_mat
      res$constraint_vector <- cons_vec
    }
    res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight, "tform" = tform), model_control)
  }
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  caseconres <- new_caseconres(res)
  caseconres
}

#' Fully runs a joint poisson regression model, returning the model and results
#'
#' \code{PoisRunJoint} uses a list of formula, data.table, and list of controls to prepare and
#' run a Colossus poisson regression function on a joint dataset
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Poisson Wrapper Functions
#' @examples
#' library(data.table)
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
#'   "Ending_Age" = c(30, 45, 57, 47, 36, 60, 55),
#'   "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0),
#'   "Flu_Status" = c(0, 1, 0, 0, 1, 0, 1),
#'   "a" = c(0, 1, 1, 0, 1, 0, 1),
#'   "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
#'   "c" = c(10, 11, 10, 11, 12, 9, 11),
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula_list <- list(Pois(Ending_Age, Cancer_Status) ~ plinear(d, 0),
#'   Pois(Ending_Age, Flu_Status) ~ loglinear(d, 0),
#'   "shared" = Pois(Ending_Age) ~ loglinear(a, b, c, 0)
#' )
#' res <- PoisRunJoint(formula_list, df, control = control)
PoisRunJoint <- function(model, df, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, observed_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), norm = "null", ...) {
  func_t_start <- Sys.time()
  if (is(model, "poismodel")) {
    # using already prepped formula and data
    poismodel <- copy(model)
    calls <- poismodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
  } else if (is.list(model)) {
    # using a list of formula
    res <- get_form_joint(model, df)
    poismodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be list of formula or poismodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if (!missing(a_n)) {
    poismodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    poismodel$keep_constant <- keep_constant
  }
  # Checks that the current coxmodel is valid
  validate_poissurv(poismodel, df)
  poismodel <- validate_formula(poismodel, df, control$verbose)
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  a_n <- poismodel$a_n
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (length(unique(term_n)) == 1) {
    modelform <- "M"
  } else if (modelform == "GMIX") {
    model_control[["gmix_term"]] <- poismodel$gmix_term
    model_control[["gmix_theta"]] <- poismodel$gmix_theta
  }
  if (all(poismodel$strata != "NONE")) {
    model_control["strata"] <- TRUE
  }
  if (ncol(cons_mat) > 1) {
    model_control["constraint"] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  int_count <- 0.0
  if (any(norm_weight != 1.0)) {
    int_avg_weight <- 0.0
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        int_avg_weight <- int_avg_weight + norm_weight[i]
        int_count <- int_count + 1
      }
    }
    if (int_count > 0) {
      if (control$verbose >= 3) {
        message("Note: Threshold max step adjusted to match new weighting")
      }
      control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  if (int_count > 0) {
    control$thres_step_max <- control$thres_step_max * (int_avg_weight / int_count)
  }
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight, "tform" = tform), model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  poisres <- new_poisres(res)
  poisres
}

#' Generic relative risk calculation function
#'
#' \code{RelativeRisk} Generic relative risk calculation function
#' @param x result object from a regression, class coxres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
RelativeRisk <- function(x, df, ...) {
  UseMethod("RelativeRisk", x)
}

#' Generic relative risk calculation function, default option
#'
#' \code{RelativeRisk.default} Generic relative risk calculation function, by default nothing happens
#' @param x result object from a regression, class coxres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
RelativeRisk.default <- function(x, df, ...) {
  return(x)
}

#' Calculates hazard ratios for a reference vector
#'
#' \code{coxres.RelativeRisk} uses a cox result object and data, to evaluate
#' relative risk in the data using the risk model from the result
#'
#' @param x result object from a regression, class coxres
#' @param ... extended to match any future parameters needed
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Cox Analysis Functions
#' @examples
#' library(data.table)
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
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df,
#'   a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)),
#'   control = control
#' )
#' res_risk <- RelativeRisk(res, df)
RelativeRisk.coxres <- function(x, df, a_n = c(), ...) {
  #
  coxmodel <- x$model
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  modelform <- coxmodel$modelform
  #
  calls <- coxmodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  ce <- c(time1, time2, event0)
  val <- Check_Trunc(df, ce)
  if (any(val$ce != ce)) {
    df <- val$df
    ce <- val$ce
    time1 <- ce[1]
    time2 <- ce[1]
  }
  #
  object <- validate_coxres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  control <- object$control
  model_control <- object$modelcontrol
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  res <- Cox_Relative_Risk(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  res
}

#' Performs Cox Proportional Hazard model plots
#'
#' \code{plot.coxres} uses user provided data, time/event columns,
#' vectors specifying the model, and options to choose and save plots
#'
#' @param x result object from a regression, class coxres
#' @param ... can include the named entries for the plot_options parameter
#' @inheritParams R_template
#'
#' @return saves the plots in the current directory and returns the data used for plots
#' @family Plotting Wrapper Functions
#' @export
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
#'   "d" = c(0, 0, 0, 1, 1, 1, 1)
#' )
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df,
#'   control = control,
#'   a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4))
#' )
#' plot_options <- list(
#'   "type" = c("surv", paste(tempfile(),
#'     "run",
#'     sep = ""
#'   )), "studyid" = "UserID",
#'   "verbose" = FALSE
#' )
#' res_plot <- plot(res, df, plot_options)
plot.coxres <- function(x, df, plot_options, a_n = c(), ...) {
  #
  coxmodel <- x$model
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  modelform <- coxmodel$modelform
  strat_col <- coxmodel$strata
  #
  calls <- coxmodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  ce <- c(time1, time2, event0)
  val <- Check_Trunc(df, ce)
  df <- val$df
  if (any(val$ce != ce)) {
    ce <- val$ce
    time1 <- ce[1]
    time2 <- ce[1]
    x$model$start_age <- time1
    x$model$end_age <- time2
  }
  #
  object <- validate_coxres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  control <- object$control
  model_control <- object$modelcontrol
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  #
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- c("verbose", "type", "age_unit", "strat_haz", "strat_col", "martingale", "km", "time_lims", "cov_cols", "studyid") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(plot_options)) {
    plot_options <- extraArgs
  } else if (is.list(plot_options)) {
    plot_options <- c(plot_options, extraArgs)
  } else {
    stop("Error: control argument must be a list")
  }
  if (all(coxmodel$strata != "NONE")) {
    plot_options[["strat_haz"]] <- TRUE
    plot_options$strat_col <- "_strata_col"
    #
    df$"_strata_col" <- format(df[, strat_col[1], with = FALSE]) # defining a strata column
    for (i in seq_len(length(strat_col) - 1)) {
      df$"_strata_col" <- paste(df$"_strata_col", format(df[, strat_col[i + 1], with = FALSE]), sep = "_") # interacting with any other strata columns
    }
    df$"_strata_col" <- factor(df$"_strata_col") # converting to a factor
  }
  res <- RunCoxPlots(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, plot_options = plot_options, model_control = model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  res
}

#' Fully runs a cox or fine-gray regression model with multiple column realizations, returning the model and results
#'
#' \code{CoxRunMulti} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus cox or fine-gray regression function
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "t0" = c(18, 20, 18, 19, 21, 20, 18),
#'   "t1" = c(30, 45, 57, 47, 36, 60, 55),
#'   "lung" = c(0, 0, 1, 0, 1, 0, 0),
#'   "dose" = c(0, 1, 1, 0, 1, 0, 1)
#' )
#' set.seed(3742)
#' df$rand <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
#' names <- c("dose", "rand")
#' realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
#' realization_index <- c("rand")
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiter" = 1,
#'   "halfmax" = 2, "epsilon" = 1e-6,
#'   "deriv_epsilon" = 1e-6, "step_max" = 1.0,
#'   "thres_step_max" = 100.0,
#'   "verbose" = 0, "ties" = "breslow", "double_step" = 1
#' )
#' formula <- Cox(t0, t1, lung) ~ loglinear(dose, rand, 0) + multiplicative()
#' res <- CoxRun(formula, df, control = control)
CoxRunMulti <- function(model, df, a_n = list(c(0)), keep_constant = c(0), realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), gradient_control = list(), single = FALSE, observed_info = FALSE, fma = FALSE, mcml = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...) {
  func_t_start <- Sys.time()
  if (is(model, "coxmodel")) {
    # using already prepped formula and data
    coxmodel <- copy(model)
    calls <- coxmodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    #
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    coxmodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or coxmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if ("CONST" %in% coxmodel$names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (!missing(a_n)) {
    coxmodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    coxmodel$keep_constant <- keep_constant
  }
  # Checks that the current coxmodel is valid
  validate_coxsurv(coxmodel, df)
  if (!coxmodel$null) {
    coxmodel <- validate_formula(coxmodel, df, control$verbose)
  }
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  a_n <- coxmodel$a_n
  modelform <- coxmodel$modelform
  strat_col <- coxmodel$strata
  cens_weight <- coxmodel$weight
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (coxmodel$null) {
    stop()
    model_control["null"] <- TRUE
    #
    names <- c("CONST")
    term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0)
  } else {
    # check for basic and linear_err
    if (length(unique(term_n)) == 1) {
      modelform <- "M"
      if (all(unique(tform) == c("loglin"))) {
        model_control[["basic"]] <- TRUE
      } else if (identical(sort(unique(tform)), c("loglin", "plin")) && sum(tform == "plin") == 1) {
        model_control[["linear_err"]] <- TRUE
      }
    } else if (modelform == "GMIX") {
      model_control[["gmix_term"]] <- coxmodel$gmix_term
      model_control[["gmix_theta"]] <- coxmodel$gmix_theta
    }
  }
  if (all(coxmodel$strata != "NONE")) {
    model_control[["strata"]] <- TRUE
    #
    df$"_strata_col" <- format(df[, strat_col[1], with = FALSE]) # defining a strata column
    for (i in seq_len(length(strat_col) - 1)) {
      df$"_strata_col" <- paste(df$"_strata_col", format(df[, strat_col[i + 1], with = FALSE]), sep = "_") # interacting with any other strata columns
    }
    df$"_strata_col" <- factor(df$"_strata_col") # converting to a factor
  }
  if (coxmodel$weight != "NONE") {
    model_control[["cr"]] <- TRUE
  }
  if (ncol(cons_mat) > 1) {
    model_control[["constraint"]] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  if (fma != mcml) {
    model_control["mcml"] <- mcml
    fma <- !mcml
  } else {
    if (fma) {
      stop("Error: Do not select both fma and mcml, only pick one")
    } else {
      model_control["mcml"] <- mcml
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCoxRegression_Omnibus_Multidose(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = "_strata_col", cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
  res$model <- coxmodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  if (fma) {
    coxres <- new_coxresfma(res)
  } else {
    coxres <- new_coxres(res)
  }
  coxres
}

#' Fully runs a poisson regression model with multiple column realizations, returning the model and results
#'
#' \code{PoisRunMulti} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus poisson regression function
#'
#' @param ... can include the named entries for the control list parameter
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the regression results
#' @export
#' @family Poisson Wrapper Functions
#' @examples
#' library(data.table)
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "t0" = c(18, 20, 18, 19, 21, 20, 18),
#'   "t1" = c(30, 45, 57, 47, 36, 60, 55),
#'   "lung" = c(0, 0, 1, 0, 1, 0, 0),
#'   "dose" = c(0, 1, 1, 0, 1, 0, 1)
#' )
#' set.seed(3742)
#' df$rand <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand0 <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand1 <- floor(runif(nrow(df), min = 0, max = 5))
#' df$rand2 <- floor(runif(nrow(df), min = 0, max = 5))
#' names <- c("dose", "rand")
#' realization_columns <- matrix(c("rand0", "rand1", "rand2"), nrow = 1)
#' realization_index <- c("rand")
#' control <- list(
#'   "ncores" = 1, "lr" = 0.75, "maxiter" = 1,
#'   "halfmax" = 2, "epsilon" = 1e-6,
#'   "deriv_epsilon" = 1e-6, "step_max" = 1.0,
#'   "thres_step_max" = 100.0,
#'   "verbose" = 0, "ties" = "breslow", "double_step" = 1
#' )
#' formula <- Pois(t1, lung) ~ loglinear(CONST, dose, rand, 0) + multiplicative()
#' res <- PoisRun(formula, df, control = control)
PoisRunMulti <- function(model, df, a_n = list(c(0)), keep_constant = c(0), realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), gradient_control = list(), single = FALSE, observed_info = FALSE, fma = FALSE, mcml = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...) {
  func_t_start <- Sys.time()
  if (is(model, "poismodel")) {
    # using already prepped formula and data
    poismodel <- copy(model)
    calls <- poismodel$expres_calls
    df <- ColossusExpressionCall(calls, df)
    #
  } else if (is(model, "formula")) {
    # using a formula class
    res <- get_form(model, df)
    poismodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Error: Incorrect type used for formula, '%s', must be formula or poismodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguments to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(control)) {
    control <- ColossusControl(...)
  } else if (is.list(control)) {
    if (length(extraArgs)) {
      control <- c(control[!(names(control) %in% names(extraArgs))], extraArgs)
    }
    control_args <- intersect(names(control), names(formals(ColossusControl)))
    control <- do.call(ColossusControl, control[control_args])
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  if ("CONST" %in% poismodel$names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (!missing(a_n)) {
    poismodel$a_n <- a_n # assigns the starting parameter values if given
  }
  if (!missing(keep_constant)) { # assigns the paramter constant values if given
    poismodel$keep_constant <- keep_constant
  }
  # Checks that the current poismodel is valid
  validate_poissurv(poismodel, df)
  poismodel <- validate_formula(poismodel, df, control$verbose)
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  # Pull out the actual model vectors and values
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  a_n <- poismodel$a_n
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  # ------------------------------------------------------------------------------ #
  # We want to create the previously used model_control list, based on the input
  model_control <- list()
  if (length(unique(term_n)) == 1) {
    modelform <- "M"
  } else if (modelform == "GMIX") {
    model_control[["gmix_term"]] <- poismodel$gmix_term
    model_control[["gmix_theta"]] <- poismodel$gmix_theta
  }
  if (all(poismodel$strata != "NONE")) {
    model_control["strata"] <- TRUE
  }
  if (ncol(cons_mat) > 1) {
    model_control["constraint"] <- TRUE
  }
  if (!missing(gradient_control)) {
    model_control["gradient"] <- TRUE
    for (nm in names(gradient_control)) {
      model_control[nm] <- gradient_control[nm]
    }
  }
  if (fma != mcml) {
    model_control["mcml"] <- mcml
    fma <- !mcml
  } else {
    if (fma) {
      stop("Error: Do not select both fma and mcml, only pick one")
    } else {
      model_control["mcml"] <- mcml
    }
  }
  model_control["single"] <- single
  model_control["observed_info"] <- observed_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "observed_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunPoisRegression_Omnibus_Multidose(df, pyr0 = pyr0, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = strat_col, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  func_t_end <- Sys.time()
  res$RunTime <- func_t_end - func_t_start
  # ------------------------------------------------------------------------------ #
  if (fma) {
    poisres <- new_poisresfma(res)
  } else {
    poisres <- new_poisres(res)
  }
  poisres
}

#' Generic likelihood boundary calculation function
#'
#' \code{LikelihoodBound} Generic likelihood boundary calculation function
#' @param x result object from a regression, class coxres or poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
LikelihoodBound <- function(x, df, curve_control = list(), control = list(), ...) {
  UseMethod("LikelihoodBound", x)
}

#' Generic likelihood boundary calculation function, default option
#'
#' \code{LikelihoodBound} Generic likelihood boundary calculation function, by default nothing happens
#' @param x result object from a regression, class coxres or poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
LikelihoodBound.default <- function(x, df, curve_control = list(), control = list(), ...) {
  return(x)
}

#' Calculates the likelihood boundary for a completed cox model
#'
#' \code{LikelihoodBound.coxres} solves the confidence interval for a cox model, starting at the optimum point and
#' iteratively optimizing end-points of intervals.
#'
#' @param x result object from a regression, class coxres
#' @param ... can include the named entries for the curve_control list parameter
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
LikelihoodBound.coxres <- function(x, df, curve_control = list(), control = list(), ...) {
  coxmodel <- x$model
  norm <- x$norm
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  modelform <- coxmodel$modelform
  cons_mat <- as.matrix(c(0))
  cons_vec <- c(0)
  strat_col <- coxmodel$strata
  cens_weight <- coxmodel$weight
  #
  calls <- coxmodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  ce <- c(time1, time2, event0)
  val <- Check_Trunc(df, ce)
  df <- val$df
  if (any(val$ce != ce)) {
    ce <- val$ce
    time1 <- ce[1]
    time2 <- ce[1]
    x$model$start_age <- time1
    x$model$end_age <- time2
  }
  #
  object <- validate_coxres(x, df)
  #
  a_n <- object$beta_0
  if (missing(control)) {
    control <- object$control
  }
  #
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control <- do.call(ColossusControl, control[control_args])
  #
  model_control <- object$modelcontrol
  #
  if (model_control[["constraint"]]) {
    cons_mat <- res$constraint_matrix
    cons_vec <- res$constraint_vector
  }
  #
  model_control["log_bound"] <- TRUE
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- c("bisect", "qchi", "para_number", "manual", "search_mult", "maxstep", "step_size") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(curve_control)) {
    curve_control <- extraArgs
  } else if (is.list(curve_control)) {
    curve_control <- c(curve_control, extraArgs)
  } else {
    stop("Error: control argument must be a list")
  }
  #
  model_control <- c(model_control, curve_control)
  if (!("para_number" %in% names(model_control))) {
    model_control["para_number"] <- 1
  } else {
    if (model_control["para_number"] > length(names)) {
      stop("Error: The paranumber used was too large, please use a number between 1 and the number of model elements.")
    } else if (model_control["para_number"] < 1) {
      stop("Error: The paranumber used was less than 1, please use a number between 1 and the number of model elements.")
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  if (all(strat_col != "NONE")) {
    #
    df$"_strata_col" <- format(df[, strat_col[1], with = FALSE]) # defining a strata column
    for (i in seq_len(length(strat_col) - 1)) {
      df$"_strata_col" <- paste(df$"_strata_col", format(df[, strat_col[i + 1], with = FALSE]), sep = "_") # interacting with any other strata columns
    }
    df$"_strata_col" <- factor(df$"_strata_col") # converting to a factor
  }
  #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  if (any(norm_weight != 1.0)) {
    int_avg_weight <- 0.0
    int_count <- 0.0
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        int_avg_weight <- int_avg_weight + norm_weight[i]
        int_count <- int_count + 1
      }
    }
    if (int_count > 0) {
      if (control$verbose >= 3) {
        message("Note: Threshold max step adjusted to match new weighting")
      }
      control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
    }
  }
  #
  if ("bisect" %in% names(model_control)) {
    res <- CoxCurveSolver(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "_strata_col", cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "bisection"
  } else {
    res <- RunCoxRegression_Omnibus(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = "_strata_col", cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "Venzon-Moolgavkar"
  }
  res$model <- coxmodel
  res$beta_0 <- object$beta_0
  res$para_number <- model_control$para_number
  res$coxres <- object
  #
  if (tolower(norm) == "null") {
    # nothing changes
  } else if (tolower(norm) %in% c("max", "mean")) {
    # weight by the maximum value
    res$Parameter_Limits <- res$Parameter_Limits / norm_weight[model_control$para_number]
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        res$Lower_Values[i] <- res$Lower_Values[i] * norm_weight[i]
        res$Upper_Values[i] <- res$Upper_Values[i] * norm_weight[i]
      } else {
        res$Lower_Values[i] <- res$Lower_Values[i] / norm_weight[i]
        res$Upper_Values[i] <- res$Upper_Values[i] / norm_weight[i]
      }
    }
  } else {
    stop(gettextf(
      "Error: Normalization arguement '%s' not valid.",
      norm
    ), domain = NA)
  }
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  coxres <- new_coxresbound(res)
  coxres <- validate_coxresbound(coxres, df)
  coxres
}

#' Calculates the likelihood boundary for a completed Poisson model
#'
#' \code{LikelihoodBound.poisres} solves the confidence interval for a Poisson model, starting at the optimum point and
#' iteratively optimizing end-points of intervals.
#'
#' @param x result object from a regression, class poisres
#' @param ... can include the named entries for the curve_control list parameter
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Poisson Wrapper Functions
LikelihoodBound.poisres <- function(x, df, curve_control = list(), control = list(), ...) {
  poismodel <- x$model
  norm <- x$norm
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  modelform <- poismodel$modelform
  cons_mat <- as.matrix(c(0))
  cons_vec <- c(0)
  strat_col <- poismodel$strata
  #
  calls <- poismodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  object <- validate_poisres(x, df)
  #
  a_n <- object$beta_0
  if (missing(control)) {
    control <- object$control
  }
  #
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control <- do.call(ColossusControl, control[control_args])
  #
  model_control <- object$modelcontrol
  #
  if (model_control[["constraint"]]) {
    cons_mat <- res$constraint_matrix
    cons_vec <- res$constraint_vector
  }
  #
  model_control["log_bound"] <- TRUE
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- c("bisect", "qchi", "para_number", "manual", "search_mult", "maxstep", "step_size") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(curve_control)) {
    curve_control <- extraArgs
  } else if (is.list(curve_control)) {
    curve_control <- c(curve_control, extraArgs)
  } else {
    stop("Error: control argument must be a list")
  }
  #
  model_control <- c(model_control, curve_control)
  if (!("para_number" %in% names(model_control))) {
    model_control["para_number"] <- 1
  } else {
    if (model_control["para_number"] > length(names)) {
      stop("Error: The paranumber used was too large, please use a number between 1 and the number of model elements.")
    } else if (model_control["para_number"] < 1) {
      stop("Error: The paranumber used was less than 1, please use a number between 1 and the number of model elements.")
    }
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat, "tform" = tform), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  if (any(norm_weight != 1.0)) {
    int_avg_weight <- 0.0
    int_count <- 0.0
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        int_avg_weight <- int_avg_weight + norm_weight[i]
        int_count <- int_count + 1
      }
    }
    if (int_count > 0) {
      if (control$verbose >= 3) {
        message("Note: Threshold max step adjusted to match new weighting")
      }
      control$thres_step_max <- control$thres_step_max / (int_avg_weight / int_count)
    }
  }
  #
  if ("bisect" %in% names(model_control)) {
    res <- PoissonCurveSolver(df, pyr0, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "bisection"
  } else {
    res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
    res$method <- "Venzon-Moolgavkar"
  }
  res$model <- poismodel
  res$beta_0 <- object$beta_0
  res$para_number <- model_control$para_number
  res$poisres <- object
  #
  if (tolower(norm) == "null") {
    # nothing changes
  } else if (tolower(norm) %in% c("max", "mean")) {
    # weight by the maximum value
    res$Parameter_Limits <- res$Parameter_Limits / norm_weight[model_control$para_number]
    for (i in seq_along(names)) {
      if (grepl("_int", tform[i])) {
        res$Lower_Values[i] <- res$Lower_Values[i] * norm_weight[i]
        res$Upper_Values[i] <- res$Upper_Values[i] * norm_weight[i]
      } else {
        res$Lower_Values[i] <- res$Lower_Values[i] / norm_weight[i]
        res$Upper_Values[i] <- res$Upper_Values[i] / norm_weight[i]
      }
    }
  } else {
    stop(gettextf(
      "Error: Normalization arguement '%s' not valid.",
      norm
    ), domain = NA)
  }
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  poisres <- new_poisresbound(res)
  poisres <- validate_poisresbound(poisres, df)
  poisres
}


#' Generic background/excess event calculation function
#'
#' \code{EventAssignment} Generic background/excess event calculation function
#' @param x result object from a regression, class poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
EventAssignment <- function(x, df, ...) {
  UseMethod("EventAssignment", x)
}

#' Predicts how many events are due to baseline vs excess
#'
#' \code{EventAssignment} Generic background/excess event calculation function, by default nothing happens
#' @param x result object from a regression, class poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
EventAssignment.default <- function(x, df, ...) {
  return(x)
}

#' Predicts how many events are due to baseline vs excess for a completed poisson model
#'
#' \code{EventAssignment.poisres} uses user provided data, person-year/event columns, vectors specifying the model,
#' and options to calculate background and excess events for a solved Poisson regression
#'
#' @param x result object from a regression, class poisres
#' @param assign_control control list for bounds calculated
#' @param ... can include the named entries for the assign_control list parameter
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Poisson Wrapper Functions
EventAssignment.poisres <- function(x, df, assign_control = list(), control = list(), a_n = c(), ...) {
  poismodel <- x$model
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  #
  calls <- poismodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  object <- validate_poisres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  if (missing(control)) {
    control <- object$control
  }
  #
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control <- do.call(ColossusControl, control[control_args])
  #
  model_control <- object$modelcontrol
  #
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- c("bound", "check_num", "z") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(assign_control)) {
    assign_control <- extraArgs
  } else if (is.list(assign_control)) {
    assign_control <- c(assign_control, extraArgs)
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  check_num <- 1
  z <- 2
  if (length(assign_control) > 0) {
    assign_control$bound <- TRUE
    if ("check_num" %in% names(assign_control)) {
      check_num <- assign_control$check_num
    }
    if ("z" %in% names(assign_control)) {
      z <- assign_control$z
    }
  } else {
    assign_control$bound <- FALSE
  }
  if (!assign_control$bound) {
    # Just a basic event assignment
    res <- RunPoissonEventAssignment(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
  } else {
    a_n <- object$beta_0
    stdev <- object$Standard_Deviation
    #
    # running a boundary solution
    e_mid <- RunPoissonEventAssignment(
      df, pyr0, event0, names, term_n,
      tform, keep_constant, a_n, modelform,
      control, strat_col,
      model_control
    )
    if (length(names) == 1) {
      # There is only one parameter, so we don't need to reoptimize
      a_n <- object$beta_0
      a_n[check_num] <- a_n[check_num] - z * stdev[check_num]
      e_low <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
      a_n <- object$beta_0
      a_n[check_num] <- a_n[check_num] + z * stdev[check_num]
      e_high <- RunPoissonEventAssignment(
        df, pyr0, event0, names,
        term_n, tform, keep_constant,
        a_n, modelform,
        control, strat_col,
        model_control
      )
    } else {
      # We need to shift the parameter, fix it, and then optimize before getting cases
      keep_constant[check_num] <- 1
      #
      model_control <- object$modelcontrol
      norm <- object$norm
      if (!model_control$null) {
        if (model_control[["constraint"]]) {
          cons_mat <- object$constraint_matrix
          cons_vec <- object$constraint_vector
        }
      }
      # Start with low
      a_n <- object$beta_0
      a_n[check_num] <- a_n[check_num] - z * stdev[check_num]
      # Get the new optimum values
      if (model_control[["constraint"]]) {
        low_res <- PoisRun(object, df, control = control, norm = norm, cons_mat = cons_mat, cons_vec = cons_vec, keep_constant = keep_constant, a_n = a_n)
      } else {
        low_res <- PoisRun(object, df, control = control, norm = norm, keep_constant = keep_constant, a_n = a_n)
      }
      a_n <- low_res$beta_0
      e_low <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
      # Now the high
      a_n <- object$beta_0
      a_n[check_num] <- a_n[check_num] + z * stdev[check_num]
      # Get the new optimum values
      if (model_control[["constraint"]]) {
        high_res <- PoisRun(object, df, control = control, norm = norm, cons_mat = cons_mat, cons_vec = cons_vec, keep_constant = keep_constant, a_n = a_n)
      } else {
        high_res <- PoisRun(object, df, control = control, norm = norm, keep_constant = keep_constant, a_n = a_n)
      }
      a_n <- high_res$beta_0
      e_high <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
    }
    res <- list(
      "lower_limit" = e_low, "midpoint" = e_mid,
      "upper_limit" = e_high
    )
    res$parameter_info <- c(names[check_num], tform[check_num], term_n[check_num])
    names(res$parameter_info) <- c("Column", "Subterm", "term_number")
  }
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  res
}

#' Predicts how many events are due to baseline vs excess for a completed poisson likelihood boundary regression
#'
#' \code{EventAssignment.poisresbound} uses user provided data, person-year/event columns, vectors specifying the model,
#' and options to calculate background and excess events for a solved Poisson regression
#'
#' @param x result object from a regression, class poisres
#' @param assign_control control list for bounds calculated
#' @param ... can include the named entries for the assign_control list parameter
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Poisson Wrapper Functions
EventAssignment.poisresbound <- function(x, df, assign_control = list(), control = list(), a_n = c(), ...) {
  poisres <- x$poisres
  poismodel <- poisres$model
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  #
  calls <- poismodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  object <- validate_poisres(poisres, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  if (missing(control)) {
    control <- object$control
  }
  #
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control <- do.call(ColossusControl, control[control_args])
  #
  model_control <- object$modelcontrol
  #
  extraArgs <- list(...) # gather additional arguments
  if (length(extraArgs)) {
    controlargs <- c("bound", "check_num", "z") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Error: Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(assign_control)) {
    assign_control <- extraArgs
  } else if (is.list(assign_control)) {
    assign_control <- c(assign_control, extraArgs)
  } else {
    stop("Error: control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  #
  check_num <- x$para_number
  assign_control$bound <- TRUE
  if (!assign_control$bound) {
    # Just a basic event assignment
    res <- RunPoissonEventAssignment(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
  } else {
    a_n <- object$beta_0
    Parameter_Limits <- x$Parameter_Limits
    #
    # running a boundary solution
    e_mid <- RunPoissonEventAssignment(
      df, pyr0, event0, names, term_n,
      tform, keep_constant, a_n, modelform,
      control, strat_col,
      model_control
    )
    if (length(names) == 1) {
      # There is only one parameter, so we don't need to reoptimize
      a_n <- object$beta_0
      a_n[check_num] <- Parameter_Limits[1]
      e_low <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
      a_n <- object$beta_0
      a_n[check_num] <- Parameter_Limits[2]
      e_high <- RunPoissonEventAssignment(
        df, pyr0, event0, names,
        term_n, tform, keep_constant,
        a_n, modelform,
        control, strat_col,
        model_control
      )
    } else {
      # We need to shift the parameter, fix it, and then optimize before getting cases
      keep_constant[check_num] <- 1
      #
      model_control <- object$modelcontrol
      norm <- object$norm
      if (!model_control$null) {
        if (model_control[["constraint"]]) {
          cons_mat <- object$constraint_matrix
          cons_vec <- object$constraint_vector
        }
      }
      # Start with low
      a_n <- x$Lower_Values
      #      a_n <- object$beta_0
      #      a_n[check_num] <- Parameter_Limits[1]
      # Get the new optimum values
      if (model_control[["constraint"]]) {
        low_res <- PoisRun(object, df, control = control, norm = norm, cons_mat = cons_mat, cons_vec = cons_vec, keep_constant = keep_constant, a_n = a_n)
      } else {
        low_res <- PoisRun(object, df, control = control, norm = norm, keep_constant = keep_constant, a_n = a_n)
      }
      a_n <- low_res$beta_0
      e_low <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
      # Now the high
      a_n <- x$Upper_Values
      #      a_n <- object$beta_0
      #      a_n[check_num] <- Parameter_Limits[2]
      # Get the new optimum values
      if (model_control[["constraint"]]) {
        high_res <- PoisRun(object, df, control = control, norm = norm, cons_mat = cons_mat, cons_vec = cons_vec, keep_constant = keep_constant, a_n = a_n)
      } else {
        high_res <- PoisRun(object, df, control = control, norm = norm, keep_constant = keep_constant, a_n = a_n)
      }
      a_n <- high_res$beta_0
      e_high <- RunPoissonEventAssignment(
        df, pyr0, event0, names, term_n,
        tform, keep_constant, a_n, modelform,
        control, strat_col,
        model_control
      )
    }
    res <- list(
      "lower_limit" = e_low, "midpoint" = e_mid,
      "upper_limit" = e_high
    )
  }
  res$parameter_info <- c(names[check_num], tform[check_num], term_n[check_num])
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  res
}

#' Generic Residual calculation function
#'
#' \code{Residual} Generic Residual calculation function
#' @param x result object from a regression, class coxres or poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
Residual <- function(x, df, ...) {
  UseMethod("Residual", x)
}

#' Generic Residual calculation function, default option
#'
#' \code{Residual.default} Generic Residual calculation function, by default nothing happens
#' @param x result object from a regression, class coxres or poisres
#' @param ... extended for other necessary parameters
#' @inheritParams R_template
#' @export
Residual.default <- function(x, df, ...) {
  return(x)
}

#' Calculates the Residuals for a completed poisson model
#'
#' \code{Residual.poisres} uses user provided data, person-year/event columns, vectors specifying the model,
#' and options to calculate residuals for a solved Poisson regression
#'
#' @param x result object from a regression, class poisres
#' @param pearson boolean to calculate pearson residuals
#' @param deviance boolean to calculate deviance residuals
#' @param ... can include the named entries for the assign_control list parameter
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Poisson Wrapper Functions
Residual.poisres <- function(x, df, control = list(), a_n = c(), pearson = FALSE, deviance = FALSE, ...) {
  poismodel <- x$model
  pyr0 <- poismodel$person_year
  event0 <- poismodel$event
  names <- poismodel$names
  term_n <- poismodel$term_n
  tform <- poismodel$tform
  keep_constant <- poismodel$keep_constant
  modelform <- poismodel$modelform
  strat_col <- poismodel$strata
  #
  calls <- poismodel$expres_calls
  df <- ColossusExpressionCall(calls, df)
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (any(grepl(":intercept", names))) {
    # one of the columns has a :intercept flag
    for (name in names[grepl(":intercept", names)]) {
      if (!(name %in% names(df))) {
        # this isn't a preexisting column
        new_col <- substr(name, 1, nchar(name) - 10)
        df[, name] <- df[, new_col, with = FALSE]
      }
    }
  }
  object <- validate_poisres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  if (missing(control)) {
    control <- object$control
  }
  #
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control <- do.call(ColossusControl, control[control_args])
  #
  model_control <- object$modelcontrol
  #
  if ((pearson == deviance) && (pearson)) {
    stop("Error: Both pearson and deviance cannot be used at once, select only one")
  }
  model_control$pearson <- pearson
  model_control$deviance <- deviance
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  thread_0 <- setDTthreads(control$ncores) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Residual(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
  # ------------------------------------------------------------------------------ #
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  res
}
