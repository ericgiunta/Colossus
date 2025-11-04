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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
  if (coxmodel$strata != "NONE") {
    model_control[["strata"]] <- TRUE
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
  if (!coxmodel$null) {
    norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
    a_n <- norm_res$a_n
    cons_mat <- norm_res$cons_mat
    norm_weight <- norm_res$norm_weight
    df <- norm_res$df
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, cens_weight, model_control, cons_mat, cons_vec)
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
    res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight), model_control)
  }
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
    model_control["gmix_term"] <- coxmodel$gmix_term
    model_control["gmix_theta"] <- coxmodel$gmix_theta
  }
  if (poismodel$strata != "NONE") {
    model_control["strata"] <- TRUE
  }
  if (cons_vec != c(0)) {
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
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight), model_control)
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
    model_control["gmix_term"] <- coxmodel$gmix_term
    model_control["gmix_theta"] <- coxmodel$gmix_theta
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
  if (cons_vec != c(0)) {
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
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  # ------------------------------------------------------------------------------ #
  res <- RunLogisticRegression_Omnibus(df, trial0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, model_control, cons_mat, cons_vec)
  res$model <- logitmodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight), model_control)
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
      model_control["gmix_term"] <- caseconmodel$gmix_term
      model_control["gmix_theta"] <- caseconmodel$gmix_theta
    }
  }
  if (caseconmodel$strata != "NONE") {
    model_control[["strata"]] <- TRUE
  }
  if (time1 != time2) {
    model_control[["time_risk"]] <- TRUE
  }
  if (cons_vec != c(0)) {
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
  if (!caseconmodel$null) {
    norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
    a_n <- norm_res$a_n
    cons_mat <- norm_res$cons_mat
    norm_weight <- norm_res$norm_weight
    df <- norm_res$df
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCaseControlRegression_Omnibus(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, cens_weight, model_control, cons_mat, cons_vec)
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
    res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight), model_control)
  }
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
    model_control["gmix_term"] <- coxmodel$gmix_term
    model_control["gmix_theta"] <- coxmodel$gmix_theta
  }
  if (poismodel$strata != "NONE") {
    model_control["strata"] <- TRUE
  }
  if (cons_vec != c(0)) {
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
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  res$norm <- norm
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  res <- apply_norm(df, norm, names, FALSE, list("output" = res, "norm_weight" = norm_weight), model_control)
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
#'   "halfmax" = 1
#' )
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df,
#'   a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)),
#'   control = control
#' )
#' RelativeRisk(res, df)
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
  Cox_Relative_Risk(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)
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
#'   "ncores" = 2, "lr" = 0.75, "maxiters" = c(1, 1),
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
#' plot(res, df, plot_options)
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
  if (coxmodel$strata != "NONE") {
    plot_options[["strat_haz"]] <- TRUE
    plot_options$strat_col <- strat_col
  }
  RunCoxPlots(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, plot_options = plot_options, model_control = model_control)
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
#'   "ncores" = 2, "lr" = 0.75, "maxiter" = 1,
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
  if (coxmodel$strata != "NONE") {
    model_control[["strata"]] <- TRUE
  }
  if (coxmodel$weight != "NONE") {
    model_control[["cr"]] <- TRUE
  }
  if (cons_vec != c(0)) {
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
  res <- RunCoxRegression_Omnibus_Multidose(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, realization_columns = realization_columns, realization_index = realization_index, control = control, strat_col = strat_col, cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
  res$model <- coxmodel
  res$modelcontrol <- model_control
  res$control <- control
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
#'   "ncores" = 2, "lr" = 0.75, "maxiter" = 1,
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
    model_control["gmix_term"] <- coxmodel$gmix_term
    model_control["gmix_theta"] <- coxmodel$gmix_theta
  }
  if (poismodel$strata != "NONE") {
    model_control["strata"] <- TRUE
  }
  if (cons_vec != c(0)) {
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
  #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
  #
  if ("bisect" %in% names(model_control)) {
    res <- CoxCurveSolver(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "bisection"
  } else {
    res <- RunCoxRegression_Omnibus(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "Venzon-Moolgavkar"
  }
  res$model <- coxmodel
  res$beta_0 <- object$beta_0
  #
  if (tolower(norm) == "null") {
    # nothing changes
  } else if (tolower(norm) %in% c("max", "mean")) {
    # weight by the maximum value
    res$Parameter_Limits <- res$Parameter_Limits / norm_weight[model_control$para_number]
  } else {
    stop(gettextf(
      "Error: Normalization arguement '%s' not valid.",
      norm
    ), domain = NA)
  }
  #
  coxres <- new_coxresbound(res)
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
  object <- validate_poisres(x, df)
  #
  a_n <- object$beta_0
  if (missing(control)) {
    control <- object$control
  }
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
  #
  norm_res <- apply_norm(df, norm, names, TRUE, list("a_n" = a_n, "cons_mat" = cons_mat), model_control)
  a_n <- norm_res$a_n
  cons_mat <- norm_res$cons_mat
  norm_weight <- norm_res$norm_weight
  df <- norm_res$df
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
  #
  if (tolower(norm) == "null") {
    # nothing changes
  } else if (tolower(norm) %in% c("max", "mean")) {
    # weight by the maximum value
    res$Parameter_Limits <- res$Parameter_Limits / norm_weight[model_control$para_number]
  } else {
    stop(gettextf(
      "Error: Normalization arguement '%s' not valid.",
      norm
    ), domain = NA)
  }
  #
  poisres <- new_poisresbound(res)
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
  object <- validate_poisres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  if (missing(control)) {
    control <- object$control
  }
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
  #
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
    # running a boundary solution
    res <- RunPoissonEventAssignment_bound(df, pyr0, event0, x, keep_constant, modelform, check_num, z, control, strat_col, model_control)
  }
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
  object <- validate_poisres(x, df)
  #
  if (missing(a_n)) {
    a_n <- object$beta_0
  }
  if (missing(control)) {
    control <- object$control
  }
  model_control <- object$modelcontrol
  #
  if ((pearson == deviance) && (pearson)) {
    stop("Error: Both pearson and deviance cannot be used at once, select only one")
  }
  model_control$pearson <- pearson
  model_control$deviance <- deviance
  #
  res <- RunPoissonRegression_Residual(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control)
  res
}
