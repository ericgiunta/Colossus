#' Fully runs a cox or fine-gray regression model, returning the model and results
#'
#' \code{CoxRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus cox or fine-gray regression function
#'
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
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df, a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)))
CoxRun <- function(model, table, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, oberved_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...) {
  if (class(model) == "coxmodel") {
    # using already prepped formula and data
    coxmodel <- copy(model)
    df <- copy(table)
    #
  } else if (class(model) == "formula") {
    # using a formula class
    res <- get_form(model, table)
    coxmodel <- res$model
    df <- res$data
    # } else if (class(formula) == "function") {
    #   # using a formula class
    #   res <- get_form(as.formula(formula), table)
    #   coxmodel <- res$model
    #   df <- res$data
  } else {
    print(model)
    stop(gettextf(
      "Incorrect type used for formula, '%s', must be formula or coxmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguements to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
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
    stop("control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
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
      } else if (all(sort(unique(tform)) == c("loglin", "plin")) && sum(tform == "plin") == 1) {
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
  model_control["single"] <- single
  model_control["oberved_info"] <- oberved_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "oberved_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, cens_weight, model_control, cons_mat, cons_vec)
  res$model <- coxmodel
  res$modelcontrol <- model_control
  res$control <- control
  if (model_control[["constraint"]]) {
    res$constraint_matrix <- cons_mat
    res$constraint_vector <- cons_vec
  }
  # ------------------------------------------------------------------------------ #
  coxres <- new_coxres(res)
  coxres
}

#' Fully runs a poisson regression model, returning the model and results
#'
#' \code{PoisRun} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus poisson regression function
#'
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
#' formula <- Pois(Ending_Age, Cancer_Status) ~ loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- PoisRun(formula, df, a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)))
PoisRun <- function(model, table, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, oberved_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...) {
  if (class(model) == "poismodel") {
    # using already prepped formula and data
    poismodel <- copy(model)
    df <- copy(table)
  } else if (class(model) == "formula") {
    # using a formula class
    res <- get_form(model, table)
    poismodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Incorrect type used for formula, '%s', must be formula or poismodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguements to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
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
    stop("control argument must be a list")
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
  if (poismodel$null) {
    stop("Null model is not valid for poisson regression")
  } else {
    # check for basic and linear_err
    if (length(unique(term_n)) == 1) {
      modelform <- "M"
    } else if (modelform == "GMIX") {
      model_control["gmix_term"] <- coxmodel$gmix_term
      model_control["gmix_theta"] <- coxmodel$gmix_theta
    }
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
  model_control["oberved_info"] <- oberved_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "oberved_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  poisres <- new_poisres(res)
  poisres
}

#' Fully runs a joint poisson regression model, returning the model and results
#'
#' \code{PoisRunJoint} uses a list of formula, data.table, and list of controls to prepare and
#' run a Colossus poisson regression function on a joint dataset
#'
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
#'   "Flu_Status"    = c(0, 1, 0, 0, 1, 0, 1),
#'   "a" = c(0, 1, 1, 0, 1, 0, 1),
#'   "b" = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
#'   "c" = c(10, 11, 10, 11, 12, 9, 11),
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' formula_list <- list(Pois(Ending_Age, Cancer_Status) ~ plinear(d, 0),
#'                      Pois(Ending_Age, Flu_Status)    ~ loglinear(d, 0),
#'                  "shared" = Pois(Ending_Age)    ~ loglinear(a, b, c, 0)
#' )
#' res <- PoisRunJoint(formula_list, df)
PoisRunJoint <- function(model, table, a_n = list(c(0)), keep_constant = c(0), control = list(), gradient_control = list(), single = FALSE, oberved_info = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...){
  if (class(model) == "poismodel") {
    # using already prepped formula and data
    poismodel <- copy(model)
    df <- copy(table)
  } else if (is.list(model)) {
    # using a list of formula
    res <- get_form_joint(model, table)
    poismodel <- res$model
    df <- res$data
  } else {
    stop(gettextf(
      "Incorrect type used for formula, '%s', must be list of formula or poismodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguements to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
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
    stop("control argument must be a list")
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
  if (poismodel$null) {
    stop("Null model is not valid for poisson regression")
  } else {
    # check for basic and linear_err
    if (length(unique(term_n)) == 1) {
      modelform <- "M"
    } else if (modelform == "GMIX") {
      model_control["gmix_term"] <- coxmodel$gmix_term
      model_control["gmix_theta"] <- coxmodel$gmix_theta
    }
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
  model_control["oberved_info"] <- oberved_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "oberved_info"
  )
  for (nm in control_def_names) {
    if (!(nm %in% names(model_control))) {
      model_control[nm] <- FALSE
    }
  }
  # ------------------------------------------------------------------------------ #
  res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
  res$model <- poismodel
  res$modelcontrol <- model_control
  res$control <- control
  # ------------------------------------------------------------------------------ #
  poisres <- new_poisres(res)
  poisres
}

#' @export
RelativeRisk <- function(object, df) {
  UseMethod("RelativeRisk", object)
}

#' @export
RelativeRisk.default <- function(object, df) {
  return(object)
}

#' Calculates hazard ratios for a reference vector
#'
#' \code{coxres.RelativeRisk} uses a cox result object and data, to evaluate
#' relative risk in the data using the risk model from the result
#'
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
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df, a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)))
#' RelativeRisk(res, df)
RelativeRisk.coxres <- function(object, df) {
  #
  coxmodel <- object$model
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  modelform <- coxmodel$modelform
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
  if (any(val$ce != ce)){
      df <- val$df
      ce <- val$ce
      time1 <- ce[1]
      time2 <- ce[1]
  }
  #
  object <- validate_coxres(object, df)
  #
  a_n <- object$beta_0
  control <- object$control
  model_control <- object$modelcontrol
  # Cox_Relative_Risk <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), model_control = list())
  Cox_Relative_Risk(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, model_control = model_control)
}

#' Performs Cox Proportional Hazard model plots
#'
#' \code{plot.coxres} uses user provided data, time/event columns,
#' vectors specifying the model, and options to choose and save plots
#'
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
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' res <- CoxRun(formula, df, a_n = list(c(1.1, -0.1, 0.2, 0.5), c(1.6, -0.12, 0.3, 0.4)))
#' plot_options <- list(
#'   "type" = c("surv", paste(tempfile(),
#'     "run",
#'     sep = ""
#'   )), "studyid" = "UserID",
#'   "verbose" = FALSE
#' )
#' plot(res, df, plot_options)
plot.coxres <- function(object, df, plot_options, ...) {
  #
  coxmodel <- object$model
  time1 <- coxmodel$start_age
  time2 <- coxmodel$end_age
  event0 <- coxmodel$event
  names <- coxmodel$names
  term_n <- coxmodel$term_n
  tform <- coxmodel$tform
  keep_constant <- coxmodel$keep_constant
  modelform <- coxmodel$modelform
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
  if (any(val$ce != ce)){
      df <- val$df
      ce <- val$ce
      time1 <- ce[1]
      time2 <- ce[1]
  }
  #
  object <- validate_coxres(object, df)
  #
  a_n <- object$beta_0
  control <- object$control
  model_control <- object$modelcontrol
  #
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- c("verbose", "type", "age_unit", "strat_haz", "strat_col", "martingale", "km", "time_lims", "cov_cols", "studyid") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(plot_options)) {
    plot_options <- extraArgs
  } else if (is.list(plot_options)) {
    plot_options <- c(plot_options, extraArgs)
  } else {
    stop("control argument must be a list")
  }
  # RunCoxPlots(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), plot_options = list(), model_control = list())
  RunCoxPlots(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, plot_options = plot_options, model_control = model_control)
}

#' Fully runs a cox or fine-gray regression model with multiple column realizations, returning the model and results
#'
#' \code{CoxRunMulti} uses a formula, data.table, and list of controls to prepare and
#' run a Colossus cox or fine-gray regression function
#'
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
#' res <- CoxRun(formula, df, a_n = a_n, control = control)
CoxRunMulti <- function(model, table, a_n = list(c(0)), keep_constant = c(0), realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), gradient_control = list(), single = FALSE, oberved_info = FALSE, fma = FALSE, mcml = FALSE, cons_mat = as.matrix(c(0)), cons_vec = c(0), ...) {
  if (class(model) == "coxmodel") {
    # using already prepped formula and data
    coxmodel <- copy(model)
    df <- copy(table)
    #
  } else if (class(model) == "formula") {
    # using a formula class
    res <- get_form(model, table)
    coxmodel <- res$model
    df <- res$data
  } else {
    print(model)
    stop(gettextf(
      "Incorrect type used for formula, '%s', must be formula or coxmodel class",
      class(model)
    ))
  }
  # ------------------------------------------------------------------------------ #
  # we want to let the user add in control arguements to their call
  # code copied from survival/R/coxph.R github and modified for our purpose
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- names(formals(ColossusControl)) # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
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
    stop("control argument must be a list")
  }
  # ------------------------------------------------------------------------------ #
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
      } else if (all(sort(unique(tform)) == c("loglin", "plin")) && sum(tform == "plin") == 1) {
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
  } else {
    if (fma) {
      stop("Do not select both fma and mcml, only pick one")
    } else {
      model_control["mcml"] <- mcml
    }
  }
  model_control["single"] <- single
  model_control["oberved_info"] <- oberved_info
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "oberved_info"
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
  # ------------------------------------------------------------------------------ #
  if (fma) {
    coxres <- new_coxresfma(res)
  } else {
    coxres <- new_coxres(res)
  }
  coxres
}

#' @export
LikelihoodBound <- function(object, df, curve_control = list(), control = list()) {
  UseMethod("LikelihoodBound", object)
}

#' @export
LikelihoodBound.default <- function(object, df, curve_control = list(), control = list()) {
  return(object)
}

#' Calculates the likelihood boundary for a completed cox model
#'
#' \code{LikelihoodBound.coxres} solves the confidence interval for a cox model, starting at the optimum point and
#' iteratively optimizing end-points of intervals.
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
LikelihoodBound.coxres <- function(object, df, curve_control = list(), control = list(), ...) {
  coxmodel <- object$model
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
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  ce <- c(time1, time2, event0)
  val <- Check_Trunc(df, ce)
  if (any(val$ce != ce)){
      df <- val$df
      ce <- val$ce
      time1 <- ce[1]
      time2 <- ce[1]
  }
  #
  object <- validate_coxres(object, df)
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
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- c("bisect", "qchi", "para_number", "manual", "search_mult", "maxstep", "step_size") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(curve_control)) {
    curve_control <- extraArgs
  } else if (is.list(curve_control)) {
    curve_control <- c(curve_control, extraArgs)
  } else {
    stop("control argument must be a list")
  }
  #
  model_control <- c(model_control, curve_control)
  #
  if ("bisect" %in% names(model_control)) {
    res <- CoxCurveSolver(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "bisection"
  } else {
    res <- RunCoxRegression_Omnibus(df, time1 = time1, time2 = time2, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, cens_weight = cens_weight, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "Venzon-Moolgavkar"
  }
  coxres <- new_coxresbound(res)
  coxres
}

#' Calculates the likelihood boundary for a completed Poisson model
#'
#' \code{LikelihoodBound.poisres} solves the confidence interval for a Poisson model, starting at the optimum point and
#' iteratively optimizing end-points of intervals.
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
LikelihoodBound.poisres <- function(object, df, curve_control = list(), control = list(), ...) {
  poismodel <- object$model
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
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  object <- validate_poisres(object, df)
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
  extraArgs <- list(...) # gather additional arguements
  if (length(extraArgs)) {
    controlargs <- c("bisect", "qchi", "para_number", "manual", "search_mult", "maxstep", "step_size") # names used in control function
    indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L) # check for any mismatched names
    if (any(indx == 0L)) {
      stop(gettextf(
        "Argument '%s' not matched",
        names(extraArgs)[indx == 0L]
      ), domain = NA)
    }
  }
  if (missing(curve_control)) {
    curve_control <- extraArgs
  } else if (is.list(curve_control)) {
    curve_control <- c(curve_control, extraArgs)
  } else {
    stop("control argument must be a list")
  }
  #
  model_control <- c(model_control, curve_control)
  #
  if ("bisect" %in% names(model_control)) {
    res <- PoissonCurveSolver(df, pyr0, event0 = event0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n, modelform = modelform, control = control, strat_col = strat_col, model_control = model_control, cons_mat = cons_mat, cons_vec = cons_vec)
    res$method <- "bisection"
  } else {
    res <- RunPoissonRegression_Omnibus(df, pyr0, event0, names, term_n, tform, keep_constant, a_n, modelform, control, strat_col, model_control, cons_mat, cons_vec)
    res$method <- "Venzon-Moolgavkar"
  }
  poisres <- new_poisresbound(res)
  poisres
}
