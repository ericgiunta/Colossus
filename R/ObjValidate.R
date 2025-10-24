# CLASS OBJECTS
# ---------------------------------------------------------------------------------------------------------------------------------------------------------- #

new_coxmodel <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "coxmodel"
  )
}

new_coxres <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "coxres"
  )
}

new_coxresfma <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "coxresfma"
  )
}

new_coxresmcml <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "coxresmcml"
  )
}

new_coxresbound <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "coxresbound"
  )
}

new_poismodel <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "poismodel"
  )
}

new_poisres <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "poisres"
  )
}

new_poisresbound <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "poisresbound"
  )
}

new_poisresfma <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "poisresfma"
  )
}

new_poisresmcml <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "poisresmcml"
  )
}

new_caseconmodel <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "caseconmodel"
  )
}

new_caseconres <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "caseconres"
  )
}

new_logitmodel <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "logitmodel"
  )
}

new_logitres <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "logitres"
  )
}


validate_formula <- function(x, df, verbose = FALSE) {
  verbose <- Check_Verbose(verbose)
  if (any(x$term_n != round(x$term_n))) {
    stop(paste("Error: term_n expects integer values, atleast one value was noninteger",
      sep = ""
    ))
  }
  if (min(x$term_n) != 0) {
    if (verbose >= 2) {
      warning(paste("Warning: term_n expects nonnegative integer values and a minimum of 0, minimum value was ",
        min(x$term_n),
        ". Minimum value set to 0, others shifted by ",
        -1 * min(x$term_n),
        sep = ""
      ))
    }
    x$term_n <- x$term_n - min(x$term_n)
  }
  if (length(sort(unique(x$term_n))) != length(min(x$term_n):max(x$term_n))) {
    stop(paste("Error: term_n expects no missing integer values. Term numbers range from ", min(x$term_n),
      " to ", max(x$term_n), " but term_n has ",
      length(unique(x$term_n)), " unique values instead of ",
      length(min(x$term_n):max(x$term_n)),
      sep = ""
    ))
  }
  if (length(x$term_n) < length(x$names)) {
    if (verbose >= 2) {
      warning(paste("Warning: Terms used: ", length(x$term_n),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
    x$term_n <- c(x$term_n, rep(0, length(x$names) -
      length(x$term_n)))
  } else if (length(x$term_n) > length(x$names)) {
    if (verbose >= 2) {
      warning(paste("Warning: Terms used: ", length(x$term_n),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
    x$term_n <- x$term_n[seq_len(length(x$names))]
  }
  # --------------------------------------------------------------------- #
  if (length(x$keep_constant) < length(x$names)) {
    x$keep_constant <- c(x$keep_constant, rep(0, length(x$names) -
      length(x$keep_constant)))
  } else if (length(x$keep_constant) > length(x$names)) {
    x$keep_constant <- x$keep_constant[seq_len(length(x$names))]
  }
  if (min(x$keep_constant) < 0) {
    stop(paste("Error: keep_constant expects 0/1 values, minimum value was ",
      min(x$keep_constant),
      sep = ""
    ))
  }
  if (max(x$keep_constant) > 1) {
    stop(paste("Error: keep_constant expects 0/1 values, maximum value was ",
      max(x$keep_constant),
      sep = ""
    ))
  }
  if (any(x$keep_constant != round(x$keep_constant))) {
    stop(paste("Error: keep_constant expects 0/1 values, atleast one value was noninteger",
      sep = ""
    ))
  }
  # --------------------------------------------------------------------- #
  if (length(x$tform) < length(x$names)) {
    if (verbose >= 2) {
      warning(paste("Warning: Term types used: ", length(x$tform),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
    x$tform <- c(x$tform, rep("loglin", length(x$names) -
      length(x$tform)))
  } else if (length(x$tform) > length(x$names)) {
    if (verbose >= 2) {
      warning(paste("Warning: Term types used: ", length(x$tform),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
    x$tform <- x$tform[seq_len(length(x$names))]
  }
  # --------------------------------------------------------------------- #
  x$tform <- tolower(x$tform)
  for (i in seq_along(x$tform)) {
    if (x$tform[i] %in% c("plin", "plinear", "product-linear")) {
      x$tform[i] <- "plin"
    } else if (x$tform[i] %in% c("lin", "linear")) {
      x$tform[i] <- "lin"
    } else if (x$tform[i] %in% c("loglin", "loglinear", "log-linear")) {
      x$tform[i] <- "loglin"
    }
  }
  tform_order <- c(
    "loglin", "lin", "plin", "loglin_slope", "loglin_top",
    "lin_slope", "lin_int", "quad_slope",
    "step_slope", "step_int",
    "lin_quad_slope", "lin_quad_int", "lin_exp_slope",
    "lin_exp_int", "lin_exp_exp_slope"
  )
  tform_iden <- match(x$tform, tform_order)
  if (any(is.na(tform_iden))) {
    stop(paste("Error: Missing tform option ", x$tform[is.na(tform_iden)],
      ", ",
      sep = ""
    ))
  }
  a <- x$tform
  for (i in seq_len(length(a))) {
    if (i < length(a)) {
      if ((a[i] == "loglin_slope")) {
        if (a[i + 1] != "loglin_top") {
          stop("Error: loglin_top missing")
        }
      } else if ((a[i] == "lin_slope")) {
        if (a[i + 1] != "lin_int") {
          stop("Error: lin_int missing")
        }
      } else if ((a[i] == "step_slope")) {
        if (a[i + 1] != "step_int") {
          stop("Error: step_int missing")
        }
      } else if ((a[i] == "lin_quad_slope")) {
        if (a[i + 1] != "lin_quad_int") {
          stop("Error: lin_quad_int missing")
        }
      } else if ((a[i] == "lin_exp_slope")) {
        if (a[i + 1] != "lin_exp_int") {
          stop("Error: lin_exp_int missing")
        }
      } else if ((a[i] == "lin_exp_int")) {
        if (a[i + 1] != "lin_exp_exp_slope") {
          stop("Error: lin_exp_exp_slope missing")
        }
      }
    }
    if (i > 1) {
      if ((a[i] == "lin_int")) {
        if (a[i - 1] != "lin_slope") {
          stop("Error: lin_slope missing")
        }
      } else if ((a[i] == "step_int")) {
        if (a[i - 1] != "step_slope") {
          stop("Error: step_slope missing")
        }
      } else if ((a[i] == "lin_quad_int")) {
        if (a[i - 1] != "lin_quad_slope") {
          stop("Error: lin_quad_slope missing")
        }
      } else if ((a[i] == "lin_exp_int")) {
        if (a[i - 1] != "lin_exp_slope") {
          stop("Error: lin_exp_slope missing")
        }
      } else if ((a[i] == "lin_exp_exp_slope")) {
        if (a[i - 1] != "lin_exp_int") {
          stop("Error: lin_exp_int missing")
        }
      }
    }
  }
  # --------------------------------------------------------------------- #
  if (((typeof(x$a_n) == "list") && (length(x$a_n) == 1)) || (typeof(x$a_n) != "list")) {
    if (typeof(x$a_n) == "list") {
      x$a_n <- x$a_n[[1]]
    }
    if (length(x$a_n) < length(x$names)) {
      if (verbose >= 2) {
        warning(paste("Warning: Parameters used: ",
          length(x$a_n), ", Covariates used: ",
          length(x$names), ", Remaining filled with 0.01",
          sep = ""
        ))
      }
      x$a_n <- c(x$a_n, rep(0.01, length(x$names) - length(x$a_n)))
    } else if (length(x$a_n) > length(x$names)) {
      stop(paste("Error: Parameters used: ", length(x$a_n),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
  } else {
    a_0 <- x$a_n[[1]]
    for (a_i in x$a_n) {
      if (length(a_i) != length(a_0)) {
        stop(paste("Error: Parameters used in first option: ",
          length(a_0),
          ", Parameters used in different option: ",
          length(a_i),
          ", please fix parameter length",
          sep = ""
        ))
      }
    }
    if (length(a_0) < length(x$names)) {
      if (verbose >= 2) {
        warning(paste("Warning: Parameters used: ", length(a_0),
          ", Covariates used: ", length(x$names),
          ", Remaining filled with 0.01",
          sep = ""
        ))
      }
      for (i in seq_len(length(x$a_n))) {
        x$a_n[[i]] <- c(
          x$a_n[[i]],
          rep(0.01, length(x$names) - length(x$a_n[[i]]))
        )
      }
    } else if (length(a_0) > length(x$names)) {
      stop(paste("Error: Parameters used: ", length(a_0),
        ", Covariates used: ", length(x$names),
        sep = ""
      ))
    }
  }
  # --------------------------------------------------------------------- #
  name_check <- unique(x$names)
  if (!all(name_check %in% names(df))) {
    stop("Error: Atleast one model covariate not in the data")
  }
  # --------------------------------------------------------------------- #
  term_tot <- max(x$term_n) + 1
  x$modelform <- toupper(x$modelform)
  acceptable <- toupper(c(
    "multiplicative", "multiplicative-excess", "additive", "product-additive",
    "product-additive-excess", "a", "pa", "pae", "m", "me",
    "gmix", "geometric-mixture", "gmix-r", "relative-geometric-mixture",
    "gmix-e", "excess-geometric-mixture"
  ))
  if (x$modelform %in% acceptable) {
    if (x$modelform %in% c("ME", "MULTIPLICATIVE", "MULTIPLICATIVE-EXCESS")) {
      x$modelform <- "M"
    } else if (x$modelform == "ADDITIVE") {
      x$modelform <- "A"
    } else if (x$modelform == "PRODUCT-ADDITIVE") {
      x$modelform <- "PA"
    } else if (x$modelform == "PRODUCT-ADDITIVE-EXCESS") {
      x$modelform <- "PAE"
    } else if (x$modelform %in% c("GMIX-R", "RELATIVE-GEOMETRIC-MIXTURE")) {
      x$gmix_term <- rep(0, term_tot)
      x$modelform <- "GMIX"
    } else if (x$modelform %in% c("GMIX-E", "EXCESS-GEOMETRIC-MIXTURE")) {
      x$gmix_term <- rep(1, term_tot)
      x$modelform <- "GMIX"
    }
  } else {
    stop(paste("Error: Model formula ", x$modelform,
      " not in acceptable list",
      sep = ""
    ))
  }
  if (x$modelform == "GMIX") {
    if (length(x$gmix_term) != term_tot) {
      stop(paste("Error: Terms used:", term_tot,
        ", Terms with gmix types available:",
        length(x$gmix_term),
        sep = " "
      ))
    }
  }
  # --------------------------------------------------------------------- #
  x
}

validate_coxsurv <- function(x, df) {
  if (!is(x, "coxmodel")) {
    stop("Error: Non cox formula used in cox regression")
  }
  if (x$start_age == x$end_age) {
    stop(
      "Error: The starting and ending interval times were set to the same column, they must be different"
    )
  }
  if (x$event == "") {
    stop(
      "Error: The event column must not be empty"
    )
  }
  if (!(x$start_age %in% names(df))) {
    stop(
      "Error: Interval start column not in the data"
    )
  }
  if (!(x$end_age %in% names(df))) {
    stop(
      "Error: Interval end column not in the data"
    )
  }
  if (!(x$event %in% names(df))) {
    stop(
      "Error: Event column not in the data"
    )
  }
}

validate_poissurv <- function(x, df) {
  if (!is(x, "poismodel")) {
    stop("Error: Non Poisson formula used in Poisson regression")
  }
  if (x$event == "") {
    stop(
      "Error: The event column must not be empty"
    )
  }
  if (!(x$person_year %in% names(df))) {
    stop(
      "Error: Person-Year column not in the data"
    )
  }
  if (!(x$event %in% names(df))) {
    stop(
      "Error: Event column not in the data"
    )
  }
}

validate_caseconsurv <- function(x, df) {
  if (!is(x, "caseconmodel")) {
    stop("Error: Non Case-Control formula used in Case_Control regression")
  }
  if (x$start_age == x$end_age) {
    if (x$start != "NONE") {
      stop(
        "Error: The starting and ending interval times were set to the same column, they must be different or both '%trunc%'"
      )
    }
  } else {
    if (!(x$start_age %in% names(df))) {
      stop(
        "Error: Interval start column not in the data"
      )
    }
    if (!(x$end_age %in% names(df))) {
      stop(
        "Error: Interval end column not in the data"
      )
    }
  }
  if (x$event == "") {
    stop(
      "Error: The event column must not be empty"
    )
  }
  if (!(x$event %in% names(df))) {
    stop(
      "Error: Event column not in the data"
    )
  }
}

validate_logitsurv <- function(x, df) {
  if (!is(x, "logitmodel")) {
    stop("Error: Non logistic formula used in logistic regression")
  }
  if (x$event == "") {
    stop(
      "Error: The event column must not be empty"
    )
  }
  if (x$trials == "") {
    stop(
      "Error: The trials column must not be empty"
    )
  }
  if (!(x$trials %in% names(df))) {
    stop(
      "Error: Interval start column not in the data"
    )
  }
  if (!(x$event %in% names(df))) {
    stop(
      "Error: Event column not in the data"
    )
  }
  if (any(df[, x$event, with = FALSE] > df[, x$trials, with = FALSE])) {
    stop("Error: In atleast one row, the number of events was larger than the number of trials")
  }
}

validate_coxres <- function(x, df) {
  coxmodel <- x$model
  null <- coxmodel$null
  control <- x$control
  verbose <- control$verbose
  #
  validate_coxsurv(coxmodel, df)
  if (!null) {
    coxmodel <- validate_formula(coxmodel, df, verbose)
  }
  x$model <- coxmodel
  x
}

validate_poisres <- function(x, df) {
  poismodel <- x$model
  control <- x$control
  verbose <- control$verbose
  #
  validate_poissurv(poismodel, df)
  poismodel <- validate_formula(poismodel, df, verbose)
  x$model <- poismodel
  x
}

coxmodel <- function(start_age = "",
                     end_age = "",
                     event = "",
                     strata = "",
                     weight = "",
                     null = FALSE,
                     term_n = c(),
                     tform = c(),
                     names = c(),
                     modelform = "",
                     gmix_term = c(),
                     gmix_theta = 0,
                     a_n = list(),
                     keep_constant = c(),
                     df = data.table(),
                     expres_calls = list(),
                     verbose = FALSE) {
  cox_obj <- list(
    "start_age" = start_age, "end_age" = end_age, "event" = event, "strata" = strata, "weight" = weight, "null" = null,
    "term_n" = term_n, "tform" = tform, "names" = names, "a_n" = a_n, "keep_constant" = keep_constant, "modelform" = modelform,
    "gmix_term" = gmix_term, "gmix_theta" = gmix_theta, "expres_calls" = expres_calls
  )
  cox_obj <- new_coxmodel(cox_obj)
  validate_coxsurv(cox_obj, df)
  if (!null) {
    cox_obj <- validate_formula(cox_obj, df, verbose)
  }
  cox_obj
}

poismodel <- function(person_year = "",
                      event = "",
                      strata = "",
                      term_n = c(),
                      tform = c(),
                      names = c(),
                      modelform = "",
                      gmix_term = c(),
                      gmix_theta = 0,
                      a_n = list(),
                      keep_constant = c(),
                      df = data.table(),
                      expres_calls = list(),
                      verbose = FALSE) {
  pois_obj <- list(
    "person_year" = person_year, "event" = event, "strata" = strata,
    "term_n" = term_n, "tform" = tform, "names" = names, "a_n" = a_n, "keep_constant" = keep_constant, "modelform" = modelform,
    "gmix_term" = gmix_term, "gmix_theta" = gmix_theta, "expres_calls" = expres_calls
  )
  pois_obj <- new_poismodel(pois_obj)
  validate_poissurv(pois_obj, df)
  pois_obj <- validate_formula(pois_obj, df, verbose)
  pois_obj
}

caseconmodel <- function(start_age = "",
                         end_age = "",
                         event = "",
                         strata = "",
                         null = FALSE,
                         term_n = c(),
                         tform = c(),
                         names = c(),
                         modelform = "",
                         gmix_term = c(),
                         gmix_theta = 0,
                         a_n = list(),
                         keep_constant = c(),
                         df = data.table(),
                         expres_calls = list(),
                         verbose = FALSE) {
  casecon_obj <- list(
    "start_age" = start_age, "end_age" = end_age, "event" = event, "strata" = strata, "null" = null,
    "term_n" = term_n, "tform" = tform, "names" = names, "a_n" = a_n, "keep_constant" = keep_constant, "modelform" = modelform,
    "gmix_term" = gmix_term, "gmix_theta" = gmix_theta, "expres_calls" = expres_calls
  )
  casecon_obj <- new_caseconmodel(casecon_obj)
  validate_caseconsurv(casecon_obj, df)
  if (!null) {
    casecon_obj <- validate_formula(casecon_obj, df, verbose)
  }
  casecon_obj
}

logitmodel <- function(trials = "",
                       event = "",
                       strata = "",
                       term_n = c(),
                       tform = c(),
                       names = c(),
                       modelform = "",
                       gmix_term = c(),
                       gmix_theta = 0,
                       a_n = list(),
                       keep_constant = c(),
                       df = data.table(),
                       expres_calls = list(),
                       verbose = FALSE) {
  logit_obj <- list(
    "trials" = trials, "event" = event, "strata" = strata,
    "term_n" = term_n, "tform" = tform, "names" = names, "a_n" = a_n, "keep_constant" = keep_constant, "modelform" = modelform,
    "gmix_term" = gmix_term, "gmix_theta" = gmix_theta, "expres_calls" = expres_calls
  )
  logit_obj <- new_logitmodel(logit_obj)
  validate_logitsurv(logit_obj, df)
  logit_obj <- validate_formula(logit_obj, df, verbose)
  logit_obj
}

# ------------------------------------------------------------------------ #
ColossusControl <- function(verbose = 0,
                            lr = 0.75,
                            maxiter = 20,
                            halfmax = 5,
                            epsilon = 1e-4,
                            deriv_epsilon = 1e-4,
                            step_max = 1.0,
                            thres_step_max = 1.0,
                            ties = "breslow",
                            ncores = as.numeric(detectCores())) {
  verbose <- Check_Verbose(verbose)
  control <- list(
    "verbose" = verbose,
    "lr" = lr,
    "maxiter" = maxiter,
    "halfmax" = halfmax,
    "epsilon" = epsilon,
    "deriv_epsilon" = deriv_epsilon,
    "step_max" = step_max,
    "thres_step_max" = thres_step_max,
    "ties" = ties,
    "ncores" = ncores
  )
  control_def <- list(
    "verbose" = 0, "lr" = 0.75, "maxiter" = 20,
    "halfmax" = 5, "epsilon" = 1e-4,
    "deriv_epsilon" = 1e-4, "step_max" = 1.0,
    "thres_step_max" = 100.0,
    "ties" = "breslow",
    "ncores" = as.numeric(detectCores())
  )
  #
  names(control) <- tolower(names(control))
  for (nm in names(control_def)) {
    if (nm %in% names(control)) {
      if (nm == "ncores") {
        if (control$ncores > control_def$ncores) {
          stop(paste("Error: Cores Requested:", control["ncores"],
            ", Cores Available:", control_def["ncores"],
            sep = " "
          )) # nocov
        }
      } else if (nm == "verbose") {
        control$verbose <- Check_Verbose(control$verbose)
      }
    } else {
      control[nm] <- control_def[nm]
    }
  }
  control["ties"] <- tolower(control["ties"])
  control_min <- list(
    "verbose" = 0, "lr" = 0.0, "maxiter" = -1,
    "halfmax" = 0, "epsilon" = 0.0,
    "deriv_epsilon" = 0.0, "step_max" = 0.0,
    "thres_step_max" = 0.0
  )
  for (nm in names(control_min)) {
    if (control[[nm]] < control_min[[nm]]) {
      control[nm] <- control_min[nm]
    }
  }
  control_int <- list(
    "verbose" = 0, "maxiter" = -1,
    "halfmax" = 0
  )
  for (nm in names(control_int)) {
    control[nm] <- as.integer(control[nm])
  }
  #
  control
}
