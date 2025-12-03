#' Interprets a Poisson joint formula and makes necessary changes to data
#'
#' \code{get_form_joint} uses two event formula, a shared formula,
#' and data.table, to fully describe the model for a joint Poisson model.
#'
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the updated data
#' @family Formula Interpretation
#' @export
get_form_joint <- function(formula_list, df) {
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
  if (is.list(formula_list)) {
    if ("shared" %in% names(formula_list)) {
      formula_shared <- formula_list$shared
      formula_list <- formula_list[names(formula_list) != "shared"]
    } else {
      formula_shared <- .~.
    }
  } else if (is(formula_list, "formula")) {
    formula_list <- list(formula_list)
    formula_shared <- .~.
  } else {
    stop("Error: Joint formula list wasn't a list or a formula")
  }
  model_list <- list()
  for (formula_i in seq_along(formula_list)) {
    #
    formula <- formula_list[[formula_i]]
    surv_obj <- format(formula[[2]])
    model_obj <- paste(format(formula[[3]]), collapse = " ")
    surv_obj <- gsub(" ", "", surv_obj)
    model_obj <- gsub(" ", "", model_obj)
    res <- get_form_list(surv_obj, model_obj, df)
    model <- res$model
    df <- res$data
    if (!grepl("pois", model$surv_model_type)) {
      stop("Error: Atleast one model was not a poisson model")
    }
    model_temp <- list(model)
    model_list <- c(model_list, model_temp)
  }
  #
  model_obj <- paste(format(formula_shared[[3]]), collapse = " ")
  if (model_obj != ".") {
    res <- get_form_risk(model_obj, df)
    model_share <- res$model
    df <- res$data
  } else {
    model_share <- list(tform = c(), term_n = c(), names = c(), modelform = "M", gmix_theta = 0, gmix_term = c())
  }
  #
  # Now we need to check that the values are similar
  model_1 <- model_list[[1]]
  for (model_i in 2:length(model_list)) {
    model_2 <- model_list[[model_i]]
    # The pyr should all be the same
    if (model_1$pyr != model_2$pyr) {
      stop(paste("Error: The joint models need to use the same person-year column. Instead they use ", model_1$pyr, " and ", model_2$pyr, ".", sep = ""))
    }
    # The strata should match
    if (model_1$strata != model_2$strata) {
      stop(paste("Error: The joint models need to use the same stratification.", sep = ""))
    }
    # The modelform should match
    if (model_1$modelform != model_2$modelform) {
      stop(paste("Error: The joint models need to use the same modelform. Instead they use ", model_1$modelform, " and ", model_2$modelform, ".", sep = ""))
    }
    if (model_1$gmix_theta != model_2$gmix_theta) {
      stop(paste("Error: The joint models need to use the same geometric mixture theta value. Instead they use ", model_1$gmix_theta, " and ", model_2$gmix_theta, ".", sep = ""))
    }
  }
  if (length(model_share$tform) != 0) {
    if (model_1$modelform != model_share$modelform) {
      stop(paste("Error: The joint models and the shared model need to use the same modelform. Instead they use ", model_1$modelform, " and ", model_share$modelform, ".", sep = ""))
    }
    if (model_1$gmix_theta != model_share$gmix_theta) {
      stop(paste("Error: The joint models and the shared model need to use the same geometric mixture theta value. Instead they use ", model_1$gmix_theta, " and ", model_share$gmix_theta, ".", sep = ""))
    }
  } else {
    model_share$modelform <- model_1$modelform
    model_share$gmix_theta <- model_1$gmix_theta
  }
  # Now we pull everything out and get the joint model built
  # Pull it all out
  events <- c()
  name_list <- list(
    "shared" = model_share$names
  )
  term_n_list <- list(
    "shared" = model_share$term_n
  )
  tform_list <- list(
    "shared" = model_share$tform
  )
  for (model_i in seq_along(model_list)) {
    model <- model_list[[model_i]]
    event_temp <- model$event
    name_temp <- list()
    termn_temp <- list()
    tform_temp <- list()
    name_temp[[event_temp]] <- model$names
    termn_temp[[event_temp]] <- model$term_n
    tform_temp[[event_temp]] <- model$tform
    #
    events <- c(events, event_temp)
    term_n_list <- c(term_n_list, termn_temp)
    name_list <- c(name_list, name_temp)
    tform_list <- c(tform_list, tform_temp)
  }
  keep_constant_list <- list()
  a_n_list <- list()
  # Combine it
  val <- Joint_Multiple_Events(
    df, events, name_list,
    term_n_list, tform_list,
    keep_constant_list, a_n_list
  )
  df <- val$df
  names <- val$names
  term_n <- val$term_n
  tform <- val$tform
  keep_constant <- val$keep_constant
  a_n <- val$a_n
  #
  pyr <- model_1$pyr
  strata <- model_1$strata
  modelform <- model_1$modelform
  #
  gmix_term <- model_1$gmix_term
  gmix_theta <- model_1$gmix_theta
  if (modelform == "GMIX") {
    # We need to combine the results
    gmix_term_1 <- model_1$gmix_term
    gmix_term_2 <- model_2$gmix_term
    gmix_term_s <- model_share$gmix_term
    #
    len_1 <- length(gmix_term_1)
    len_2 <- length(gmix_term_2)
    len_s <- length(gmix_term_s)
    if (len_1 < len_2) {
      if (len_s < len_1) {
        # 2 is longest
        if (gmix_term_s != gmix_term_2[1:len_s]) {
          stop("Error: Second model and shared model have different geometric mixture term values.")
        }
        if (gmix_term_1 != gmix_term_2[1:len_1]) {
          stop("Error: Second model and first model have different geometric mixture term values.")
        }
        gmix_term <- gmix_term_2
      }
    }
  }
  #
  model <- poismodel(pyr, "events", strata, term_n, tform, names, modelform, gmix_term, gmix_theta, a_n, keep_constant, df)
  list(
    "model" = model, "data" = df
  )
}


#' Interprets a Colossus formula and makes necessary changes to data
#'
#' \code{get_form} uses a formula and data.table, to fully describe the model
#' for a Colossus regression function.
#'
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the updated data
#' @family Formula Interpretation
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
#'   "d" = c(0, 0, 0, 1, 1, 1, 1),
#'   "e" = c(0, 0, 1, 0, 0, 0, 1)
#' )
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
#'   loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' model <- get_form(formula, df)
get_form <- function(formula, df) {
  if (length(lapply(strsplit(format(formula), ""), function(x) which(x == "~"))[[1]]) != 1) {
    stop("Error: The formula contained multiple '~', invalid formula")
  }
  surv_obj <- format(formula[[2]])
  model_obj <- paste(format(formula[[3]]), collapse = " ")
  surv_obj <- gsub(" ", "", surv_obj)
  model_obj <- gsub(" ", "", model_obj)
  res <- get_form_list(surv_obj, model_obj, df)
  model <- res$model
  df <- res$data
  #
  surv_model_type <- model$surv_model_type
  tstart <- model$tstart
  tend <- model$tend
  pyr <- model$pyr
  trials <- model$trials
  event <- model$event
  strata <- model$strata
  weight <- model$weight
  #
  term_n <- model$term_n
  tform <- model$tform
  names <- model$names
  modelform <- model$modelform
  gmix_term <- model$gmix_term
  gmix_theta <- model$gmix_theta
  null <- model$null
  #
  expres_calls <- model$expres_calls
  #
  if (grepl("cox", surv_model_type)) {
    model <- coxmodel(tstart, tend, event, strata, weight, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df, expres_calls)
  } else if (grepl("finegray", surv_model_type)) {
    model <- coxmodel(tstart, tend, event, strata, weight, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df, expres_calls)
  } else if (grepl("pois", surv_model_type)) {
    model <- poismodel(pyr, event, strata, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df, expres_calls)
    if (all(strata != "NONE")) {
      Check_Strata_Model(term_n, tform, modelform, gmix_term, gmix_theta)
    } # verifies that a stratified model can be used
  } else if ((grepl("casecon", surv_model_type)) || (grepl("case_con", surv_model_type))) {
    model <- caseconmodel(tstart, tend, event, strata, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df, expres_calls)
  } else if ((grepl("logit", surv_model_type)) || (grepl("logistic", surv_model_type))) {
    model <- logitmodel(trials, event, strata, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df, expres_calls)
  } else {
    stop("Error: Bad survival model type passed")
  }
  list(
    "model" = model, "data" = df
  )
}

#' Interprets a Colossus formula and makes necessary changes to data, returns list not class
#'
#' \code{get_form_list} uses a formula and data.table, to fully describe the model
#' for a Colossus regression function and returns a list
#'
#' @param surv_obj output from get_form_surv, list of survival values
#' @param model_obj output from get_form_risk, list of risk factor model values
#' @inheritParams R_template
#'
#' @noRd
#' @return returns the list of model values
#' @family Formula Interpretation
get_form_list <- function(surv_obj, model_obj, df) {
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        df <- setDT(df)
      },
      error = function(e) {
        df <- data.table(df)
      }
    )
  }
  surv_obj <- gsub(" ", "", surv_obj)
  model_obj <- gsub(" ", "", model_obj)
  #
  surv_list <- get_form_surv(surv_obj, df)
  df <- surv_list$data
  surv_vals <- surv_list$surv
  surv_model_type <- surv_vals$surv_model_type
  tstart <- surv_vals$tstart
  tend <- surv_vals$tend
  pyr <- surv_vals$pyr
  event <- surv_vals$event
  strata <- surv_vals$strata
  weight <- surv_vals$weight
  trials <- surv_vals$trials
  #
  model_vals <- get_form_risk(model_obj, df)
  model <- model_vals$model
  term_n <- model$term_n
  tform <- model$tform
  names <- model$names
  modelform <- model$modelform
  gmix_term <- model$gmix_term
  gmix_theta <- model$gmix_theta
  null <- model$null
  df <- model_vals$data
  expres_calls <- model$expres_calls
  #
  list(
    "model" = list(
      surv_model_type = surv_model_type,
      tstart = tstart,
      tend = tend,
      pyr = pyr,
      trials = trials,
      event = event,
      strata = strata,
      weight = weight,
      term_n = term_n,
      tform = tform,
      names = names,
      modelform = modelform,
      gmix_term = gmix_term,
      gmix_theta = gmix_theta,
      null = null,
      expres_calls = expres_calls
    ),
    "data" = df
  )
}

#' Interprets the survival side of a formula
#'
#' \code{get_form_surv} interprets the LHS of a formula
#' @param surv_obj output from get_form_surv, list of survival values
#' @inheritParams R_template
#' @noRd
#' @return returns the left hand side of the formula reading
#' @family Formula Interpretation
get_form_surv <- function(surv_obj, df) {
  surv_obj <- gsub(" ", "", surv_obj)
  surv_obj <- gsub("\"", "", surv_obj) # remove literal strings if needed
  surv_model_type <- "NONE"
  tstart <- "NONE"
  tend <- "NONE"
  pyr <- "NONE"
  event <- "NONE"
  strata <- "NONE"
  weight <- "NONE"
  trials <- "NONE"
  #
  # split the survival model type
  second_split <- lapply(strsplit(surv_obj, ""), function(x) which(x == "("))[[1]]
  if (length(second_split) == 0) {
    stop("Error: The left hand side did not contain (")
  }
  surv_type <- tolower(substr(surv_obj, 1, second_split - 1))
  surv_paras <- substr(surv_obj, second_split + 1, nchar(surv_obj) - 1)
  surv_paras <- nested_split(surv_paras)
  # assign survival model values
  surv_model_type <- surv_type
  surv_para_list <- list()
  for (i in seq_along(surv_paras)) {
    para_cur <- surv_paras[i]
    para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]]
    if (length(para_break) == 0) {
      # no name, just add to list
      surv_para_list[[i]] <- para_cur
    } else {
      item_name <- substr(para_cur, 1, para_break - 1)
      item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
      surv_para_list[[item_name]] <- item_value
    }
  }
  if (surv_type %in% c("cox", "coxph")) {
    res <- do.call(ColossusCoxSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
  } else if (surv_type %in% c("cox_strata", "coxph_strata")) {
    res <- do.call(ColossusCoxStrataSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    strata <- res$strata
  } else if (surv_type %in% c("finegray", "fine_gray", "fg")) {
    res <- do.call(ColossusFineGraySurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    weight <- res$weight
  } else if (surv_type %in% c("finegray_strata", "fine_gray_strata", "fg_strata")) {
    res <- do.call(ColossusFineGrayStrataSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    strata <- res$strata
    weight <- res$weight
  } else if (surv_type %in% c("poisson", "pois")) {
    res <- do.call(ColossusPoisSurv, surv_para_list)
    pyr <- res$pyr
    event <- res$event
    if (res$strata != "NULL") {
      stop("Error: Too many columns passed to non-stratified Poisson model")
    }
  } else if (surv_type %in% c("poisson_strata", "pois_strata")) {
    res <- do.call(ColossusPoisSurv, surv_para_list)
    pyr <- res$pyr
    event <- res$event
    strata <- res$strata
    if (length(res$strata) == 1) {
      if (res$strata == "NULL") {
        stop("Error: Too few columns passed to non-stratified Poisson model")
      }
    }
  } else if (surv_type %in% c("casecon", "casecontrol", "case_control")) {
    res <- do.call(ColossusCaseConSurv, surv_para_list)
    event <- res$event
  } else if (surv_type %in% c("casecon_time", "casecontrol_time", "case_control_time")) {
    res <- do.call(ColossusCaseConTimeSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
  } else if (surv_type %in% c("casecon_strata", "casecontrol_strata", "case_control_strata")) {
    res <- do.call(ColossusCaseConStrataSurv, surv_para_list)
    event <- res$event
    strata <- res$strata
  } else if (surv_type %in% c("casecon_strata_time", "casecontrol_strata_time", "case_control_strata_time", "casecon_time_strata", "casecontrol_time_strata", "case_control_time_strata")) {
    res <- do.call(ColossusCaseConTimeStrataSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    strata <- res$strata
  } else if (surv_type %in% c("logit", "logistic")) {
    res <- do.call(ColossusLogitSurv, surv_para_list)
    trials <- res$trials
    event <- res$event
    if ("CONST" == trials) {
      if ("CONST" %in% names(df)) {
        # fine
      } else {
        df$CONST <- 1
      }
    }
  } else {
    stop("Error: Invalid survival model type")
  }
  #
  if (tstart == "%trunc%") {
    if (tend == "%trunc%") {
      stop("Error: Both endpoints are truncated, not acceptable")
    }
    tmin <- min(df[, get(tend)]) - 1
    if (!("right_trunc" %in% names(df))) {
      df[, ":="(right_trunc = tmin)]
    }
    tstart <- "right_trunc"
  } else if (tend == "%trunc%") {
    tmax <- max(df[, get(tstart)]) + 1
    if (!("left_trunc" %in% names(df))) {
      df[, ":="(left_trunc = tmax)]
    }
    tend <- "left_trunc"
  }
  #
  list(
    "surv" = list(
      surv_model_type = surv_model_type,
      tstart = tstart,
      tend = tend,
      trials = trials,
      pyr = pyr,
      event = event,
      strata = strata,
      weight = weight
    ),
    "data" = df
  )
}

#' Interprets the risk factor side of a formula
#'
#' \code{get_form_risk} interprets the RHS of a formula
#' @param model_obj output from get_form_risk, list of risk factor model values
#' @inheritParams R_template
#' @noRd
#' @return returns the right hand side of a formula reading
#' @family Formula Interpretation
get_form_risk <- function(model_obj, df) {
  #
  model_obj <- gsub(" ", "", model_obj)
  term_n <- c()
  tform <- c()
  names <- c()
  gmix_theta <- 0
  gmix_term <- c()
  modelform <- "NONE"
  null <- FALSE
  expres_calls <- list() # storing any variables transformed, ie factor, ns, bs, etc. to be recalled when the data is reloaded
  #
  right_model_terms <- strsplit(model_obj, "\\+")[[1]]
  for (term_i in seq_along(right_model_terms)) {
    # seperate the term type or model-formula from parameters
    third_split <- lapply(strsplit(right_model_terms[term_i], ""), function(x) which(x == "("))[[1]]
    if (length(third_split) == 0) {
      stop(paste('Error: right hand side element "', right_model_terms[term_i], '" did not contain (', sep = ""))
    }
    model_type <- tolower(substr(right_model_terms[term_i], 1, third_split - 1))
    model_type <- gsub("_", "-", model_type)
    tform_acceptable <- c(
      "plin", "lin", "loglin", "loglin-dose", "lin-dose",
      "lin-quad-dose", "lin-exp-dose", "plinear", "product-linear", "linear",
      "loglinear", "log-linear", "loglinear-dose", "log-linear-dose", "linear-dose", "linear-piecewise",
      "quadratic", "quad", "quad-dose", "quadratic-dose",
      "step-dose", "step-piecewise",
      "linear-quadratic-dose", "linear-quadratic-piecewise",
      "linear-exponential-dose", "linear-exponential-piecewise"
    )
    modelform_acceptable <- c(
      "multiplicative", "multiplicative-excess", "additive", "product-additive",
      "product-additive-excess", "a", "pa", "pae", "m", "me",
      "gmix", "geometric-mixture", "gmix-r", "relative-geometric-mixture",
      "gmix-e", "excess-geometric-mixture"
    )
    if (model_type == "null") {
      null <- TRUE
    } else if (model_type %in% tform_acceptable) {
      model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
      model_paras <- nested_split(model_paras)
      last_entry <- model_paras[length(model_paras)]
      if (is.na(suppressWarnings(as.integer(last_entry)))) {
        # the last element isn't an integer
        term_num <- 0
      } else {
        model_paras <- model_paras[1:(length(model_paras) - 1)]
        term_num <- as.integer(last_entry)
      }
      for (subterm_i in seq_along(model_paras)) {
        # add element
        # check if the element is a function
        if (grepl("\\(", model_paras[subterm_i])) {
          # Some function is being used
          if (substr(model_paras[subterm_i], 1, 7) == "factor(") {
            # baseline is set by using factor(column;baseline=level)
            factor_args <- substr(model_paras[subterm_i], 8, nchar(model_paras[subterm_i]) - 1)
            factor_args <- nested_split(factor_args)
            ##
            factor_arg_list <- list()
            for (i in seq_along(factor_args)) {
              para_cur <- factor_args[i]
              para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]]
              if (length(para_break) == 0) {
                # no name, just add to list
                factor_arg_list[[i]] <- para_cur
              } else {
                item_name <- substr(para_cur, 1, para_break - 1)
                item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
                factor_arg_list[[item_name]] <- parse_literal_string(item_value)
              }
            }
            repeat_list <- c(list("_exp_type" = "factor"), copy(factor_arg_list)) # The arguements needed to repeat the processing
            # Either the item is named x, or is it the first
            if ("x" %in% names(factor_arg_list)) {
              factor_col <- factor_arg_list$x
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            } else {
              factor_col <- factor_arg_list[[1]]
              names(factor_arg_list)[[1]] <- "x"
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            }
            if (!(factor_col %in% names(df))) {
              stop(paste("Error: Column: ", factor_col, " not in data", sep = ""))
            }
            xtemp <- do.call(factor, factor_arg_list)
            if (!("levels" %in% names(repeat_list))) { # using the levels will recreate the same factoring
              repeat_list[["levels"]] <- levels(xtemp)
            }
            df[[factor_col]] <- xtemp
            val <- factorize(df, factor_col)
            df <- val$df
            col_name <- val$cols
            level_ref <- paste(factor_col, levels(xtemp)[1], sep = "_")
            col_name <- col_name[col_name != level_ref]
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
          } else if (substr(model_paras[subterm_i], 1, 2) == "I(") {
            factor_args <- substr(model_paras[subterm_i], 3, nchar(model_paras[subterm_i]) - 1)
            repeat_list <- c(list("_exp_type" = "power"), list(factor_args))
            vals <- strsplit(factor_args, "\\^")[[1]]
            if (length(vals) != 2) {
              stop("Error: I() currently only available for I(var^n)")
            }
            col <- vals[1]
            raised <- vals[2]
            if (!(col %in% names(df))) {
              stop(paste("Error: Column: ", col, " not in data", sep = ""))
            }
            options(warn = -1)
            if (all(vapply(raised, function(x) grepl("^[\\-]{0,1}[0-9]*\\.{0,1}[0-9]*$", x), logical(1))) || all(vapply(raised, function(x) grepl("^[\\-]{0,1}[0-9]+e[\\-]{0,1}[0-9]+$", x), logical(1)))) {
              options(warn = 0) # checks for an integer, decimal, decimal places or scientific notation
              raised <- as.numeric(raised)
            } else {
              stop("Error: Column was not raised to a numeric power")
            }
            col_name <- factor_args
            df[[col_name]] <- df[[col]]^raised
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
          } else if (substr(model_paras[subterm_i], 1, 3) == "ns(") {
            # natural cubic spline
            factor_args <- substr(model_paras[subterm_i], 4, nchar(model_paras[subterm_i]) - 1)
            factor_args <- nested_split(factor_args)
            ##
            factor_arg_list <- list()
            for (i in seq_along(factor_args)) {
              para_cur <- factor_args[i]
              para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]]
              if (length(para_break) == 0) {
                # no name, just add to list
                factor_arg_list[[i]] <- para_cur
              } else {
                item_name <- substr(para_cur, 1, para_break - 1)
                item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
                factor_arg_list[[item_name]] <- parse_literal_string(item_value)
              }
            }
            repeat_list <- c(list("_exp_type" = "ns"), copy(factor_arg_list)) # The arguements needed to repeat the processing
            # Either the item is named x, or is it the first
            if ("x" %in% names(factor_arg_list)) {
              factor_col <- factor_arg_list$x
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            } else {
              factor_col <- factor_arg_list[[1]]
              names(factor_arg_list)[[1]] <- "x"
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            }
            if (system.file(package = "splines") == "") {
              stop("Error: Attempted to use ns(), but splines not detected on system.")
            }
            if (!(factor_col %in% names(df))) {
              stop(paste("Error: Column: ", factor_col, " not in data", sep = ""))
            }
            xtemp <- do.call(splines::ns, factor_arg_list)
            #
            ns_att <- attributes(xtemp)
            # Applying the same knots, boundary, and intercept gives the same transformation
            if (!("knots" %in% names(repeat_list))) {
              repeat_list[["knots"]] <- ns_att$knots
            }
            if (!("Boundary.knots" %in% names(repeat_list))) {
              repeat_list[["Boundary.knots"]] <- ns_att$Boundary.knots
            }
            if (!("intercept" %in% names(repeat_list))) {
              repeat_list[["intercept"]] <- ns_att$intercept
            }
            #
            col_name <- c()
            for (i in seq_len(ncol(xtemp))) {
              x_col <- paste(factor_col, "_ns", i, sep = "")
              if (!(x_col %in% names(df))) {
                df[[x_col]] <- xtemp[, i]
              }
              col_name <- c(col_name, x_col)
            }
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
            ##
          } else if (substr(model_paras[subterm_i], 1, 3) == "bs(") {
            # b-spline for polynomial spline
            factor_args <- substr(model_paras[subterm_i], 4, nchar(model_paras[subterm_i]) - 1)
            factor_args <- nested_split(factor_args)
            ##
            factor_arg_list <- list()
            for (i in seq_along(factor_args)) {
              para_cur <- factor_args[i]
              para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]]
              if (length(para_break) == 0) {
                # no name, just add to list
                factor_arg_list[[i]] <- para_cur
              } else {
                item_name <- substr(para_cur, 1, para_break - 1)
                item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
                factor_arg_list[[item_name]] <- parse_literal_string(item_value)
              }
            }
            repeat_list <- c(list("_exp_type" = "bs"), copy(factor_arg_list)) # The arguements needed to repeat the processing
            # Either the item is named x, or is it the first
            if ("x" %in% names(factor_arg_list)) {
              factor_col <- factor_arg_list$x
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            } else {
              factor_col <- factor_arg_list[[1]]
              names(factor_arg_list)[[1]] <- "x"
              factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
            }
            if (system.file(package = "splines") == "") {
              stop("Error: Attempted to use bs(), but splines not detected on system.")
            }
            if (!(factor_col %in% names(df))) {
              stop(paste("Error: Column: ", factor_col, " not in data", sep = ""))
            }
            xtemp <- do.call(splines::bs, factor_arg_list)
            #
            bs_att <- attributes(xtemp)
            # Applying the same degree, knots, boundary, and intercept gives the same transformation
            if (!("degree" %in% names(repeat_list))) {
              repeat_list[["degree"]] <- bs_att$degree
            }
            if (!("knots" %in% names(repeat_list))) {
              repeat_list[["knots"]] <- bs_att$knots
            }
            if (!("Boundary.knots" %in% names(repeat_list))) {
              repeat_list[["Boundary.knots"]] <- bs_att$Boundary.knots
            }
            if (!("intercept" %in% names(repeat_list))) {
              repeat_list[["intercept"]] <- bs_att$intercept
            }
            col_name <- c()
            for (i in seq_len(ncol(xtemp))) {
              x_col <- paste(factor_col, "_bs", i, sep = "")
              if (!(x_col %in% names(df))) {
                df[[x_col]] <- xtemp[, i]
              }
              col_name <- c(col_name, x_col)
            }
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
            ##
          } else {
            # ----------------------------------------------------------------------------------- #
            # generic function call
            para_split <- lapply(strsplit(model_paras[subterm_i], ""), function(x) which(x == "("))[[1]]
            exp_name <- substr(model_paras[subterm_i], 1, para_split - 1)
            factor_args <- substr(model_paras[subterm_i], para_split + 1, nchar(model_paras[subterm_i]) - 1)
            factor_args <- nested_split(factor_args)
            ##
            factor_arg_list <- list()
            for (i in seq_along(factor_args)) {
              para_cur <- factor_args[i]
              para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]]
              if (length(para_break) == 0) {
                # no name, just add to list
                factor_arg_list[[i]] <- para_cur
              } else {
                item_name <- substr(para_cur, 1, para_break - 1)
                item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
                factor_arg_list[[item_name]] <- parse_literal_string(item_value)
              }
            }
            repeat_list <- c(list("_exp_type" = exp_name), copy(factor_arg_list)) # The arguements needed to repeat the processing
            # Either the item is named x, or is it the first
            factor_col <- factor_arg_list[[1]]
            factor_vals <- unlist(factor_arg_list, use.names = FALSE)
            for (i in seq_along(factor_vals)) {
              if (factor_vals[[i]] %in% names(df)) {
                factor_arg_list[[i]] <- copy(df[[factor_vals[[i]]]])
              }
            }
            if (!(factor_col %in% names(df))) {
              stop(paste("Error: Column: ", factor_col, " not in data", sep = ""))
            }
            xtemp <- do.call(exp_name, factor_arg_list)
            #
            col_name <- c()
            if (is.atomic(xtemp)) {
              x_col <- paste(factor_col, "_", exp_name, sep = "")
              if (!(x_col %in% names(df))) {
                df[[x_col]] <- xtemp
              }
              col_name <- c(col_name, x_col)
            } else {
              for (i in seq_len(ncol(xtemp))) {
                x_col <- paste(factor_col, "_", exp_name, i, sep = "")
                if (!(x_col %in% names(df))) {
                  df[[x_col]] <- xtemp[, i]
                }
                col_name <- c(col_name, x_col)
              }
            }
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
          }
        } else if (grepl("\\*", model_paras[subterm_i])) {
          # interaction element
          repeat_list <- c(list("_exp_type" = "interaction"), model_paras[subterm_i]) # The arguements needed to repeat the processing
          #
          # split into the columns
          cols <- strsplit(model_paras[subterm_i], "\\*")[[1]]
          for (col in cols) {
            if (!(col %in% names(df))) {
              # good, it is in there
              stop(paste("Error: Interaction column missing: ", col, sep = ""))
            }
          }

          recur_interact <- function(x, y) {
            if (length(y) == 1) {
              return(paste(x, y[[1]], sep = ":"))
            } else {
              for (i in y[[1]]) {
                y0 <- paste(x, i, sep = ":")
                return(recur_interact(y0, y[2:length(y)]))
              }
            }
          }

          interact_tables <- do.call(c, lapply(seq_along(cols), combn, x = cols, simplify = FALSE))
          col_name <- c()
          for (i_table in interact_tables) {
            vals <- unlist(i_table, use.names = FALSE)
            if (length(vals) == 1) {
              if (is.factor(df[[vals]])) {
                # factor
                val <- factorize(df, vals)
                df <- val$df
                fac_names <- val$cols
                level_ref <- paste(vals, levels(df[[vals]])[1], sep = "_")
                col_name <- c(col_name, fac_names[fac_names != level_ref])
              } else {
                # not factor
                col_name <- c(col_name, vals)
              }
            } else {
              # there are multiple to combine
              # get the levels for each element
              element_levels <- list()
              for (term_i in seq_along(vals)) {
                factor_col <- vals[term_i]
                if (is.factor(df[[factor_col]])) {
                  i_levels <- paste(factor_col, levels(df[[factor_col]]), sep = "_")
                  level_ref <- paste(factor_col, levels(df[[factor_col]])[1], sep = "_")
                  i_levels <- i_levels[i_levels != level_ref]
                  element_levels[[term_i]] <- i_levels
                } else {
                  element_levels[[term_i]] <- factor_col
                }
              }
              combs <- c()
              for (i in element_levels[[1]]) {
                y0 <- i
                combs <- recur_interact(y0, element_levels[2:length(element_levels)])
                for (comb in combs) {
                  #
                  entries <- strsplit(comb, ":")[[1]]
                  df[[comb]] <- 1
                  for (j in entries) {
                    df[[comb]] <- df[[comb]] * df[[j]]
                  }
                  col_name <- c(col_name, comb)
                }
              }
            }
          }
          expres_calls[[length(expres_calls) + 1]] <- repeat_list
        } else {
          # check if the column is actually a factor
          element_col <- model_paras[subterm_i]
          if (is.null(levels(df[[element_col]]))) {
            col_name <- c(element_col)
            if ("CONST" == model_paras[subterm_i]) {
              if ("CONST" %in% names(df)) {
                # fine
              } else {
                df$CONST <- 1
              }
            }
          } else {
            # it is a factor
            factor_arg_list <- list()
            factor_arg_list[[1]] <- element_col
            repeat_list <- c(list("_exp_type" = "factor"), copy(factor_arg_list))
            repeat_list[["levels"]] <- levels(df[[element_col]])
            expres_calls[[length(expres_calls) + 1]] <- repeat_list
            val <- factorize(df, element_col)
            df <- val$df
            col_name <- val$cols
            level_ref <- paste(element_col, levels(df[[element_col]])[1], sep = "_")
            col_name <- col_name[col_name != level_ref]
          }
        }
        if (length(col_name) == 0) {
          stop("Error: Subterm missing element, was a single-valued factor used?")
        }
        col_name <- unique(col_name)
        # convert subterm formula
        model_terms <- c(model_type)
        if (model_type %in% c("plin", "plinear", "product-linear")) {
          model_terms <- c("plin")
        } else if (model_type %in% c("lin", "linear")) {
          model_terms <- c("lin")
        } else if (model_type %in% c("loglin", "loglinear", "log-linear")) {
          model_terms <- c("loglin")
        } else if (model_type %in% c("loglin-dose", "loglinear-dose", "log-linear-dose")) {
          model_terms <- c("loglin_slope", "loglin_top")
        } else if (model_type %in% c("lin-dose", "linear-dose", "linear-piecewise")) {
          model_terms <- c("lin_slope", "lin_int")
        } else if (model_type %in% c("quadratic", "quad", "quad-dose", "quadratic-dose")) {
          model_terms <- c("quad_slope")
        } else if (model_type %in% c("step-dose", "step-piecewise")) {
          model_terms <- c("step_slope", "step_int")
        } else if (model_type %in% c("lin-quad-dose", "linear-quadratic-dose", "linear-quadratic-piecewise")) {
          model_terms <- c("lin_quad_slope", "lin_quad_int")
        } else if (model_type %in% c("lin-exp-dose", "linear-exponential-dose", "linear-exponential-piecewise")) {
          model_terms <- c("lin_exp_slope", "lin_exp_int", "lin_exp_exp_slope")
        } else {
          stop(paste("Error: Unknown subterm type used, ", model_type, sep = ""))
        }
        for (col in col_name) {
          for (model_term in model_terms) {
            if (grepl("_int", model_term)) {
              # we want to create a second column, for the intercept, to be normalized differently
              new_col <- paste(col, ":intercept", sep = "")
              if (!(new_col %in% names(df))) {
                df[, new_col] <- df[, col, with = FALSE]
              }
              names <- c(names, new_col)
            } else {
              names <- c(names, col)
            }
            tform <- c(tform, model_term)
            term_n <- c(term_n, term_num)
          }
        }
      }
    } else if (model_type %in% modelform_acceptable) {
      if (model_type %in% c("m", "multiplicative")) {
        model_type <- "M"
      } else if (model_type %in% c("me", "multiplicative-excess")) {
        model_type <- "ME"
      } else if (model_type %in% c("a", "additive")) {
        model_type <- "A"
      } else if (model_type %in% c("pa", "product-additive")) {
        model_type <- "PA"
      } else if (model_type %in% c("pae", "product-additive-excess")) {
        model_type <- "PAE"
      } else if (model_type %in% c("gm", "gmix", "geometric-mixture")) {
        model_type <- "GMIX"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        model_paras <- gsub("\"", "", model_paras) # remove literal strings if needed
        #
        if (is.na(model_paras) || model_paras == "") {
          # Nothing, just give the defaults
          gmix_theta <- 0.5
          gmix_term <- c()
        } else {
          # There is something, but what?
          if (grepl(",", model_paras)) {
            # Multiple items
            model_paras <- tolower(strsplit(model_paras, ",")[[1]])
            # we need to start by checking if the first thing is a number
            if (all(vapply(model_paras[1], function(x) grepl("^[\\-]{0,1}[0-9]*\\.{0,1}[0-9]*$", x), logical(1))) || all(vapply(model_paras[1], function(x) grepl("^[\\-]{0,1}[0-9]+e[\\-]{0,1}[0-9]+$", x), logical(1)))) {
              # Good, the first item is a number
              gmix_theta <- as.numeric(model_paras[1])
              para_else <- model_paras[2:length(model_paras)]
              # Check if they are all valid
              if (all(vapply(para_else, function(x) grepl(x, "er"), logical(1)))) {
                gmix_term <- ifelse(para_else == "e", 1, 0)
              } else {
                # Remaining entry was wrong
                stop("Error: Gmix term had an invalid option after the theta value. Please only use 'e/r'")
              }
            } else {
              # wasn't a number
              gmix_theta <- 0.5
              para_else <- model_paras[1:length(model_paras)]
              # Check if they are all valid
              if (all(vapply(para_else, function(x) grepl(x, "er"), logical(1)))) {
                gmix_term <- ifelse(para_else == "e", 1, 0)
              } else {
                # Remaining entry was wrong
                stop("Error: Gmix term had an invalid option. Please only use a number or 'e/r'")
              }
            }
          } else {
            # Single Item
            if (all(vapply(model_paras, function(x) grepl("^[\\-]{0,1}[0-9]*\\.{0,1}[0-9]*$", x), logical(1))) || all(vapply(model_paras, function(x) grepl("^[\\-]{0,1}[0-9]+e[\\-]{0,1}[0-9]+$", x), logical(1)))) {
              # It is a number, make it the gmix_theta
              gmix_theta <- as.numeric(model_paras)
              gmix_term <- c()
            } else {
              # it is not a number, is it a 'e' or 'r'?
              if (grepl(tolower(model_paras), "er")) {
                # it is an option!
                gmix_theta <- 0.5
                gmix_term <- ifelse(tolower(model_paras) == "e", 1, 0)
              } else {
                # Doesn't match anything
                stop("Error: Gmix term had an invalid option. Please only use a number or 'e/r'")
              }
            }
          }
        }
        #
      } else if (model_type %in% c("gmix-r", "relative-geometric-mixture")) {
        model_type <- "GMIX-R"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        if (is.na(model_paras) || model_paras == "") {
          # Nothing, just give the defaults
          gmix_theta <- 0.5
        } else {
          model_paras <- tolower(strsplit(model_paras, ",")[[1]])
          gmix_theta <- as.numeric(model_paras[1])
        }
      } else if (model_type %in% c("gmix-e", "excess-geometric-mixture")) {
        model_type <- "GMIX-E"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        if (is.na(model_paras) || model_paras == "") {
          # Nothing, just give the defaults
          gmix_theta <- 0.5
        } else {
          model_paras <- tolower(strsplit(model_paras, ",")[[1]])
          gmix_theta <- as.numeric(model_paras[1])
        }
      }
      if (modelform == "NONE") {
        modelform <- model_type
      } else {
        stop("Error: modelform defined twice")
      }
    } else {
      stop(paste("Error: Unknown option encountered, ", model_type, sep = ""))
    }
  }
  if (!null) {
    term_tot <- max(term_n) + 1
    if (modelform == "NONE") {
      modelform <- "ME"
    } else if (modelform == "GMIX-R") {
      modelform <- "GMIX"
      gmix_term <- rep(0, term_tot)
    } else if (modelform == "GMIX-E") {
      modelform <- "GMIX"
      gmix_term <- rep(1, term_tot)
    }
    if (length(gmix_term) < term_tot) {
      gmix_term <- c(gmix_term, rep(1.0, term_tot - length(gmix_term)))
    } else if (length(gmix_term) > term_tot) {
      stop("Error: The gmix option was used with more values than terms")
    }
  }
  #
  if (length(tform) > 1) {
    og_index <- seq_along(tform)
    all_remove <- rep(0, length(tform))
    for (u_tform in unique(tform)) {
      for (u_term_n in unique(term_n)) {
        new_name <- names[(tform == u_tform) & (term_n == u_term_n)]
        new_index <- og_index[(tform == u_tform) & (term_n == u_term_n)]
        if (length(new_name) > 1) {
          df_temp <- df[, new_name, with = FALSE]
          removed <- c(0)
          for (i in 2:length(names(df_temp))) {
            new_rank <- qr(df_temp[, 1:i])$rank
            if (new_rank == i - sum(removed)) {
              # all good
              removed <- c(removed, 0)
            } else {
              # the column was linearlly dependent
              removed <- c(removed, 1)
              all_remove[new_index[i]] <- 1
            }
          }
        }
      }
    }
    term_n <- term_n[all_remove == 0]
    tform <- tform[all_remove == 0]
    names <- names[all_remove == 0]
  }
  #
  model <- list(
    "term_n" = term_n,
    "tform" = tform,
    "names" = names,
    "modelform" = modelform,
    "gmix_term" = gmix_term,
    "gmix_theta" = gmix_theta,
    "null" = null,
    "expres_calls" = expres_calls
  )
  list(
    "model" = model, "data" = df
  )
}

#' Applies a list of function calls from a model evaluation
#'
#' \code{ColossusExpressionCall} Uses a list of function call info and a data.table, and applies the function calls
#'
#' @param calls List from 'model$expres' calls formatted as the type of function call, and then the arguments
#' @param df dataset, to be modified by the function calls
#' @noRd
#' @return returns the transformed data
#' @family Formula Interpretation
ColossusExpressionCall <- function(calls, df) {
  for (call in calls) {
    if (call[["_exp_type"]] == "factor") {
      factor_arg_list <- call[names(call) != "_exp_type"]
      if ("x" %in% names(factor_arg_list)) {
        factor_col <- factor_arg_list$x
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      } else {
        factor_col <- factor_arg_list[[1]]
        names(factor_arg_list)[[1]] <- "x"
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      }
      xtemp <- do.call(factor, factor_arg_list)
      df[[factor_col]] <- xtemp
      val <- factorize(df, factor_col)
      df <- val$df
    } else if (call[["_exp_type"]] == "power") {
      factor_args <- call[names(call) != "_exp_type"][[1]]
      vals <- strsplit(factor_args, "\\^")[[1]]
      col <- vals[1]
      raised <- vals[2]
      raised <- as.numeric(raised)
      df[[factor_args]] <- df[[col]]^raised
    } else if (call[["_exp_type"]] == "ns") {
      # natural cubic spline
      factor_arg_list <- call[names(call) != "_exp_type"]
      if ("x" %in% names(factor_arg_list)) {
        factor_col <- factor_arg_list$x
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      } else {
        factor_col <- factor_arg_list[[1]]
        names(factor_arg_list)[[1]] <- "x"
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      }
      if (system.file(package = "splines") == "") {
        stop("Error: Attempted to use ns(), but splines not detected on system.")
      }
      xtemp <- do.call(splines::ns, factor_arg_list)
      #
      for (i in seq_len(ncol(xtemp))) {
        x_col <- paste(factor_col, "_ns", i, sep = "")
        if (!(x_col %in% names(df))) {
          df[[x_col]] <- xtemp[, i]
        }
      }
      ##
    } else if (call[["_exp_type"]] == "bs") {
      # b-spline for polynomial spline
      factor_arg_list <- call[names(call) != "_exp_type"]
      if ("x" %in% names(factor_arg_list)) {
        factor_col <- factor_arg_list$x
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      } else {
        factor_col <- factor_arg_list[[1]]
        names(factor_arg_list)[[1]] <- "x"
        factor_arg_list[["x"]] <- copy(df[[factor_arg_list$x]])
      }
      if (system.file(package = "splines") == "") {
        stop("Error: Attempted to use bs(), but splines not detected on system.")
      }
      xtemp <- do.call(splines::bs, factor_arg_list)
      #
      for (i in seq_len(ncol(xtemp))) {
        x_col <- paste(factor_col, "_bs", i, sep = "")
        if (!(x_col %in% names(df))) {
          df[[x_col]] <- xtemp[, i]
        }
      }
      ##
    } else if (call["_exp_type"] == "interaction") {
      # split into the columns
      cols <- strsplit(call[[2]], "\\*")[[1]]
      recur_interact <- function(x, y) {
        if (length(y) == 1) {
          return(paste(x, y[[1]], sep = ":"))
        } else {
          for (i in y[[1]]) {
            y0 <- paste(x, i, sep = ":")
            return(recur_interact(y0, y[2:length(y)]))
          }
        }
      }
      interact_tables <- do.call(c, lapply(seq_along(cols), combn, x = cols, simplify = FALSE))
      for (i_table in interact_tables) {
        vals <- unlist(i_table, use.names = FALSE)
        if (length(vals) == 1) {
          if (is.factor(df[[vals]])) {
            # factor
            val <- factorize(df, vals)
            df <- val$df
          }
        } else {
          # there are multiple to combine
          # get the levels for each element
          element_levels <- list()
          for (term_i in seq_along(vals)) {
            factor_col <- vals[term_i]
            if (is.factor(df[[factor_col]])) {
              i_levels <- paste(factor_col, levels(df[[factor_col]]), sep = "_")
              level_ref <- paste(factor_col, levels(df[[factor_col]])[1], sep = "_")
              i_levels <- i_levels[i_levels != level_ref]
              element_levels[[term_i]] <- i_levels
            } else {
              element_levels[[term_i]] <- factor_col
            }
          }
          combs <- c()
          for (i in element_levels[[1]]) {
            y0 <- i
            combs <- recur_interact(y0, element_levels[2:length(element_levels)])
            for (comb in combs) {
              #
              entries <- strsplit(comb, ":")[[1]]
              df[[comb]] <- 1
              for (j in entries) {
                df[[comb]] <- df[[comb]] * df[[j]]
              }
            }
          }
        }
      }
    } else {
      # --------------------------------------------------------------------- #
      exp_name <- call[["_exp_type"]]
      factor_arg_list <- call[names(call) != "_exp_type"]
      # Either the item is named x, or is it the first
      factor_col <- factor_arg_list[[1]]
      factor_vals <- unlist(factor_arg_list, use.names = FALSE)
      for (i in seq_along(factor_vals)) {
        if (factor_vals[[i]] %in% names(df)) {
          factor_arg_list[[i]] <- copy(df[[factor_vals[[i]]]])
        }
      }

      xtemp <- do.call(exp_name, factor_arg_list)
      #
      if (is.atomic(xtemp)) {
        x_col <- paste(factor_col, "_", exp_name, sep = "")
        if (!(x_col %in% names(df))) {
          df[[x_col]] <- xtemp
        }
      } else {
        for (i in seq_len(ncol(xtemp))) {
          x_col <- paste(factor_col, "_", exp_name, i, sep = "")
          if (!(x_col %in% names(df))) {
            df[[x_col]] <- xtemp[, i]
          }
        }
      }
    }
  }
  df
}

#' Interprets basic cox survival formula RHS
#'
#' \code{ColossusCoxSurv} assigns and interprets interval columns for cox model.
#' This functions is called using the arguments for Cox in the right-hand side of
#' the formula. Uses an interval start time, end time, and event status. These are
#' expected to be in order or named: tstart, tend, and event.
#' The Fine-Gray and Stratified versions use strata and weight named options
#' or the last two entries.
#'
#' @param ... entries for a cox survival object, tstart, tend, and event. Either in order or named. If unnamed and two entries, tend and event are assumed.
#'
#' @export
#' @return returns list with interval endpoints and event
#' @family Formula Interpretation
ColossusCoxSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Error: Too few entries in Cox survival object")
  }
  if (length(args) > 3) {
    stop("Error: Too many entries in Cox survival object")
  }
  tstart <- "%trunc%"
  tend <- "%trunc%"
  event <- "NULL"
  indx <- pmatch(argName[argName != ""], c("tstart", "tend", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Error: Argument '%s' not matched in survival object",
      argName[argName != ""][indx == 0L]
    ), domain = NA)
  }
  if (all(argName == "")) {
    # If none have names, assume (start, end, event) or (end, event)
    if (length(args) == 2) {
      tstart <- "%trunc%"
      tend <- args[[1]]
      event <- args[[2]]
    } else if (length(args) == 3) {
      tstart <- args[[1]]
      tend <- args[[2]]
      event <- args[[3]]
    } else {
      stop("Error: Incorrect number of arguments to survival object")
    }
  } else if (any(argName == "")) {
    # start by directly assigning what is available
    if ("tstart" %in% argName) {
      tstart <- args$tstart
    }
    if ("tend" %in% argName) {
      tend <- args$tend
    }
    if ("event" %in% argName) {
      event <- args$event
    }
    # now determine what is left
    # start by guessing by position
    if (length(args) == 2) {
      # either end/event or start/event
      if (event != "NULL") {
        # assume end/event
        tstart <- "%trunc%"
        tend <- args[[1]]
      } else {
        event <- args[[2]]
      }
    } else if (length(args) == 3) {
      if ((tstart == "%trunc%") && (argName[1] == "")) {
        tstart <- args[[1]]
      }
      if ((tend == "%trunc%") && (argName[2] == "")) {
        tend <- args[[2]]
      }
      if ((event == "NULL") && (argName[3] == "")) {
        event <- args[[3]]
      }
    } else {
      stop("Error: Incorrect number of arguments to survival object")
    }
  } else {
    if ("tstart" %in% argName) {
      tstart <- args$tstart
    }
    if ("tend" %in% argName) {
      tend <- args$tend
    }
    if ("event" %in% argName) {
      event <- args$event
    }
  }
  list("tstart" = tstart, "tend" = tend, "event" = event)
}

#' Interprets basic cox survival formula RHS with strata
#'
#' \code{ColossusCoxStrataSurv} assigns and interprets interval columns for cox model with stratification
#'
#' @param ... entries for a cox survival object with strata, tstart, tend, event, and strata. Either in order or named. If unnamed, the last is assumed to be strata and the rest are passed to ColossusCoxSurv.
#'
#' @noRd
#' @return returns list with interval endpoints, event, and strata column
#' @family Formula Interpretation
ColossusCoxStrataSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 3) {
    stop("Error: Too few entries in Cox Strata survival object")
  }
  if (length(args) > 4) {
    stop("Error: Too many entries in Cox Strata survival object")
  }
  strata <- "NULL"
  # Is stata a named entry?
  if ("strata" %in% argName) {
    strata <- parse_literal_string(args$strata)
    res <- do.call(ColossusCoxSurv, args[names(args) != "strata"])
  } else if (all(argName == "")) {
    strata <- parse_literal_string(args[[length(args)]])
    res <- do.call(ColossusCoxSurv, args[seq_len(length(args) - 1)])
  } else {
    if (argName[length(args)] == "") {
      strata <- parse_literal_string(args[[length(args)]])
      res <- do.call(ColossusCoxSurv, args[seq_len(length(args) - 1)])
    } else {
      stop("Error: Final entry of Cox Strata object was not named correctly")
    }
  }
  res[["strata"]] <- strata
  res
}

#' Interprets basic fine-gray survival formula RHS
#'
#' \code{ColossusFineGraySurv} assigns and interprets interval columns for fine-gray model
#'
#' @param ... entries for a Fine_Gray object, tstart, tend, event, and weighting. Either in order or named. If unnamed, the last is assumed to be weighting and the rest are passed to ColossusCoxSurv.
#'
#' @noRd
#' @return returns list with interval endpoints, event, and weighting columns
#' @family Formula Interpretation
ColossusFineGraySurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 3) {
    stop("Error: Too few entries in FineGray survival object")
  }
  if (length(args) > 4) {
    stop("Error: Too many entries in FineGray survival object")
  }
  weight <- "NULL"
  # Is stata a named entry?
  if ("weight" %in% argName) {
    weight <- args$weight
    res <- do.call(ColossusCoxSurv, args[names(args) != "weight"])
  } else if (all(argName == "")) {
    weight <- args[[length(args)]]
    res <- do.call(ColossusCoxSurv, args[seq_len(length(args) - 1)])
  } else {
    if (argName[length(args)] == "") {
      weight <- args[[length(args)]]
      res <- do.call(ColossusCoxSurv, args[seq_len(length(args) - 1)])
    } else {
      stop("Error: Final entry of FineGray object was not named correctly")
    }
  }
  res["weight"] <- weight
  res
}

#' Interprets basic fine-gray survival formula RHS with strata
#'
#' \code{ColossusFineGraySurv} assigns and interprets interval columns for fine-gray model and stratification
#'
#' @param ... entries for a Fine_Gray object with strata, tstart, tend, event, strata, and weighting. Either in order or named. If unnamed, the last is assumed to be weighting and the rest are passed to ColossusCoxStrataSurv.
#'
#' @noRd
#' @return returns list with interval endpoints, event, strata, and weighting columns
#' @family Formula Interpretation
ColossusFineGrayStrataSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 4) {
    stop("Error: Too few entries in FineGray Strata survival object")
  }
  if (length(args) > 5) {
    stop("Error: Too many entries in FineGray Strata survival object")
  }
  weight <- "NULL"
  # Is stata a named entry?
  if ("weight" %in% argName) {
    weight <- args$weight
    res <- do.call(ColossusCoxStrataSurv, args[names(args) != "weight"])
  } else if (all(argName == "")) {
    weight <- args[[length(args)]]
    res <- do.call(ColossusCoxStrataSurv, args[seq_len(length(args) - 1)])
  } else {
    if (argName[length(args)] == "") {
      weight <- args[[length(args)]]
      res <- do.call(ColossusCoxStrataSurv, args[seq_len(length(args) - 1)])
    } else {
      stop("Error: Final entry of FineGray object was not named correctly")
    }
  }
  res["weight"] <- weight
  res
}

# --------------------------------------------------------------------------------- #
#' Interprets basic poisson survival formula RHS
#'
#' \code{ColossusPoisSurv} assigns and interprets interval columns for poisson model.
#' This functions is called using the arguments for Poisson or Poisson_Strata in the
#' right-hand side of the formula. Uses an person-year column, number of events, and
#' any strata columns. The first two are expected to be in order or named: pyr and event.
#' Anything beyond the event name is assumed to be strata if Poisson_Strata is used.
#'
#' @param ... entries for a Poisson object with or without strata, pyr, event, and any strata columns. Either in order or named. The first two are assumed to be pyr and event, the rest assumed to be strata columns
#'
#' @export
#' @return returns list with duration, strata if used, and event
#' @family Formula Interpretation
ColossusPoisSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Error: Too few entries in Poisson survival object")
  }
  pyr <- "NULL"
  event <- "NULL"
  strata <- "NULL"
  indx <- pmatch(argName[argName != ""], c("pyr", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Error: Argument '%s' not matched in survival object",
      argName[argName != ""][indx == 0L]
    ), domain = NA)
  }
  if (all(argName == "")) {
    # If none have names, assume (pyr, event, and all strata)
    if (length(args) == 2) {
      pyr <- args[[1]]
      event <- args[[2]]
    } else {
      pyr <- args[[1]]
      event <- args[[2]]
      strata <- unlist(args[3:length(args)], use.names = FALSE)
    }
  } else if (any(argName == "")) {
    # start by directly assigning what is available
    if ("pyr" %in% argName) {
      pyr <- args$pyr
    }
    if ("event" %in% argName) {
      event <- args$event
    }
    # now determine what is left
    # start by guessing by position
    if (length(args) == 2) {
      if (event != "NULL") {
        pyr <- args[[1]]
      } else {
        event <- args[[2]]
      }
    } else {
      if ((pyr == "NULL") && (argName[1] == "")) {
        pyr <- args[[1]]
      }
      if ((event == "NULL") && (argName[2] == "")) {
        event <- args[[2]]
      }
      strata <- unlist(args[3:length(args)], use.names = FALSE)
    }
  } else {
    pyr <- args$pyr
    event <- args$event
  }
  list("pyr" = pyr, "event" = event, "strata" = strata)
}

# --------------------------------------------------------------------------------- #
#' Interprets basic case-control survival formula RHS with no grouping
#'
#' \code{ColossusCaseConSurv} assigns and interprets interval columns for case-control model without grouping.
#'
#' @param ... entries for a Case-Control object, only the event column
#'
#' @noRd
#' @return returns list with event
#' @family Formula Interpretation
ColossusCaseConSurv <- function(...) {
  # expects one arguement, as the event column
  args <- list(...)
  argName <- names(args)
  if (length(args) < 1) {
    stop("Error: Too few entries in Case-Control survival object")
  }
  if (length(args) > 1) {
    stop("Error: Too many entries in Case-Control survival object")
  }
  if ((argName[[1]] == "event") || is.null(argName[[1]])) {
    event <- args[[1]]
  } else {
    stop("Error: Entry must be unnamed or named 'event'")
  }
  list("event" = event)
}

#' Interprets basic case-control survival formula RHS with grouping by risk group
#'
#' \code{ColossusCaseConTimeSurv} assigns and interprets interval columns for case-control model with grouping by risk group.
#'
#' @param ... entries for a Case-Control object, with risk group info. Similar to basic Cox survival, interval start, end, and the event column. Either named or in order.
#'
#' @noRd
#' @return returns list with event
#' @family Formula Interpretation
ColossusCaseConTimeSurv <- function(...) {
  # expects time and event arguements
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Error: Too few entries in Case-Control survival object matched on time")
  }
  if (length(args) > 3) {
    stop("Error: Too many entries in Case-Control survival object matched on time")
  }
  tstart <- "%trunc%"
  tend <- "%trunc%"
  event <- "NULL"
  indx <- pmatch(argName[argName != ""], c("tstart", "tend", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Error: Argument '%s' not matched in survival object",
      argName[argName != ""][indx == 0L]
    ), domain = NA)
  }
  if (all(argName == "")) {
    # If none have names, assume (start, end, event) or (end, event)
    if (length(args) == 2) {
      tstart <- "%trunc%"
      tend <- args[[1]]
      event <- args[[2]]
    } else if (length(args) == 3) {
      tstart <- args[[1]]
      tend <- args[[2]]
      event <- args[[3]]
    } else {
      stop("Error: Incorrect number of arguments to survival object")
    }
  } else if (any(argName == "")) {
    # start by directly assigning what is available
    if ("tstart" %in% argName) {
      tstart <- args$tstart
    }
    if ("tend" %in% argName) {
      tend <- args$tend
    }
    if ("event" %in% argName) {
      event <- args$event
    }
    # now determine what is left
    # start by guessing by position
    if (length(args) == 2) {
      # either end/event or start/event
      if (event != "NULL") {
        # assume end/event
        tstart <- "%trunc%"
        tend <- args[[1]]
      } else {
        event <- args[[2]]
      }
    } else if (length(args) == 3) {
      if ((tstart == "%trunc%") && (argName[1] == "")) {
        tstart <- args[[1]]
      }
      if ((tend == "%trunc%") && (argName[2] == "")) {
        tend <- args[[2]]
      }
      if ((event == "NULL") && (argName[3] == "")) {
        event <- args[[3]]
      }
    } else {
      stop("Error: Incorrect number of arguments to survival object")
    }
  } else {
    if ("tstart" %in% argName) {
      tstart <- args$tstart
    }
    if ("tend" %in% argName) {
      tend <- args$tend
    }
    if ("event" %in% argName) {
      event <- args$event
    }
  }
  list("tstart" = tstart, "tend" = tend, "event" = event)
}

#' Interprets basic case-control survival formula RHS with grouping by strata
#'
#' \code{ColossusCaseConStrataSurv} assigns and interprets interval columns for case-control model with grouping by strata.
#'
#' @param ... entries for a Case-Control object, with strata. Expects an event column and strata column, either named or in order.
#'
#' @noRd
#' @return returns list with event
#' @family Formula Interpretation
ColossusCaseConStrataSurv <- function(...) {
  # expects two arguements, the event and strata columns
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Error: Too few entries in Case-Control survival object matched on strata")
  }
  if (length(args) > 2) {
    stop("Error: Too many entries in Case-Control survival object matched on strata")
  }
  indx <- pmatch(argName[argName != ""], c("strata", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Error: Argument '%s' not matched in survival object",
      argName[argName != ""][indx == 0L]
    ), domain = NA)
  }
  #
  if ("event" %in% argName) {
    event <- args$event
    strata <- parse_literal_string(args[names(args) != "event"][[1]])
  } else if ("strata" %in% argName) {
    strata <- parse_literal_string(args$strata)
    event <- args[names(args) != "strata"][[1]]
  } else {
    event <- args[[1]]
    strata <- parse_literal_string(args[[2]])
  }
  #
  list("event" = event, "strata" = strata)
}

#' Interprets basic case-control survival formula RHS with grouping by strata and risk group
#'
#' \code{ColossusCaseConTimeStrataSurv} assigns and interprets interval columns for case-control model with grouping by strata and risk group.
#'
#' @param ... entries for a Case-Control object, with strata and interval info. Similar to the stratified Cox survival object. Expects interval start, end, event, and strata in order or named. Removes strata and calls the standard Case-Control survival object grouped by risk group.
#'
#' @noRd
#' @return returns list with event
#' @family Formula Interpretation
ColossusCaseConTimeStrataSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 3) {
    stop("Error: Too few entries in Case-Control survival object matched on strata and time")
  }
  if (length(args) > 4) {
    stop("Error: Too many entries in Case-Control survival object matched on strata and time")
  }
  strata <- "NULL"
  # Is stata a named entry?
  if ("strata" %in% argName) {
    strata <- parse_literal_string(args$strata)
    res <- do.call(ColossusCaseConTimeSurv, args[names(args) != "strata"])
  } else if (all(argName == "")) {
    strata <- parse_literal_string(args[[length(args)]])
    res <- do.call(ColossusCaseConTimeSurv, args[seq_len(length(args) - 1)])
  } else {
    if (argName[length(args)] == "") {
      strata <- parse_literal_string(args[[length(args)]])
      res <- do.call(ColossusCaseConTimeSurv, args[seq_len(length(args) - 1)])
    } else {
      stop("Error: Final entry of Case Control Strata and Time object was not named correctly")
    }
  }
  res[["strata"]] <- strata
  res
}

#' Interprets basic logistic survival formula RHS with no grouping
#'
#' \code{ColossusLogitSurv} assigns and interprets columns for trials and events in logistic model with no grouping.
#'
#' @param ... entries for a Logistic object, trials and events. trials not provided assumes one trial per row.
#'
#' @export
#' @return returns list with event
#' @family Formula Interpretation
ColossusLogitSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  #
  if (length(args) < 1) {
    stop("Error: Too few entries in Logistic survival object")
  }
  if (length(args) > 2) {
    stop("Error: Too many entries in Logistic survival object")
  }
  trials <- "NULL"
  events <- "NULL"
  indx <- pmatch(argName[argName != ""], c("trials", "event", "events"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Error: Argument '%s' not matched in survival object",
      argName[argName != ""][indx == 0L]
    ), domain = NA)
  }
  if (("event" %in% argName) && ("events" %in% argName)) {
    stop("Error: Cannot provide both 'event' and 'events'")
  }
  if (all(argName == "")) {
    if (length(args) == 1) {
      # should just be the events
      # assume only one trial
      events <- args[[1]]
      trials <- "CONST"
    } else {
      # should be events and then cases
      events <- args[[2]]
      trials <- args[[1]]
    }
  } else {
    # atleast one is named
    if (any(argName == "")) {
      # one is named
      if ("trials" %in% argName) {
        trials <- args$trials
        events <- args[names(args) != "trials"][[1]]
      }
      if ("event" %in% argName) {
        events <- args$event
        trials <- args[names(args) != "event"][[1]]
      }
      if ("events" %in% argName) {
        events <- args$events
        trials <- args[names(args) != "events"][[1]]
      }
    } else {
      trials <- args$trials
      events <- args[names(args) != "trials"][[1]]
    }
  }
  list("event" = events, "trials" = trials)
}
