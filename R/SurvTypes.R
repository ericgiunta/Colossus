#' Interprets a Colossus formula and makes necessary changes to data
#'
#' \code{get_form} uses a formula and data.table, to fully describe the model
#' for a Colossus regression function.
#'
#' @inheritParams R_template
#'
#' @return returns a class fully describing the model and the updated data
#' @export
#' @family Formula Interpretation
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
#' formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~ loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()
#' model <- get_form(formula, df)
get_form <- function(formula, df) {
  surv_obj <- format(formula[[2]])
  model_obj <- paste(format(formula[[3]]), collapse = " ")
  surv_obj <- gsub(" ", "", surv_obj)
  model_obj <- gsub(" ", "", model_obj)
  term_n <- c()
  tform <- c()
  names <- c()
  gmix_theta <- 0
  gmix_term <- c()
  modelform <- "NONE"
  surv_model_type <- "NONE"
  tstart <- "NONE"
  tend <- "NONE"
  pyr <- "NONE"
  event <- "NONE"
  strata <- "NONE"
  weight <- "NONE"
  null <- FALSE
  #
  # split the survival model type
  second_split <- lapply(strsplit(surv_obj, ""), function(x) which(x == "("))[[1]]
  if (length(second_split) == 0) {
    stop("Error: The left hand side did not contain (")
  }
  surv_type <- tolower(substr(surv_obj, 1, second_split - 1))
  surv_paras <- substr(surv_obj, second_split + 1, nchar(surv_obj) - 1)
  surv_paras <- strsplit(surv_paras, ",")[[1]]
  # assign survival model values
  surv_model_type <- surv_type
  surv_para_list <- list()
  for (i in 1:length(surv_paras)) {
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
  if (surv_type %in% c("cox","coxph")) {
    res <- do.call(ColossusCoxSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
  } else if (surv_type %in% c("cox_strata","coxph_strata")) {
    res <- do.call(ColossusCoxStrataSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    strata <- res$strata
  } else if (surv_type %in% c("finegray","fg")) {
    res <- do.call(ColossusFineGraySurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    weight <- res$weight
  } else if (surv_type %in% c("finegray_strata", "fg_strata")) {
    res <- do.call(ColossusFineGrayStrataSurv, surv_para_list)
    tstart <- res$tstart
    tend <- res$tend
    event <- res$event
    strata <- res$strata
    weight <- res$weight
  } else if (surv_type %in% c("poisson","pois")) {
    res <- do.call(ColossusPoisSurv, surv_para_list)
    pyr <- res$pyr
    event <- res$event
  } else if (surv_type %in% c("poisson_strata", "pois_strata")) {
    res <- do.call(ColossusPoisSurv, surv_para_list)
    pyr <- res$pyr
    event <- res$event
    strata <- res$strata
  } else {
    stop("Error: Invalid survival model type")
  }
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
      model_paras <- strsplit(model_paras, ",")[[1]]
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
        if (substr(model_paras[subterm_i], 1, 7) == "factor(") {
          # baseline is set by using factor(column;baseline=level)
          col_name <- substr(model_paras[subterm_i], 8, nchar(model_paras[subterm_i]) - 1)
          base_split <- lapply(strsplit(col_name, ""), function(x) which(x == ";"))[[1]]
          if (length(base_split) == 1) {
            col_names <- substr(col_name, 1, base_split - 1)
            base_line <- substr(col_name, base_split + 1 + 9, nchar(col_name))
            base_level <- paste(col_names, base_line, sep = "_")
            #
            val <- factorize(df, col_names)
            df <- val$df
            col_name <- val$cols
            if (base_level %in% names(df)) {
              df <- df[, (base_level) := NULL]
              col_name <- col_name[!col_name == base_level]
            } else {
              warning(paste("Warning: Baseline level ", base_level, " not found.", sep = ""))
            }
          } else {
            val <- factorize(df, col_name)
            df <- val$df
            col_name <- val$cols
          }
        } else {
          col_name <- c(model_paras[subterm_i])
          if ("CONST" == model_paras[subterm_i]) {
            if ("CONST" %in% names(df)) {
              # fine
            } else {
              df$CONST <- 1
            }
          }
        }
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
            names <- c(names, col)
            tform <- c(tform, model_term)
            term_n <- c(term_n, term_num)
          }
        }
      }
    } else if (model_type %in% modelform_acceptable) {
      if (model_type %in% c("m", "me", "multiplicative", "multiplicative-excess")) {
        model_type <- "M"
      } else if (model_type %in% c("a", "additive")) {
        model_type <- "A"
      } else if (model_type %in% c("pa", "product-additive")) {
        model_type <- "PA"
      } else if (model_type %in% c("pae", "product-additive-excess")) {
        model_type <- "PAE"
      } else if (model_type %in% c("gmix", "geometric-mixture")) {
        model_type <- "GMIX"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        model_paras <- tolower(strsplit(model_paras, ",")[[1]])
        gmix_theta <- as.numeric(model_paras[1])
        gmix_term <- ifelse(model_paras[2:length(model_paras)] == "e", 1, 0)
      } else if (model_type %in% c("gmix-r", "relative-geometric-mixture")) {
        model_type <- "GMIX-R"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        model_paras <- tolower(strsplit(model_paras, ",")[[1]])
        gmix_theta <- as.numeric(model_paras[1])
      } else if (model_type %in% c("gmix-e", "excess-geometric-mixture")) {
        model_type <- "GMIX-E"
        model_paras <- substr(right_model_terms[term_i], third_split + 1, nchar(right_model_terms[term_i]) - 1)
        model_paras <- tolower(strsplit(model_paras, ",")[[1]])
        gmix_theta <- as.numeric(model_paras[1])
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
      modelform <- "M"
    } else if (modelform == "GMIX-R") {
      modelform <- "GMIX"
      gmix_term <- rep(0, term_tot)
    } else if (modelform == "GMIX-E") {
      modelform <- "GMIX"
      gmix_term <- rep(1, term_tot)
    }
  }
  if (grepl("cox", surv_model_type)) {
    model <- coxmodel(tstart, tend, event, strata, weight, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df)
  } else if (grepl("finegray", surv_model_type)) {
    model <- coxmodel(tstart, tend, event, strata, weight, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df)
  } else if (grepl("pois", surv_model_type)) {
    model <- poismodel(pyr, event, strata, null, term_n, tform, names, modelform, gmix_term, gmix_theta, c(), c(), df)
  } else {
    stop("Bad survival model type passed")
  }
  list(
    "model" = model, "data" = df
  )
}

#' Interprets basic cox survival formula RHS
#'
#' \code{ColossusCoxSurv} assigns and interprets interval columns for cox model
#'
#' @inheritParams R_template
#'
#' @return returns list with interval endpoints and event
#' @family Formula Interpretation
ColossusCoxSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Too few entries in Cox survival object")
  }
  if (length(args) > 3) {
    stop("Too many entries in Cox survival object")
  }
  tstart <- "%trunc%"
  tend <- "%trunc%"
  event <- "NULL"
  indx <- pmatch(argName[argName != ""], c("tstart", "tend", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Argument '%s' not matched in survival object",
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
      stop("Incorrect number of arguements to survival object")
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
      stop("Incorrect number of arguements to survival object")
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
#' @inheritParams R_template
#'
#' @return returns list with interval endpoints, event, and strata column
#' @family Formula Interpretation
ColossusCoxStrataSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 3) {
    stop("Too few entries in Cox Strata survival object")
  }
  if (length(args) > 4) {
    stop("Too many entries in Cox Strata survival object")
  }
  strata <- "NULL"
  # Is stata a named entry?
  if ("strata" %in% argName) {
    strata <- args$strata
    res <- do.call(ColossusCoxSurv, args[names(args) != "strata"])
  } else if (all(argName == "")) {
    strata <- args[[length(args)]]
    res <- do.call(ColossusCoxSurv, args[1:length(args) - 1])
  } else {
    if (argName[length(args)] == "") {
      strata <- args[[length(args)]]
      res <- do.call(ColossusCoxSurv, args[1:length(args) - 1])
    } else {
      stop("Final entry of Cox Strata object was not named correctly")
    }
  }
  res["strata"] <- strata
  res
}

#' Interprets basic fine-gray survival formula RHS
#'
#' \code{ColossusFineGraySurv} assigns and interprets interval columns for fine-gray model
#'
#' @inheritParams R_template
#'
#' @return returns list with interval endpoints, event, and weighting columns
#' @family Formula Interpretation
ColossusFineGraySurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 3) {
    stop("Too few entries in FineGray survival object")
  }
  if (length(args) > 4) {
    stop("Too many entries in FineGray survival object")
  }
  weight <- "NULL"
  # Is stata a named entry?
  if ("weight" %in% argName) {
    weight <- args$weight
    res <- do.call(ColossusCoxSurv, args[names(args) != "weight"])
  } else if (all(argName == "")) {
    weight <- args[[length(args)]]
    res <- do.call(ColossusCoxSurv, args[1:length(args) - 1])
  } else {
    if (argName[length(args)] == "") {
      weight <- args[[length(args)]]
      res <- do.call(ColossusCoxSurv, args[1:length(args) - 1])
    } else {
      stop("Final entry of FineGray object was not named correctly")
    }
  }
  res["weight"] <- weight
  res
}

#' Interprets basic fine-gray survival formula RHS with strata
#'
#' \code{ColossusFineGraySurv} assigns and interprets interval columns for fine-gray model and stratification
#'
#' @inheritParams R_template
#'
#' @return returns list with interval endpoints, event, strata, and weighting columns
#' @family Formula Interpretation
ColossusFineGrayStrataSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 4) {
    stop("Too few entries in FineGray Strata survival object")
  }
  if (length(args) > 5) {
    stop("Too many entries in FineGray Strata survival object")
  }
  weight <- "NULL"
  # Is stata a named entry?
  if ("weight" %in% argName) {
    weight <- args$weight
    res <- do.call(ColossusCoxStrataSurv, args[names(args) != "weight"])
  } else if (all(argName == "")) {
    weight <- args[[length(args)]]
    res <- do.call(ColossusCoxStrataSurv, args[1:length(args) - 1])
  } else {
    if (argName[length(args)] == "") {
      weight <- args[[length(args)]]
      res <- do.call(ColossusCoxStrataSurv, args[1:length(args) - 1])
    } else {
      stop("Final entry of FineGray object was not named correctly")
    }
  }
  res["weight"] <- weight
  res
}

# --------------------------------------------------------------------------------- #
#' Interprets basic poisson survival formula RHS
#'
#' \code{ColossusPoisSurv} assigns and interprets interval columns for poisson model
#'
#' @inheritParams R_template
#'
#' @return returns list with duration, strata if used, and event
#' @family Formula Interpretation
ColossusPoisSurv <- function(...) {
  args <- list(...)
  argName <- names(args)
  if (length(args) < 2) {
    stop("Too few entries in Poisson survival object")
  }
  if (length(args) > 3) {
    stop("Too few entries in Poisson survival object")
  }
  pyr <- "NULL"
  event <- "NULL"
  strata <- "NULL"
  indx <- pmatch(argName[argName != ""], c("pyr", "event"), nomatch = 0L)
  if (any(indx == 0L)) {
    stop(gettextf(
      "Argument '%s' not matched in survival object",
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
