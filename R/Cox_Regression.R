#' Performs Cox Proportional Hazards regression using the omnibus function
#'
#' \code{RunCoxRegression_Omnibus} uses user provided data, time/event columns,
#'   vectors specifying the model, and options to control the convergence
#'   and starting positions. Has additional options for starting with several
#'   initial guesses, using stratification, multiplicative loglinear 1-term,
#'   competing risks, and calculation without derivatives
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @noRd
#' @family Cox Wrapper Functions
#' @importFrom rlang .data
RunCoxRegression_Omnibus <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", cens_weight = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  #
  to_remove <- c("CONST", "%trunc%")
  to_keep <- c(time1, time2, event0, names)
  if (model_control$cr == TRUE) {
    to_keep <- c(to_keep, cens_weight)
  }
  if (model_control$strata == TRUE) {
    to_keep <- c(to_keep, strat_col)
  }
  to_keep <- unique(to_keep)
  to_keep <- to_keep[!to_keep %in% to_remove]
  to_keep <- to_keep[to_keep %in% names(df)]
  df <- df[, to_keep, with = FALSE]
  tend <- Sys.time()
  #
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  ## Cox regression only uses intervals which contain an event time
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  # remove rows that end before first event
  df <- df[get(time2) >= tu[1], ]
  # remove rows that start after the last event
  df <- df[get(time1) <= tu[length(tu)], ]
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$basic == TRUE) {
    if (all(unique(tform) == c("loglin"))) {
      # good
    } else {
      if (control$verbose >= 2) {
        warning("Warning: Basic loglinear model used, but atleast one subterm was not loglin. Subterms all set to loglin")
      }
      tform <- rep("loglin", length(tform))
    }
    if (length(unique(term_n)) > 1) {
      if (control$verbose >= 2) {
        warning("Warning: Basic loglinear model used, but more than one term number used. Term numbers all set to 0")
      }
      term_n <- rep(0, length(term_n))
    }
    if (modelform != "M") {
      if (control$verbose >= 2) {
        warning("Warning: Basic loglinear model used, but multiplicative model not used. Modelform corrected")
      }
      modelform <- "M"
    }
  }
  if (model_control$linear_err == TRUE) {
    if (all(sort(unique(tform)) != c("loglin", "plin"))) {
      stop("Error: Linear ERR model used, but term formula wasn't only loglin and plin")
    }
    if (sum(tform == "plin") > 1) {
      stop("Error: Linear ERR model used, but more than one plin element was used")
    }
    if (length(unique(term_n)) > 1) {
      if (control$verbose >= 2) {
        warning("Warning: Linear ERR model used, but more than one term number used. Term numbers all set to 0")
      }
      term_n <- rep(0, length(term_n))
    }
    if (modelform != "M") {
      if (control$verbose >= 2) {
        warning("Warning: Linear ERR model used, but multiplicative model not used. Modelform corrected")
      }
      modelform <- "M"
    }
  }
  if (model_control$cr == TRUE) {
    if (cens_weight %in% names(df)) {
      # good
    } else {
      stop("Error: censoring weight column not in the dataframe.")
    }
  } else {
    df[[cens_weight]] <- 1
  }
  if (model_control$strata == FALSE) {
    data.table::setkeyv(df, c(event0, time2, time1))
    uniq <- c(0)
    ce <- c(time1, time2, event0)
  } else {
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
    ce <- c(time1, time2, event0, strat_col)
  }
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (control$verbose >= 3) {
    message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
  }
  all_names <- unique(names)
  if (!model_control$null) {
    df <- Replace_Missing(df, all_names, 0.0, control$verbose)
    # make sure any constant 0 columns are constant
    for (i in seq_along(keep_constant)) {
      if ((keep_constant[i] == 0) && (names[i] %in% names(df))) {
        if (names[i] != "CONST") {
          if (min(df[[names[i]]]) == max(df[[names[i]]])) {
            keep_constant[i] <- 1
            if (control$verbose >= 2) {
              warning(paste("Warning: element ", i,
                " with column name ", names[i],
                " was set constant",
                sep = ""
              ))
            }
          }
        }
      }
    }
    if (min(keep_constant) > 0) {
      stop("Error: Atleast one parameter must be free")
    }
  }
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  a_ns <- c()
  for (i in a_n) {
    a_ns <- c(a_ns, i)
  }
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
    e <- cox_ph_Omnibus_Bounds_transition(
      term_n, tform, a_ns,
      dfc, x_all, 0,
      modelform, control, as.matrix(df[, ce, with = FALSE]), tu,
      keep_constant, term_tot, uniq, df[[cens_weight]], model_control,
      cons_mat, cons_vec
    )
    if ("Status" %in% names(e)) {
      if (e$Status != "PASSED") {
        stop(e$Status)
      }
    }
  } else {
    res <- Check_Iters(control, a_n)
    control <- res$control
    a_n <- res$a_n
    if (model_control$null) {
      a_ns <- matrix(a_ns)
    } else {
      a_ns <- matrix(a_ns, nrow = length(control$maxiters) - 1, byrow = TRUE)
    }
    e <- cox_ph_Omnibus_transition(
      term_n, tform, a_ns, dfc, x_all, 0,
      modelform, control, as.matrix(df[, ce, with = FALSE]), tu,
      keep_constant, term_tot, uniq, df[[cens_weight]], model_control,
      cons_mat, cons_vec
    )
    if ("Status" %in% names(e)) {
      if (is.nan(e$LogLik)) {
        stop(e$Status)
      }
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Cox"
  e$modelcontrol <- model_control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}

#' Calculates hazard ratios for a reference vector
#'
#' \code{RunCoxRegression} uses user provided data,  vectors specifying the model,
#' and options to calculate relative risk for every row in the provided data
#'
#' @noRd
#' @inheritParams R_template
#' @family Plotting Wrapper Functions
#' @return returns a list of the final results
Cox_Relative_Risk <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), model_control = list()) {
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  all_names <- unique(names)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  model_control$risk_subset <- TRUE
  e <- Plot_Omnibus_transition(
    term_n, tform, a_n, dfc, x_all, 0, 0,
    modelform, control, matrix(c(0)),
    c(1), keep_constant, term_tot, c(0),
    c(0), model_control
  )
  return(e)
}

#' Performs Cox Proportional Hazard model plots
#'
#' \code{RunCoxPlots} uses user provided data, time/event columns,
#' vectors specifying the model, and options to choose and save plots
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return saves the plots in the current directory and returns the data used for plots
#' @family Plotting Wrapper Functions
RunCoxPlots <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), plot_options = list(), model_control = list()) {
  names(plot_options) <- tolower(names(plot_options))
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  plot_options$verbose <- Check_Verbose(plot_options$verbose)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  if (plot_options$verbose >= 3) {
    message("Note: Starting Plot Function") # nocov
  }
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  data.table::setkeyv(df, c(event0, time2, time1))
  base <- NULL
  plot_type <- plot_options$type
  if (plot_options$verbose >= 3) {
    message("Note: Getting Plot Info") # nocov
  }
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  if (plot_options$verbose >= 3) {
    message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
  }
  if ("type" %in% names(plot_options)) {
    # fine
  } else {
    stop("Error: Plot type not given")
  }
  if ("age_unit" %in% names(plot_options)) {
    # fine
  } else {
    plot_options$age_unit <- "unitless"
  }
  if ("strat_haz" %in% names(plot_options)) {
    if (plot_options$strat_haz) {
      if ("strat_col" %in% names(plot_options)) {
        if (plot_options$strat_col %in% names(df)) {
          # fine
          strat_col <- plot_options$strat_col
          if (!is.null(levels(df[[strat_col]]))) {
            # The column is a factor, so we can convert to numbers
            factor_lvl <- levels(df[[strat_col]])
            df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
          } else if (is(typeof(df[[strat_col]]), "character")) {
            df[[strat_col]] <- factor(df[[strat_col]])
            factor_lvl <- levels(df[[strat_col]])
            df[[strat_col]] <- as.integer(factor(df[[strat_col]], levels = factor_lvl)) - 1
          }
        } else {
          stop("Error: Stratification Column not in dataframe")
        }
      } else {
        stop("Error: Stratification Column not given")
      }
    }
  } else {
    plot_options$strat_haz <- FALSE
  }
  if ("martingale" %in% names(plot_options)) {
    if (plot_options$martingale) {
      if ("cov_cols" %in% names(plot_options)) {
        for (cov_i in seq_along(plot_options$cov_cols)) {
          dose_col <- unlist(plot_options$cov_cols,
            use.names = FALSE
          )[cov_i]
          if (dose_col %in% names(df)) {
            # fine
          } else {
            stop("Error: Covariate column " +
              dose_col + " is not in the dataframe")
          }
        }
      } else {
        stop("Error: dose column not given")
      }
      if ("studyid" %in% names(plot_options)) {
        if (plot_options$studyid %in% names(df)) {
          # fine
        } else {
          stop("Error: ID column is not in the dataframe")
        }
      } else {
        stop("Error: ID column not given")
      }
    }
  } else {
    plot_options$martingale <- FALSE
  }
  if ("km" %in% names(plot_options)) {
    if (plot_options$km) {
      if ("studyid" %in% names(plot_options)) {
        if (plot_options$studyid %in% names(df)) {
          # fine
        } else {
          stop("Error: ID column is not in the dataframe")
        }
      } else {
        stop("Error: ID column not given")
      }
    }
  }
  model_control <- Def_model_control(model_control)
  if (tolower(plot_type[1]) == "surv") {
    if ("time_lims" %in% names(plot_options)) {
      # fine
    } else {
      plot_options$time_lims <- c(min(tu), max(tu))
    }
  }
  for (iden_col in c("verbose", "martingale", "surv_curv", "strat_haz", "km")) {
    if (iden_col %in% names(plot_options)) {
      # fine
    } else {
      plot_options[iden_col] <- FALSE
    }
  }
  plot_options$verbose <- Check_Verbose(plot_options$verbose)
  control <- Def_Control(control)
  verbose <- data.table::copy(plot_options$verbose)
  maxiterc <- data.table::copy(control$maxiter)
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  control$maxiters <- c(-1, -1)
  control$guesses <- 1
  e <- RunCoxRegression_Omnibus(
    df, time1, time2, event0, names, term_n,
    tform, keep_constant, a_n, modelform,
    control, model_control
  )
  control$maxiter <- maxiterc
  b <- e$beta_0
  er <- e$Standard_Deviation
  plot_table <- list()
  if (tolower(plot_type[1]) == "surv") {
    if (verbose >= 3) {
      message("Note: starting ph_plot") # nocov
    }
    if (plot_options$strat_haz == FALSE) {
      if (verbose >= 3) {
        message("Note: nonStratified survival curve calculation") # nocov
      }
      model_control$surv <- TRUE
      e <- Plot_Omnibus_transition(
        term_n, tform, a_n, dfc, x_all, 0, 0,
        modelform, control,
        as.matrix(df[, ce, with = FALSE]), tu,
        keep_constant, term_tot, c(0), c(0),
        model_control
      )
      t <- c()
      h <- c()
      ch <- c()
      surv <- c()
      if (verbose >= 3) {
        message("Note: writing survival data") # nocov
      }
      dft <- data.table::data.table(
        "time" = tu, "base" = e$baseline,
        "basehaz" = e$standard_error
      )
      for (i in tu) {
        t <- c(t, i)
        temp <- sum(dft[time < i, base])
        ch <- c(ch, temp)
        if (length(h) == 0) {
          h <- c(temp)
        } else {
          h <- c(h, ch[length(ch)] - ch[length(ch) - 1])
        }
        surv <- c(surv, exp(-1 * temp))
      }
      age_unit <- plot_options$age_unit
      if (plot_options$martingale == TRUE) {
        plot_table <- CoxMartingale(
          verbose, df, time1, time2, event0,
          e, t, ch,
          plot_options$cov_cols,
          plot_type[2], age_unit, plot_options$studyid
        )
      }
      if (plot_options$surv_curv == TRUE) {
        plot_table <- CoxSurvival(
          t, h, ch, surv, plot_type[2], verbose,
          plot_options$time_lims, age_unit
        )
      }
    } else {
      age_unit <- plot_options$age_unit
      if (verbose >= 3) {
        message("Note: Stratified survival curve calculation") # nocov
      }
      if (plot_options$surv_curv == TRUE) {
        model_control$strata <- TRUE
        plot_table <- CoxStratifiedSurvival(
          verbose, df, event0,
          time1, time2,
          all_names, term_n, tform, a_n, er,
          modelform, control, keep_constant, plot_type,
          plot_options$strat_col, plot_options$time_lims, age_unit
        )
      }
    }
    if (plot_options$km == TRUE) {
      plot_table <- CoxKaplanMeier(
        verbose, plot_options$studyid,
        all_names, df, event0, time1, time2, tu, term_n,
        tform, a_n, er, modelform,
        control, keep_constant, plot_type, age_unit
      )
    }
  } else if (tolower(plot_type[1]) == "risk") {
    plot_table <- CoxRisk(
      verbose, df, event0, time1, time2,
      names, term_n, tform,
      a_n, modelform, control, keep_constant,
      plot_type, b, er
    )
  } else if (tolower(plot_type[1]) == "schoenfeld") {
    age_unit <- plot_options$age_unit
    plot_table <- PlotCox_Schoenfeld_Residual(
      df, time1,
      time2, event0, names, term_n,
      tform, keep_constant, a_n, modelform,
      control, age_unit, plot_type[2]
    )
  }
  return(plot_table)
}

#' Performs Cox Proportional Hazards regression using the omnibus function with multiple column realizations
#'
#' \code{RunCoxRegression_Omnibus_Multidose} uses user provided data, time/event columns,
#'       vectors specifying the model, and options to control the convergence
#'       and starting positions. Used for 2DMC column uncertainty methods.
#'       Returns optimized parameters, log-likelihood, and standard deviation for each realization.
#'       Has additional options for using stratification,
#'       multiplicative loglinear 1-term,
#'       competing risks, and calculation without derivatives
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return returns a list of the final results for each realization
#' @family Cox Wrapper Functions
#' @importFrom rlang .data
RunCoxRegression_Omnibus_Multidose <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", realization_columns = matrix(c("temp00", "temp01", "temp10", "temp11"), nrow = 2), realization_index = c("temp0", "temp1"), control = list(), strat_col = "null", cens_weight = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  #
  control <- Def_Control(control)
  model_control <- Def_model_control(model_control)
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
  }
  #
  to_remove <- c("CONST", "%trunc%")
  to_keep <- c(time1, time2, event0, names, realization_index, as.vector(realization_columns))
  if (model_control$cr == TRUE) {
    to_keep <- c(to_keep, cens_weight)
  }
  if (model_control$strata == TRUE) {
    to_keep <- c(to_keep, strat_col)
  }
  to_keep <- unique(to_keep)
  to_keep <- to_keep[!to_keep %in% to_remove]
  to_keep <- to_keep[to_keep %in% names(df)]
  df <- df[, to_keep, with = FALSE]
  #
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  ## Cox regression only uses intervals which contain an event time
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  # remove rows that end before first event
  df <- df[get(time2) >= tu[1], ]
  # remove rows that start after the last event
  df <- df[get(time1) <= tu[length(tu)], ]
  #
  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$cr == TRUE) {
    if (cens_weight %in% names(df)) {
      # good
    } else if (length(cens_weight) < nrow(df)) {
      stop("Error: censoring weight column not in the dataframe.")
    }
  } else {
    df[[cens_weight]] <- 1
  }
  if ("MCML" %in% names(model_control)) {
    # fine
  } else {
    model_control$MCML <- FALSE
  }
  if (model_control$strata == FALSE) {
    data.table::setkeyv(df, c(event0, time2, time1))
    uniq <- c(0)
    ce <- c(time1, time2, event0)
  } else {
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
      message(paste("Note:", length(uniq), " strata used",
        sep = " "
      )) # nocov
    }
    data.table::setkeyv(df, c(strat_col, event0, time2, time1))
    ce <- c(time1, time2, event0, strat_col)
  }
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  if (control$verbose >= 3) {
    message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
  }
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  # replace_missing equivalent for the realization columns
  if (length(realization_index) == length(realization_columns[, 1])) {
    # pass
  } else {
    # the number of columns per realization does not match the number of indexes provided
    stop(paste("Error:", length(realization_index),
      " column indexes provided, but ",
      length(realization_columns[, 1]),
      " rows of realizations columns provided",
      sep = " "
    ))
  }
  if (all(realization_index %in% all_names)) {
    # pass
  } else {
    stop(paste("Error: Atleast one realization column provided was not used in the model", sep = " "))
  }
  #  all_names <- unique(c(all_names, as.vector(realization_columns)))
  dose_names <- unique(as.vector(realization_columns))
  if (all(dose_names %in% names(df))) {
    # pass
  } else {
    stop(paste("Error: Atleast one realization column provided was not in the data.table", sep = " "))
  }
  dfc <- match(names, all_names)
  dose_cols <- matrix(match(realization_columns, dose_names), nrow = nrow(realization_columns))
  dose_index <- match(realization_index, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  dose_all <- as.matrix(df[, dose_names, with = FALSE])
  e <- cox_ph_multidose_Omnibus_transition(
    term_n, tform, a_n,
    as.matrix(dose_cols, with = FALSE), dose_index, dfc, x_all, dose_all,
    0, modelform, control,
    as.matrix(df[, ce, with = FALSE]), tu,
    keep_constant, term_tot, uniq, df[[cens_weight]], model_control,
    cons_mat, cons_vec
  )
  if ("Status" %in% names(e)) {
    if (e$Status != "PASSED") {
      stop(e$Status)
    }
  }
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$modelformula <- modelform
  e$Parameter_Lists$keep_constant <- keep_constant
  if (model_control$MCML) {
    e$Survival_Type <- "Cox_Multidose"
  } else {
    e$Survival_Type <- "Cox_Multidose"
  }
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}

#' Calculates the likelihood curve for a cox model directly
#'
#' \code{CoxCurveSolver} solves the confidence interval for a cox model, starting at the optimum point and
#' iteratively optimizing end-points of intervals. Intervals updated using the bisection method.
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return returns a list of the final results
#' @family Cox Wrapper Functions
#' @importFrom rlang .data
CoxCurveSolver <- function(df, time1 = "%trunc%", time2 = "%trunc%", event0 = "event", names = c("CONST"), term_n = c(0), tform = "loglin", keep_constant = c(0), a_n = c(0), modelform = "M", control = list(), strat_col = "null", cens_weight = "null", model_control = list(), cons_mat = as.matrix(c(0)), cons_vec = c(0)) {
  func_t_start <- Sys.time()
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df) # nocov
      },
      error = function(e) { # nocov
        df <- data.table(df) # nocov
      }
    )
  }
  control <- Def_Control(control)
  model_control$log_bound <- TRUE
  model_control <- Def_model_control(model_control)
  if (typeof(a_n) != "list") {
    a_n <- list(a_n)
  }
  #
  to_remove <- c("CONST", "%trunc%")
  to_keep <- c(time1, time2, event0, names)
  if (model_control$cr == TRUE) {
    to_keep <- c(to_keep, cens_weight)
  }
  if (model_control$strata == TRUE) {
    to_keep <- c(to_keep, strat_col)
  }
  to_keep <- unique(to_keep)
  to_keep <- to_keep[!to_keep %in% to_remove]
  to_keep <- to_keep[to_keep %in% names(df)]
  df <- df[, to_keep, with = FALSE]
  #
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

  if ("CONST" %in% names) {
    if ("CONST" %in% names(df)) {
      # fine
    } else {
      df$CONST <- 1
    }
  }
  if (model_control$linear_err == TRUE) {
    if (all(sort(unique(tform)) != c("loglin", "plin"))) {
      stop("Error: Linear ERR model used, but term formula wasn't only loglin and plin")
    }
    if (sum(tform == "plin") > 1) {
      stop("Error: Linear ERR model used, but more than one plin element was used")
    }
    if (length(unique(term_n)) > 1) {
      if (control$verbose >= 2) {
        warning("Warning: Linear ERR model used, but more than one term number used. Term numbers all set to 0")
      }
      term_n <- rep(0, length(term_n))
    }
    if (modelform != "M") {
      if (control$verbose >= 2) {
        warning("Warning: Linear ERR model used, but multiplicative model not used. Modelform corrected")
      }
      modelform <- "M"
    }
  }
  if (model_control$cr == TRUE) {
    if (cens_weight %in% names(df)) {
      # good
    } else {
      stop("Error: censoring weight column not in the dataframe.")
    }
  } else {
    df[[cens_weight]] <- 1
  }
  if (model_control$strata == FALSE) {
    data.table::setkeyv(df, c(event0, time2, time1))
    uniq <- c(0)
    ce <- c(time1, time2, event0)
  } else {
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
    ce <- c(time1, time2, event0, strat_col)
  }
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (control$verbose >= 3) {
    message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
  }
  all_names <- unique(names)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  # make sure any constant 0 columns are constant
  for (i in seq_along(keep_constant)) {
    if ((keep_constant[i] == 0) && (names[i] %in% names(df))) {
      if (names[i] != "CONST") {
        if (min(df[[names[i]]]) == max(df[[names[i]]])) {
          keep_constant[i] <- 1
          if (control$verbose >= 2) {
            warning(paste("Warning: element ", i,
              " with column name ", names[i],
              " was set constant",
              sep = ""
            ))
          }
        }
      }
    }
  }
  if (min(keep_constant) > 0) {
    stop("Error: Atleast one parameter must be free")
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
  if ("alpha" %in% names(model_control)) {
    model_control["qchi"] <- qchisq(1 - model_control[["alpha"]], df = 1) / 2
  } else {
    model_control["alpha"] <- 0.05
    model_control["qchi"] <- qchisq(1 - model_control[["alpha"]], df = 1) / 2
  }
  a_ns <- matrix(a_ns, nrow = length(control$maxiters) - 1, byrow = TRUE)
  para_number <- model_control$para_number
  keep_constant[para_number] <- 1
  if (min(keep_constant) == 1) {
    model_control["single"] <- TRUE
  }
  e <- cox_ph_Omnibus_CurveSearch_transition(
    term_n, tform, a_ns,
    dfc, x_all, 0,
    modelform, control, as.matrix(df[, ce, with = FALSE]), tu,
    keep_constant, term_tot, uniq, df[[cens_weight]], model_control,
    cons_mat, cons_vec
  )
  e$Parameter_Lists$names <- names
  e$Parameter_Lists$keep_constant <- keep_constant
  e$Parameter_Lists$modelformula <- modelform
  e$Survival_Type <- "Cox"
  e$modelcontrol <- model_control
  func_t_end <- Sys.time()
  e$RunTime <- func_t_end - func_t_start
  return(e)
}
