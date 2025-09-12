#' Calculates and plots martingale residuals with a named dose column
#'
#' \code{CoxMartingale} uses user provided data, columns, and identifier to create plots
#'
#' @inheritParams R_template
#' @family Plotting Functions
#' @return saves the plots in the current directory and returns a list of data-tables plotted
#' @noRd
CoxMartingale <- function(verbose, df, time1, time2, event0, e, t, ch, dnames, plot_name, age_unit, studyID) {
  IDS <- base <- res <- doses <- NULL
  if (verbose >= 3) {
    message("Note: Plotting Martingale Residuals") # nocov
  }
  time_s <- df[, get(time1)]
  time_e <- df[, get(time2)]
  ch_fun <- approxfun(x = t, y = ch, rule = 2)
  ch_e <- ch_fun(time_e)
  ch_s <- ch_fun(time_s)
  e_i <- df[, get(event0)]
  table_out <- list()
  for (cov_i in seq_len(length(dnames))) {
    dname <- dnames[cov_i]
    if (verbose >= 3) {
      message(paste("Note: Martingale Plot: ", dname, sep = "")) # nocov
    }
    if (studyID %in% names(df)) {
      dfr <- data.table::data.table(
        "Risks" = e$Risks, "ch_e" = ch_e,
        "ch_s" = ch_s, "e" = e_i,
        "IDS" = unlist(df[, studyID, with = FALSE], use.names = FALSE),
        "time" = time_e, "cov" = unlist(df[, dname, with = FALSE],
          use.names = FALSE
        )
      )
      name_temp <- names(dfr)
      for (i in seq_len(length(name_temp))) {
        if (grepl("cov", name_temp[i], fixed = TRUE)) {
          data.table::setnames(dfr, name_temp[i], c("cov"))
        } else if (grepl("IDS", name_temp[i], fixed = TRUE)) {
          data.table::setnames(dfr, name_temp[i], c("IDS"))
        }
      }
      dfr$res <- dfr$e - (dfr$ch_e - dfr$ch_s) * dfr$Risks
      Martingale_Error <- dfr[, lapply(.SD, sum), by = IDS]
      times <- dfr[, lapply(.SD, max), by = IDS]
    } else {
      dfr <- data.table::data.table(
        "Risks" = e$Risks, "ch_e" = ch_e,
        "ch_s" = ch_s, "e" = e_i, "time" = time_e,
        "cov" = unlist(df[, dname, with = FALSE], use.names = FALSE)
      )
      name_temp <- names(dfr)
      for (i in seq_len(length(name_temp))) {
        if (grepl("cov", name_temp[i], fixed = TRUE)) {
          data.table::setnames(dfr, name_temp[i], c("cov"))
        } else if (grepl("IDS", name_temp[i], fixed = TRUE)) {
          data.table::setnames(dfr, name_temp[i], c("IDS"))
        }
      }
      dfr$res <- dfr$e - (dfr$ch_e) * dfr$Risks
      Martingale_Error <- dfr
      times <- dfr
    }
    dft <- data.table::data.table(
      "cov_max" = times$cov,
      "time_max" = times$time,
      "res_sum" = Martingale_Error$res,
      "event" = Martingale_Error$e
    )
    table_out[[dname]] <- dft
  }
  if (studyID %in% names(df)) {
    dfr <- data.table::data.table(
      "Risks" = e$Risks, "ch_e" = ch_e,
      "ch_s" = ch_s, "e" = e_i,
      "IDS" = unlist(df[, studyID, with = FALSE],
        use.names = FALSE
      ),
      "time" = time_e
    )
    name_temp <- names(dfr)
    for (i in seq_len(length(name_temp))) {
      if (grepl("IDS", name_temp[i], fixed = TRUE)) {
        data.table::setnames(dfr, name_temp[i], c("IDS"))
      }
    }
    dfr$res <- dfr$e - (dfr$ch_e - dfr$ch_s) * dfr$Risks
    Martingale_Error <- dfr[, lapply(.SD, sum), by = IDS]
    times <- dfr[, lapply(.SD, max), by = IDS]
  } else {
    dfr <- data.table::data.table(
      "Risks" = e$Risks, "ch_e" = ch_e,
      "ch_s" = ch_s, "e" = e_i, "time" = time_e
    )
    name_temp <- names(dfr)
    for (i in seq_len(length(name_temp))) {
      if (grepl("IDS", name_temp[i], fixed = TRUE)) {
        data.table::setnames(dfr, name_temp[i], c("IDS"))
      }
    }
    dfr$res <- dfr$e - (dfr$ch_e) * dfr$Risks
  }
  dft <- data.table::data.table(
    "time_max" = dfr$time, "res_sum" = dfr$res,
    "event" = dfr$e
  )
  table_out[["survival_time"]] <- dft
  return(table_out)
}

#' Calculates and plots survival plots of the estimated baseline
#'
#' \code{CoxSurvival} uses user provided data, columns, and identifier to create plots
#'
#' @inheritParams R_template
#' @family Plotting Functions
#' @noRd
#' @return saves the plots in the current directory and returns a list of data-tables plotted
CoxSurvival <- function(t, h, ch, surv, plot_name, verbose, time_lims, age_unit) {
  if (verbose >= 3) {
    message("Note: Plotting Survival Curves") # nocov
  }
  table_out <- list()
  dft <- data.table::data.table("t" = t, "h" = h, "ch" = ch, "surv" = surv)
  data.table::setkeyv(dft, "t")
  dft <- dft[(t >= time_lims[1]) & (t <= time_lims[2]), ]
  table_out[["standard"]] <- dft
  Ls <- log(surv)
  Lls_u <- log(-Ls)
  Lt_u <- log(t)
  dft <- data.table::data.table("t" = Lt_u, "s" = Lls_u)
  table_out[["log"]] <- dft
  return(table_out)
}

#' Calculates and plots Kaplan-Meier survival plots
#'
#' \code{CoxKaplanMeier} uses user provided data, columns, and identifier to create plots, plots the kaplan-meier survival and log(time) vs log(-log(survival))
#'
#' @inheritParams R_template
#' @family Plotting Functions
#' @noRd
#' @return saves the plots in the current directory and returns a list of data-tables plotted
CoxKaplanMeier <- function(verbose, studyID, names, df, event0, time1, time2, tu, term_n, tform, a_n, er, modelform, control, keep_constant, plot_type, age_unit, model_control = list()) {
  if (verbose >= 3) {
    message("Note: Plotting Kaplan-Meier Curve") # nocov
  }
  model_control <- Def_model_control(model_control)
  base <- NULL
  all_names <- unique(names)
  dfc <- match(names, all_names)
  df_study <- df[, lapply(.SD, max), by = studyID]
  dfend <- df_study[get(event0) == 1, ]
  t_num <- length(unlist(unique(df_study[, studyID, with = FALSE]),
    use.names = FALSE
  ))
  # number of unique individuals in total set
  t_t <- c(0)
  n_t <- c(1)
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]),
    use.names = FALSE
  )) # all event times
  tu_s <- c(0.0, tu)
  for (i_t in seq_len(length(tu))) {
    i <- tu[i_t]
    df0 <- dfend[get(time2) <= i, ] # set of all intervals prior to this point in lower setv
    df0 <- df0[(get(time2) > tu_s[i_t]), ]
    df1 <- df_study[(get(time2) > tu_s[i_t]), ]
    t_ev <- sum(df0[, get(event0)])
    # number of intervals with event in lower set prior to the time point
    t_num <- nrow(df1)
    if (t_ev > 0) {
      temp <- (t_num - t_ev) / t_num
      t_t <- c(t_t, i)
      n_t <- c(n_t, n_t[length(n_t)] * temp)
    }
  }
  table_out <- list()
  dft <- data.table::data.table("t_t" = t_t, "n_t" = n_t)
  table_out[["kaplin-meier"]] <- dft
  return(table_out)
}

#' Calculates and plots relative risk
#'
#' \code{CoxRisk} uses user provided data, columns, and identifier to create plots of risk by covariate value for each column
#'
#' @inheritParams R_template
#' @family Plotting Functions
#' @noRd
#' @return saves the plots in the current directory and returns a list of data-tables plotted
CoxRisk <- function(verbose, df, event0, time1, time2, names, term_n, tform, a_n, modelform, control, keep_constant, plot_type, b, er, model_control = list()) {
  fir_KM <- 0
  model_control <- Def_model_control(model_control)
  dfend <- df[get(event0) == 1, ]
  ce <- c(time1, time2, event0)
  all_names <- unique(names)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  tu <- unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE)
  table_out <- list()
  for (fir_KM in seq_len(length(names))) {
    lfir <- c(names[fir_KM])
    uniq <- unlist(unique(df[, lfir, with = FALSE]), use.names = FALSE)
    der_iden <- fir_KM - 1
    model_control$risk <- TRUE
    model_control$unique_values <- length(uniq)
    e <- Plot_Omnibus_transition(
      term_n, tform, a_n, dfc, x_all, 0,
      der_iden, modelform,
      control,
      as.matrix(df[, ce, with = FALSE]), tu,
      keep_constant, term_tot, c(0), c(0),
      model_control
    )
    if ("Failure" %in% names(e)) {
      stop(e)
    }
    x <- e$x
    y <- e$y
    dft <- data.table::data.table("x" = x, "y" = y)
    table_out[[names[fir_KM]]] <- dft
  }
  return(table_out)
}

#' Calculates and plots survival curves for each unique value of the stratification column
#'
#' \code{CoxStratifiedSurvival} uses user provided data, columns, and identifier to calculate the survival fraction for each strata
#'
#' @inheritParams R_template
#' @family Plotting Functions
#' @noRd
#' @return saves the plots in the current directory and returns a list of data-tables plotted
CoxStratifiedSurvival <- function(verbose, df, event0, time1, time2, names, term_n, tform, a_n, er, modelform, control, keep_constant, plot_type, strat_col, time_lims, age_unit, model_control = list()) {
  dfend <- df[get(event0) == 1, ]
  uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
    use.names = FALSE
  ))
  for (i in seq_along(uniq)) {
    df0 <- dfend[get(strat_col) == uniq[i], ]
    tu0 <- unlist(unique(df0[, time2, with = FALSE]), use.names = FALSE)
    if (length(tu0) == 0) {
      warning(paste("Warning: no events for strata group:", uniq[i],
        sep = " "
      ))
      df <- df[get(strat_col) != uniq[i], ]
    }
  }
  uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
    use.names = FALSE
  ))
  if (control$verbose >= 3) {
    message(paste("Note:", length(uniq), " strata used", sep = " ")) # nocov
  }
  data.table::setkeyv(df, c(strat_col, event0, time2, time1))
  ce <- c(time1, time2, event0, strat_col)
  model_control <- Def_model_control(model_control)
  base <- NULL
  ce <- c(time1, time2, event0, strat_col)
  all_names <- unique(names)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  tu <- unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE)
  uniq <- sort(unlist(unique(df[, strat_col, with = FALSE]),
    use.names = FALSE
  ))
  control$maxiters <- c(-1, -1)
  control$guesses <- 1
  e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names,
    term_n, tform, keep_constant,
    a_n, modelform, control,
    strat_col = strat_col,
    model_control = list("strata" = TRUE)
  )
  a_n <- e$beta_0
  plot_name <- plot_type[2]
  tt <- c()
  tsurv <- c()
  categ <- c()
  if (verbose >= 3) {
    message(paste("Note: Starting Stratification: Calculation")) # nocov
  }
  model_control$surv <- TRUE
  model_control$strata <- TRUE
  e <- Plot_Omnibus_transition(
    term_n, tform, a_n, dfc, x_all, 0, 0,
    modelform, control,
    as.matrix(df[, ce, with = FALSE]),
    tu, keep_constant, term_tot, uniq, c(0),
    model_control
  )
  for (col_i in seq_along(uniq)) {
    if (verbose >= 3) {
      message(paste(
        "Note: Starting Stratification calculation ",
        col_i
      )) # nocov
    }
    col_u <- uniq[col_i]
    t <- c()
    h <- c()
    ch <- c()
    surv <- c()
    dft <- data.table::data.table("time" = tu, "base" = e$baseline[, col_i])
    for (i in tu) {
      if ((i <= time_lims[2]) && (i >= time_lims[1])) {
        t <- c(t, i)
        temp <- sum(dft[time < i, base])
        ch <- c(ch, temp)
        if (length(h) == 0) {
          h <- c(temp)
        } else {
          h <- c(h, ch[length(ch)] - ch[length(ch)])
        }
        surv <- c(surv, exp(-1 * temp))
      }
    }
    tt <- c(tt, t)
    tsurv <- c(tsurv, surv)
    categ <- c(categ, rep(paste(col_u), length(t)))
  }
  dft <- data.table::data.table("t" = tt, "surv" = tsurv, "cat_group" = categ)
  sbreaks <- c()
  slabels <- c()
  for (i in seq_len(length(uniq))) {
    sbreaks <- c(sbreaks, paste(uniq[i]))
    slabels <- c(slabels, paste("For ", strat_col, "=", uniq[i], sep = ""))
  }
  table_out <- list()
  table_out[["stratified_survival"]] <- dft
  return(table_out)
}

#' Calculates Schoenfeld residuals for a Cox Proportional Hazards regression and plots
#'
#' \code{RunCox_Schoenfeld_Residual} uses user provided data, time/event columns, vectors specifying the model, and options to calculate the residuals
#'
#' @inheritParams R_template
#'
#' @return saves the plots in the current directory and returns a list of data-tables plotted
#' @family Plotting Functions
#' @noRd
#' @importFrom rlang .data
PlotCox_Schoenfeld_Residual <- function(df, time1, time2, event0, names, term_n, tform, keep_constant, a_n, modelform, control, age_unit, plot_name, model_control = list()) {
  data.table::setkeyv(df, c(event0, time2, time1))
  model_control <- Def_model_control(model_control)
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (length(tu) == 0) {
    stop("Error: no events")
  }
  if (control$verbose >= 3) {
    message(paste("Note: ", length(tu), " risk groups", sep = "")) # nocov
  }
  all_names <- unique(names)
  dfc <- match(names, all_names)
  term_tot <- max(term_n) + 1
  x_all <- as.matrix(df[, all_names, with = FALSE])
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  control <- Def_Control(control)
  df <- Replace_Missing(df, all_names, 0.0, control$verbose)
  model_control$schoenfeld <- TRUE
  res_list <- Plot_Omnibus_transition(
    term_n, tform, a_n, dfc, x_all, 0, 0,
    modelform, control,
    as.matrix(df[, ce, with = FALSE]),
    tu, keep_constant, term_tot, c(0),
    c(0), model_control
  )
  res <- res_list$residual
  res_scaled <- res_list$scaled
  table_out <- list()
  for (cov in seq_len(length(a_n))) {
    if (keep_constant[cov] == 0) {
      cov_res <- cov - sum(head(keep_constant, cov))
      y <- unlist(res[, cov_res], use.names = FALSE)
      y_scale <- unlist(res_scaled[, cov_res], use.names = FALSE)
      dft <- data.table::data.table("time" = tu, "y" = y)
      dft$y_scale <- y_scale
      table_out[[names[cov]]] <- dft
    }
  }
  return(table_out)
}
