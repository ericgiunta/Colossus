#' Splits a string by commas, skipping over parenthesis sections
#'
#' \code{nested_split} splits by comma, keeps parenthesis sections together
#'
#' @noRd
#' @param total_string the complete string to split
#' @family Data Cleaning Functions
#' @return returns a vector of substrings
nested_split <- function(total_string) {
  # start by doing best
  sub_str <- strsplit(total_string, ",|(?>\\(.*?\\).*?\\K(,|$))", perl = TRUE)[[1]]
  # There may still be some areas where the nested section was split
  final_split <- c("")
  temp_count <- 0
  for (i in sub_str) {
    temp_count <- temp_count + str_count(i, "\\(") - str_count(i, "\\)") # check to see if the current selections have the same number
    if (final_split[length(final_split)] == "") {
      final_split[length(final_split)] <- i # first is just added
    } else {
      final_split[length(final_split)] <- paste(final_split[length(final_split)], i, sep = ",") # accumulate
    }
    if (temp_count == 0) {
      final_split <- c(final_split, "") # start next if the counts match
    }
  }
  if (final_split[length(final_split)] != "") {
    final_split # safety to avoid error cases
  } else {
    final_split[seq_len(length(final_split) - 1)] # the last should be empty
  }
}

#' converts a string of a vector/list to vector/list
#'
#' \code{parse_literal_string} converts the string of a vector/list to a vector/list
#'
#' @noRd
#' @param string the string to convert to vector/list
#' @family Data Cleaning Functions
#' @return returns a vector, list, or the original string
parse_literal_string <- function(string) {
  if (substr(string, 1, 2) == "c(") { # converts the string to a vector
    sub_str <- substr(string, 3, nchar(string) - 1)
    args <- nested_split(sub_str) # get the entries of the vector
    args <- gsub("\"", "", args) # remove string literals
    args <- unlist(lapply(args, parse_literal_string), use.names = F) # make sure every entry is processed to a final type
    return(args)
  } else if (substr(string, 1, 5) == "list(") { # converts string to list
    sub_str <- substr(string, 6, nchar(string) - 1)
    args <- nested_split(sub_str) # make sure every entry is processed
    #
    factor_list <- list()
    for (i in seq_along(args)) {
      para_cur <- args[i]
      para_break <- lapply(strsplit(para_cur, ""), function(x) which(x == "="))[[1]] # split into the name and value
      if (length(para_break) == 0) {
        # no name, just add to list
        factor_list[[i]] <- parse_literal_string(para_cur) # process the item to the final value
      } else {
        item_name <- substr(para_cur, 1, para_break - 1)
        item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
        item_name <- gsub("\"", "", item_name)
        item_value <- gsub("\"", "", item_value)
        factor_list[[item_name]] <- parse_literal_string(item_value) # process the item to the final value
      }
    }
    return(factor_list)
  }
  if (string != ".") { # '.' breaks the numeric search
    options(warn = -1)
    if (tolower(string) %in% c("true", "false")) { # check for boolean
      options(warn = 0)
      return(as.logical(string))
    }
    if (all(vapply(string, function(x) grepl("^[\\-]{0,1}[0-9]*\\.{0,1}[0-9]*$", x), logical(1))) || all(vapply(string, function(x) grepl("^[\\-]{0,1}[0-9]+e[\\-]{0,1}[0-9]+$", x), logical(1)))) {
      options(warn = 0) # checks for an integer, decimal, decimal places or scientific notation
      return(as.numeric(string))
    }
    options(warn = 0)
  }
  string
}

#' Automatically assigns missing values in listed columns
#'
#' \code{Replace_Missing} checks each column and fills in NA values
#'
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a filled datatable
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' df <- data.table::data.table(
#'   "UserID" = c(112, 114, 213, 214, 115, 116, 117),
#'   "Starting_Age" = c(18, 20, 18, 19, 21, 20, 18),
#'   "Ending_Age" = c(30, 45, NA, 47, 36, NA, 55),
#'   "Cancer_Status" = c(0, 0, 1, 0, 1, 0, 0)
#' )
#' df <- Replace_Missing(df, c("Starting_Age", "Ending_Age"), 70)
Replace_Missing <- function(df, name_list, msv, verbose = FALSE) {
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
  verbose <- Check_Verbose(verbose)
  if (is.na(msv)) {
    stop("Error: The missing-value replacement is also NA")
  }
  for (j in name_list) {
    #
    if (j %in% names(df)) {
      # fine
    } else {
      stop(paste("Error: ", j, " missing from column names", sep = ""))
    }
    if (sum(is.na(df[[j]]))) {
      data.table::set(df, which(is.na(df[[j]])), j, msv)
      if (verbose >= 3) {
        message(paste("Note: Column ", j, " had replaced values",
          sep = ""
        )) # nocov
      }
    }
  }
  return(df)
}

#' Automatically assigns missing control values
#'
#' \code{Def_Control} checks and assigns default values
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a filled list
Def_Control <- function(control) {
  control_def <- list(
    "verbose" = 0, "lr" = 0.75, "maxiter" = 20,
    "halfmax" = 5, "epsilon" = 1e-4,
    "deriv_epsilon" = 1e-4, "step_max" = 1.0,
    "change_all" = TRUE, "thres_step_max" = 100.0,
    "ties" = "breslow",
    "ncores" = as.numeric(detectCores())
  )
  names(control) <- tolower(names(control))
  if ((identical(Sys.getenv("TESTTHAT"), "true")) || (identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))) {
    control_def$ncores <- min(c(2, as.numeric(detectCores())))
  }
  if (Sys.getenv("ColossusOMP") == "") {
    syscheck <- System_Version()
    OpenMP <- syscheck[["OpenMP Enabled"]]
    if (!OpenMP) {
      Sys.setenv(ColossusOMP = "FALSE")
    } else {
      Sys.setenv(ColossusOMP = "TRUE")
      Sys.setenv(ColossusGCC = "TRUE")
      os <- syscheck[["Operating System"]]
      if (os == "linux") {
        cpp_compiler <- syscheck[["Default c++"]]
        if (cpp_compiler != "") {
          if (cpp_compiler == "gcc") {
            R_compiler <- syscheck[["R Compiler"]]
            if (R_compiler != "gcc") { # nocov
              Sys.setenv(ColossusGCC = "FALSE")
            }
          } else if (cpp_compiler == "clang") { # nocov
            Sys.setenv(ColossusGCC = "FALSE")
          }
        } else {
          R_compiler <- syscheck[["R Compiler"]]
          if (R_compiler != "gcc") { # nocov
            Sys.setenv(ColossusGCC = "FALSE")
          }
        }
      }
    }
  }
  if ("verbose" %in% names(control)) {
    control$verbose <- Check_Verbose(control$verbose)
  } else {
    control["verbose"] <- control_def["verbose"]
  }
  if (Sys.getenv("ColossusOMP") == "FALSE") {
    if (control["verbose"] > 1) {
      warning("Warning: OpenMP not detected, cores set to 1")
    }
    control$ncores <- 1 # nocov
  } else if ((Sys.getenv("R_COLOSSUS_NOT_CRAN") == "") && (Sys.getenv("ColossusGCC") == "FALSE")) {
    control$ncores <- 1 # nocov
    if (control["verbose"] > 1) {
      warning("Warning: linux machine not using gcc, cores set to 1. Set R_COLOSSUS_NOT_CRAN environemnt variable to skip check")
    }
  }
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
  return(control)
}

#' Automatically assigns missing model control values
#'
#' \code{Def_model_control} checks and assigns default values
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a filled list
Def_model_control <- function(control) {
  names(control) <- tolower(names(control))
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "surv",
    "schoenfeld", "risk",
    "risk_subset", "log_bound", "pearson", "deviance",
    "mcml", "observed_info", "time_risk",
    "logit_odds", "logit_ident", "logit_loglink"
  )
  for (nm in control_def_names) {
    if (nm %in% names(control)) {
      # fine
    } else {
      control[nm] <- FALSE
    }
  }
  #  if ((control["gradient"] == TRUE) & (control["constraint"] == TRUE)) {
  #    warning("Current constraints are not available with gradient descent methods. Only gradient descent will be applied")
  #    control["constaint"] <- FALSE
  #  }
  link_vec <- c(control$logit_odds, control$logit_ident, control$logit_loglink)
  if (sum(link_vec) == 0) {
    control["logit_odds"] <- TRUE
  } else if (sum(link_vec) > 1) {
    stop("Error: Multiple link functions used, only select one link function")
  }
  if (control["null"] == TRUE) {
    control["single"] <- TRUE
  }
  if ("step_size" %in% names(control)) {
    # fine
  } else {
    control["step_size"] <- 0.5
  }
  if ("unique_values" %in% names(control)) {
    # fine
  } else {
    control["unique_values"] <- 2
  }
  if ("gmix_theta" %in% names(control)) {
    # fine
  } else {
    control["gmix_theta"] <- 0.5
  }
  if ("gmix_term" %in% names(control)) {
    # fine
  } else {
    control["gmix_term"] <- c(0)
  }
  if ("conditional_threshold" %in% names(control)) {
    # fine
  } else {
    control["conditional_threshold"] <- 50
  }
  if (control[["log_bound"]]) {
    if ("qchi" %in% names(control)) {
      # fine
    } else {
      if ("alpha" %in% names(control)) {
        control["qchi"] <- qchisq(1 - control[["alpha"]], df = 1) / 2
      } else {
        control["alpha"] <- 0.05
        control["qchi"] <- qchisq(1 - control[["alpha"]], df = 1) / 2
      }
    }
    if ("para_num" %in% names(control)) {
      warning("Warning: para_num detected in model_control, did you mean para_number?")
      control["para_number"] <- control["para_num"]
    }
    if ("para_number" %in% names(control)) {
      # fine
    } else {
      control["para_number"] <- 1
    }
    if ("maxstep" %in% names(control)) {
      # fine
    } else {
      control["maxstep"] <- 10
    }
    if ("manual" %in% names(control)) {
      # fine
    } else {
      control["manual"] <- FALSE
    }
    if ("search_mult" %in% names(control)) {
      # fine
    } else {
      control["search_mult"] <- 1.0
    }
    if ("step_size" %in% names(control)) {
      # fine
    } else {
      control["search_mult"] <- 0.5
    }
  }
  if (control[["gradient"]]) {
    control_def_names <- c(
      "momentum", "adadelta", "adam"
    )
    for (nm in control_def_names) {
      if (nm %in% names(control)) {
        # fine
      } else {
        control[nm] <- FALSE
      }
    }
    # add different parameters
    if ("momentum_decay" %in% names(control)) {
      # fine
    } else {
      control["momentum_decay"] <- 0.9
    }
    if ("learning_decay" %in% names(control)) {
      # fine
    } else {
      control["learning_decay"] <- 0.999
    }
    if ("epsilon_decay" %in% names(control)) {
      # fine
    } else {
      control["epsilon_decay"] <- 1e-4
    }
    if (control["constraint"] == TRUE) {
      # needs to add constraint penalty
      if ("penalty_weight" %in% names(control)) {
        # fine
      } else {
        control["penalty_weight"] <- 1.0
      }
      if ("penalty_method" %in% names(control)) {
        # fine
      } else {
        control["penalty_method"] <- "sqr_error"
      }
    }
  }
  return(control)
}

#' Automatically checks the number of starting guesses
#'
#' \code{Check_Iters} checks the number of iterations and number of guesses, and corrects
#'
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a list with the corrected control list and a_n
Check_Iters <- function(control, a_n) {
  if ("maxiters" %in% names(control)) {
    if (length(control$maxiters) == length(a_n) + 1) {
      # all good, it matches
    } else {
      if (control$verbose >= 3) { # nocov
        message(paste("Note: Initial starts:", length(a_n),
          ", Number of iterations provided:",
          length(control$maxiters),
          ". Colossus requires one more iteration counts than number of guesses (for best guess)",
          sep = " "
        )) # nocov
      }
      if (length(control$maxiters) < length(a_n) + 1) {
        additional <- length(a_n) + 1 - length(control$maxiters)
        control$maxiters <- c(control$maxiters, rep(1, additional))
      } else {
        additional <- length(a_n) + 1
        control$maxiters <- control$maxiters[1:additional]
      }
    }
    control$guesses <- length(control$maxiters) - 1
  } else {
    control$guesses <- length(a_n)
    control$maxiters <- c(rep(1, length(a_n)), control$maxiter)
  }
  list("control" = control, "a_n" = a_n)
}

#' Calculates Full Parameter list for Special Dose Formula
#'
#' \code{Linked_Dose_Formula} Calculates all parameters for linear-quadratic and linear-exponential linked formulas
#'
#' @inheritParams R_template
#'
#' @return returns list of full parameters
#' @export
#' @examples
#' library(data.table)
#' tforms <- list("cov_0" = "quad", "cov_1" = "exp")
#' paras <- list("cov_0" = c(1, 3.45), "cov_1" = c(1.2, 4.5, 0.1))
#' full_paras <- Linked_Dose_Formula(tforms, paras)
#'
Linked_Dose_Formula <- function(tforms, paras, verbose = 0) {
  verbose <- Check_Verbose(verbose)
  full_paras <- list()
  for (nm in names(tforms)) {
    if (tforms[nm] == "quad") {
      plist <- unlist(paras[nm], use.names = FALSE)
      a0 <- plist[1]
      y <- plist[2]
      if (a0 < 0) {
        stop("Error: a0 argument was negative")
      }
      if (is.numeric(a0)) {
        # fine
      } else {
        stop("Error: a0 argument was not a number")
      }
      if (is.numeric(y)) {
        # fine
      } else {
        stop("Error: threshold argument was not a number") # nocov
      }
      a1 <- a0 / 2 / y
      b1 <- a0 * y / 2
      full_paras[[nm]] <- c(y, a0, a1, b1)
    } else if (tforms[nm] == "exp") {
      plist <- unlist(paras[nm], use.names = FALSE)
      a0 <- plist[1]
      y <- plist[2]
      b1 <- plist[3]
      if (a0 < 0) {
        stop("Error: a0 argument was negative")
      }
      if (is.numeric(a0)) {
        # fine
      } else {
        stop("Error: a0 argument was not a number") # nocov
      }
      if (is.numeric(y)) {
        # fine
      } else {
        stop("Error: threshold argument was not a number") # nocov
      }
      if (is.numeric(b1)) {
        # fine
      } else {
        stop("Error: exponential argument was not a number") # nocov
      }
      c1 <- log(a0) - log(b1) + b1 * y
      a1 <- a0 * y + exp(c1 - b1 * y)
      full_paras[[nm]] <- c(y, a0, a1, b1, c1)
    }
  }
  return(full_paras)
}

#' Calculates The Additional Parameter For a linear-exponential formula with known maximum
#'
#' \code{Linked_Lin_Exp_Para} Calculates what the additional parameter would be for a desired maximum
#'
#' @inheritParams R_template
#'
#' @return returns parameter used by Colossus
#' @export
#' @examples
#' library(data.table)
#' y <- 7.6
#' a0 <- 1.2
#' a1_goal <- 15
#' full_paras <- Linked_Lin_Exp_Para(y, a0, a1_goal)
#'
Linked_Lin_Exp_Para <- function(y, a0, a1_goal, verbose = 0) {
  verbose <- Check_Verbose(verbose)
  b1 <- 10
  lr <- 1.0
  if (a0 < 0) {
    stop("Error: a0 argument was negative")
  }
  if (a1_goal > y * a0) {
    # fine
  } else {
    stop("Error: goal is too low")
  }
  iter_i <- 0
  while (iter_i < 100) {
    iter_i <- iter_i + 1
    c1 <- log(a0 / b1) + b1 * y
    a1 <- a0 * y + exp(c1 - b1 * y)
    a_dif <- a1 - a1_goal
    if (abs(a_dif) < 1e-3) {
      break
    }
    da1 <- -1 / b1 * exp(c1 - b1 * y)
    db1 <- (a1_goal - a1) / da1
    if (-1 * db1 > b1) {
      db1 <- -0.9 * b1
    }
    b1 <- b1 + lr * db1
  }
  return(b1)
}

#' Splits a parameter into factors
#'
#' \code{factorize} uses user provided list of columns to define new parameter for each unique value and update the data.table.
#' Not for interaction terms
#'
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#' @examples
#' library(data.table)
#' a <- c(0, 1, 2, 3, 4, 5, 6)
#' b <- c(1, 2, 3, 4, 5, 6, 7)
#' c <- c(0, 1, 2, 1, 0, 1, 0)
#' df <- data.table::data.table("a" = a, "b" = b, "c" = c)
#' col_list <- c("c")
#' val <- factorize(df, col_list)
#' df <- val$df
#' new_col <- val$cols
#'
factorize <- function(df, col_list, verbose = 0) {
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
  verbose <- Check_Verbose(verbose)
  cols <- c()
  col0 <- names(df)
  tnum <- c()
  for (i in seq_len(length(col_list))) {
    col <- col_list[i]
    x <- sort(unlist(as.list(unique(df[, col, with = FALSE])),
      use.names = FALSE
    ))
    for (j in x) {
      newcol <- c(paste(col, j, sep = "_"))
      df[, newcol] <- 1 * (df[, col, with = FALSE] == j)
      cols <- c(cols, newcol)
      tnum <- c(tnum, i)
    }
  }
  cols <- setdiff(names(df), col0)
  if (length(col_list) > 1) {
    cols <- Check_Dupe_Columns(df, cols, rep(0, length(cols)), verbose, TRUE)
  }
  if (verbose >= 3) {
    message(paste("Note: Number of factors:", length(cols), sep = "")) # nocov
  }
  list("df" = df, "cols" = cols)
}

#' Defines the likelihood ratio test
#'
#' \code{Likelihood_Ratio_Test} uses two models and calculates the ratio
#'
#' @inheritParams R_template
#'
#' @return returns the score statistic
#' @export
#' @examples
#' library(data.table)
#' # In an actual example, one would run two seperate RunCoxRegression regressions,
#' #    assigning the results to e0 and e1
#' e0 <- list("name" = "First Model", "LogLik" = -120)
#' e1 <- list("name" = "New Model", "LogLik" = -100)
#' score <- Likelihood_Ratio_Test(e1, e0)
#'
Likelihood_Ratio_Test <- function(alternative_model, null_model) {
  if (("LogLik" %in% names(alternative_model)) && ("LogLik" %in% names(null_model))) {
    freedom <- length(alternative_model$beta_0) - length(null_model$beta_0)
    val <- 2 * (unlist(alternative_model["LogLik"], use.names = FALSE) - unlist(null_model["LogLik"], use.names = FALSE))
    pval <- pchisq(val, freedom)
    return(list("value" = val, "p value" = pval))
  } else {
    stop("Error: models input did not contain LogLik values")
  }
  return(NULL) # nocov
}

#' checks for duplicated column names
#'
#' \code{Check_Dupe_Columns} checks for duplicated columns, columns with the same values, and columns with single value. Currently not updated for multi-terms
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns the usable columns
Check_Dupe_Columns <- function(df, cols, term_n, verbose = 0, factor_check = FALSE) {
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
  verbose <- Check_Verbose(verbose)
  if (length(cols) > 1) {
    features_pair <- combn(cols, 2, simplify = FALSE) # list all column pairs
    terms_pair <- combn(term_n, 2, simplify = FALSE) # list all term pairs
    toRemove <- c() # init a vector to store duplicates
    for (pair_n in seq_len(length(features_pair))) {
      # put the pairs for testing into temp objects
      pair <- unlist(features_pair[pair_n])
      term <- unlist(terms_pair[pair_n])
      f1 <- pair[1]
      f2 <- pair[2]
      if (!(f1 %in% names(df))) {
        stop(paste("Error: ", f1, " not in data.table", sep = "")) # nocov
      }
      if (!(f2 %in% names(df))) {
        stop(paste("Error: ", f2, " not in data.table", sep = "")) # nocov
      }
      t1 <- term[1]
      t2 <- term[2]
      checked_factor <- TRUE
      if (factor_check) { #
        if ((is.numeric(df[[f1]])) && (is.numeric(df[[f2]]))) {
          if (sum(df[[f1]] * df[[f2]]) == 0) {
            checked_factor <- FALSE
          }
        } else if (is.numeric(df[[f1]]) != is.numeric(df[[f2]])) {
          checked_factor <- FALSE # nocov
        }
      }
      if ((t1 == t2) && (checked_factor)) {
        if (!(f1 %in% toRemove) && !(f2 %in% toRemove)) {
          if (all(df[[f1]] == df[[f2]])) { # test for duplicates
            if (verbose >= 2) {
              warning(paste("Warning: ", f1, " and ", f2,
                " are equal",
                sep = ""
              ))
            }
            toRemove <- c(toRemove, f2) # build the list of duplicates
          }
          if (min(df[[f2]]) == max(df[[f2]])) {
            if (min(df[[f2]]) == 0) {
              toRemove <- c(toRemove, f2) # remove zero values
            }
          }
        }
      }
    }
    newcol <- setdiff(cols, toRemove)
    if (length(newcol) == 1) {
      if (min(df[, newcol, with = FALSE]) == max(df[, newcol, with = FALSE])) {
        return(c()) # nocov
      } else {
        return(newcol)
      }
    }
    return(newcol)
  } else if (length(cols) == 1) {
    f1 <- cols[1]
    if (!(f1 %in% names(df))) {
      stop(paste("Error: ", f1, " not in data.table", sep = ""))
    }
    if (min(df[, cols, with = FALSE]) == max(df[, cols, with = FALSE])) {
      if (min(df[, cols, with = FALSE]) == 0) {
        return(c())
      } else {
        return(cols)
      }
    } else {
      return(cols)
    }
  } else {
    return(c())
  }
  return(c()) # nocov
}

#' Applies time duration truncation limits to create columns for Cox model
#'
#' \code{Check_Trunc} creates columns to use for truncation
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns the updated data and time period columns
Check_Trunc <- function(df, ce, verbose = 0) {
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
  verbose <- Check_Verbose(verbose)
  if (ce[1] %in% c("%trunc%", "right_trunc")) {
    if (ce[2] %in% c("%trunc%", "left_trunc")) {
      stop("Error: Both endpoints are truncated, not acceptable")
    }
    tname <- ce[2]
    tmin <- min(df[, get(tname)]) - 1
    if (!("right_trunc" %in% names(df))) {
      df[, ":="(right_trunc = tmin)]
    }
    ce[1] <- "right_trunc"
  } else if (ce[2] %in% c("%trunc%", "left_trunc")) {
    tname <- ce[1]
    tmax <- max(df[, get(tname)]) + 1
    if (!("left_trunc" %in% names(df))) {
      df[, ":="(left_trunc = tmax)]
    }
    ce[2] <- "left_trunc"
  }
  return(list("df" = df, "ce" = ce))
}

#' Applies time dependence to parameters
#'
#' \code{gen_time_dep} generates a new dataframe with time dependent covariates by applying a grid in time
#'
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#' # Adapted from the tests
#' a <- c(20, 20, 5, 10, 15)
#' b <- c(1, 2, 1, 1, 2)
#' c <- c(0, 0, 1, 1, 1)
#' df <- data.table::data.table("a" = a, "b" = b, "c" = c)
#' time1 <- "%trunc%"
#' time2 <- "a"
#' event <- "c"
#' control <- list(
#'   "lr" = 0.75, "maxiter" = -1, "halfmax" = 5, "epsilon" = 1e-9,
#'   "deriv_epsilon" = 1e-9, "step_max" = 1.0,
#'   "thres_step_max" = 100.0,
#'   "verbose" = FALSE, "ties" = "breslow", "double_step" = 1
#' )
#' grt_f <- function(df, time_col) {
#'   return((df[, "b"] * df[, get(time_col)])[[1]])
#' }
#' func_form <- c("lin")
#' df_new <- gen_time_dep(
#'   df, time1, time2, event, TRUE, 0.01, c("grt"), c(),
#'   c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 2
#' )
#' file.remove("test_new.csv")
#'
gen_time_dep <- function(df, time1, time2, event0, iscox, dt, new_names, dep_cols, func_form, fname, tform, nthreads = as.numeric(detectCores())) {
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
  dfn <- names(df)
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  dfn_same <- dfn[!(dfn %in% dep_cols)]
  dfn_dep <- c()
  for (i in seq_len(length(new_names))) {
    name0 <- paste(new_names[i], 0, sep = "_")
    name1 <- paste(new_names[i], 1, sep = "_")
    func <- func_form[i]
    df[, name0] <- lapply(func, function(f) f(df, time1))
    df[, name1] <- lapply(func, function(f) f(df, time2))
    dfn_dep <- c(dfn_dep, name0, name1)
  }
  if (length(new_names) != length(func_form)) {
    stop(paste("Error: new_names vector should be the same size as the list of functions applied",
      sep = ""
    ))
  }
  if (length(new_names) != length(tform)) {
    stop(paste("Error: new_names vector should be the same size as the list of interpolation method used",
      sep = ""
    ))
  }
  for (i in seq_len(length(tform))) {
    temp <- tform[i]
    if (temp != "lin") {
      a <- substr(temp, 1, 5)
      if (a != "step?") {
        stop(paste("Error: Interpolation method not recognized: ",
          temp,
          sep = ""
        ))
      }
    }
  }
  dfn_time <- c(time1, time2)
  dfn_event <- c(event0)
  dfn_same <- dfn_same[!(dfn_same %in% dfn_time)]
  dfn_same <- dfn_same[!(dfn_same %in% dfn_event)]
  dfend <- df[get(event0) == 1, ]
  tu <- sort(unlist(unique(dfend[, time2, with = FALSE]), use.names = FALSE))
  if (iscox) {
    df <- df[get(time2) >= min(tu), ]
    df <- df[get(time1) <= max(tu), ]
  }
  x_time <- as.matrix(df[, dfn_time, with = FALSE])
  x_dep <- as.matrix(df[, dfn_dep, with = FALSE])
  x_same <- as.matrix(df[, dfn_same, with = FALSE])
  x_event <- as.matrix(df[, dfn_event, with = FALSE])
  if (grepl(".csv", fname, fixed = TRUE)) {
    # fine
  } else {
    fname <- paste(fname, ".csv", sep = "_")
  }
  if ((identical(Sys.getenv("TESTTHAT"), "true")) || (identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))) {
    nthreads <- 2
  }
  Write_Time_Dep(
    x_time, x_dep, x_same, x_event, dt, fname,
    tform, tu, iscox, nthreads
  )
  df_new <- data.table::fread(fname,
    data.table = TRUE,
    header = FALSE, nThread = nthreads,
    col.names = c(time1, time2, new_names, dfn_same, event0)
  )
  data.table::setkeyv(df_new, c(event0, time2, time1))
  return(df_new)
}

#' Automates creating a date difference column
#'
#' \code{Date_Shift} generates a new dataframe with a column containing time difference in a given unit
#'
#' @inheritParams R_template
#' @param dcol0 list of starting month, day, and year
#' @param dcol1 list of ending month, day, and year
#' @family Data Cleaning Functions
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#' m0 <- c(1, 1, 2, 2)
#' m1 <- c(2, 2, 3, 3)
#' d0 <- c(1, 2, 3, 4)
#' d1 <- c(6, 7, 8, 9)
#' y0 <- c(1990, 1991, 1997, 1998)
#' y1 <- c(2001, 2003, 2005, 2006)
#' df <- data.table::data.table("m0" = m0, "m1" = m1, "d0" = d0, "d1" = d1, "y0" = y0, "y1" = y1)
#' df <- Date_Shift(df, c("m0", "d0", "y0"), c("m1", "d1", "y1"), "date_since")
#'
Date_Shift <- function(df, dcol0, dcol1, col_name, units = "days") {
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
  def_cols <- names(df)
  df$dt0 <- paste(df[[match(dcol0[1], names(df))]],
    df[[match(dcol0[2], names(df))]],
    df[[match(dcol0[3], names(df))]],
    sep = "-"
  )
  df$dt1 <- paste(df[[match(dcol1[1], names(df))]],
    df[[match(dcol1[2], names(df))]],
    df[[match(dcol1[3], names(df))]],
    sep = "-"
  )
  # TO NOT ENCOUNTER DAYLIGHT SAVINGS ISSUES, THE UTC TIMEZONE IS USED
  # IF NOT USED THEN RESULTS MAY HAVE TIMES OFF BY 1/24 decimals
  df[, col_name] <- difftime(
    strptime(df$dt1,
      format = "%m-%d-%Y",
      tz = "UTC"
    ),
    strptime(df$dt0, format = "%m-%d-%Y"),
    units = units,
    tz = "UTC"
  )
  def_cols <- c(def_cols, col_name)
  return(df[, def_cols, with = FALSE])
}

#' Automates creating a date since a reference column
#'
#' \code{Time_Since} generates a new dataframe with a column containing time since a reference in a given unit
#'
#' @inheritParams R_template
#' @param dcol0 list of ending month, day, and year
#' @family Data Cleaning Functions
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#' m0 <- c(1, 1, 2, 2)
#' m1 <- c(2, 2, 3, 3)
#' d0 <- c(1, 2, 3, 4)
#' d1 <- c(6, 7, 8, 9)
#' y0 <- c(1990, 1991, 1997, 1998)
#' y1 <- c(2001, 2003, 2005, 2006)
#' df <- data.table::data.table(
#'   "m0" = m0, "m1" = m1,
#'   "d0" = d0, "d1" = d1,
#'   "y0" = y0, "y1" = y1
#' )
#' tref <- strptime("3-22-1997", format = "%m-%d-%Y", tz = "UTC")
#' df <- Time_Since(df, c("m1", "d1", "y1"), tref, "date_since")
#'
Time_Since <- function(df, dcol0, tref, col_name, units = "days") {
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
  def_cols <- names(df)
  df$dt0 <- paste(df[[match(dcol0[1], names(df))]], df[[match(
    dcol0[2],
    names(df)
  )]], df[[match(dcol0[3], names(df))]], sep = "-")
  df[, col_name] <- lapply(df$dt0, function(x) {
    (difftime(
      strptime(x,
        format = "%m-%d-%Y", tz = "UTC"
      ), tref,
      units = units
    ))
  })
  def_cols <- c(def_cols, col_name)
  return(df[, def_cols, with = FALSE])
}

#' Automates creating data for a joint competing risks analysis
#'
#' \code{Joint_Multiple_Events} generates input for a regression with multiple non-independent events and models
#'
#' @inheritParams R_template
#' @param events vector of event column names
#' @param term_n_list list of vectors for term numbers for event specific or shared model elements, defaults to term 0
#' @param tform_list list of vectors for subterm types for event specific or shared model elements, defaults to loglinear
#' @param keep_constant_list list of vectors for constant elements for event specific or shared model elements, defaults to free (0)
#' @param a_n_list list of vectors for parameter values for event specific or shared model elements, defaults to term 0
#' @param name_list list of vectors for columns for event specific or shared model elements, required
#' @family Data Cleaning Functions
#' @return returns the updated dataframe and model inputs
#' @export
#' @examples
#' library(data.table)
#' a <- c(0, 0, 0, 1, 1, 1)
#' b <- c(1, 1, 1, 2, 2, 2)
#' c <- c(0, 1, 2, 2, 1, 0)
#' d <- c(1, 1, 0, 0, 1, 1)
#' e <- c(0, 1, 1, 1, 0, 0)
#' df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
#' time1 <- "t0"
#' time2 <- "t1"
#' df$pyr <- df$t1 - df$t0
#' pyr <- "pyr"
#' events <- c("e0", "e1")
#' names_e0 <- c("fac")
#' names_e1 <- c("fac")
#' names_shared <- c("t0", "t0")
#' term_n_e0 <- c(0)
#' term_n_e1 <- c(0)
#' term_n_shared <- c(0, 0)
#' tform_e0 <- c("loglin")
#' tform_e1 <- c("loglin")
#' tform_shared <- c("quad_slope", "loglin_top")
#' keep_constant_e0 <- c(0)
#' keep_constant_e1 <- c(0)
#' keep_constant_shared <- c(0, 0)
#' a_n_e0 <- c(-0.1)
#' a_n_e1 <- c(0.1)
#' a_n_shared <- c(0.001, -0.02)
#' name_list <- list("shared" = names_shared, "e0" = names_e0, "e1" = names_e1)
#' term_n_list <- list("shared" = term_n_shared, "e0" = term_n_e0, "e1" = term_n_e1)
#' tform_list <- list("shared" = tform_shared, "e0" = tform_e0, "e1" = tform_e1)
#' keep_constant_list <- list(
#'   "shared" = keep_constant_shared,
#'   "e0" = keep_constant_e0, "e1" = keep_constant_e1
#' )
#' a_n_list <- list("shared" = a_n_shared, "e0" = a_n_e0, "e1" = a_n_e1)
#' val <- Joint_Multiple_Events(
#'   df, events, name_list, term_n_list,
#'   tform_list, keep_constant_list, a_n_list
#' )
#'
Joint_Multiple_Events <- function(df, events, name_list, term_n_list = list(), tform_list = list(), keep_constant_list = list(), a_n_list = list()) {
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
  # filling missing values
  for (i in names(name_list)) {
    temp0 <- unlist(name_list[i], use.names = FALSE)
    if (i %in% names(term_n_list)) {
      temp1 <- unlist(term_n_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(paste("Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in term_n_list has ",
          length(temp1),
          " items. Omit entry in term_n_list to set to default of term 0 or add missing values",
          sep = ""
        ))
      }
    } else {
      temp <- list(rep(0, length(temp0)))
      names(temp) <- i
      term_n_list <- c(term_n_list, temp)
    }
    if (i %in% names(tform_list)) {
      temp1 <- unlist(tform_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(paste("Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in tform_list has ",
          length(temp1),
          " items. Omit entry in tform_list to set to default of 'loglin' or add missing values",
          sep = ""
        ))
      }
    } else {
      temp <- list(rep("loglin", length(temp0)))
      names(temp) <- i
      tform_list <- c(tform_list, temp)
    }
    if (i %in% names(keep_constant_list)) {
      temp1 <- unlist(keep_constant_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(paste("Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in keep_constant_list has ",
          length(temp1),
          " items. Omit entry in tform_list to set to default of 0 or add missing values",
          sep = ""
        ))
      }
    } else {
      temp <- list(rep(0, length(temp0)))
      names(temp) <- i
      keep_constant_list <- c(keep_constant_list, temp)
    }
    if (i %in% names(a_n_list)) {
      temp1 <- unlist(a_n_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(paste("Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in a_n_list has ",
          length(temp1),
          " items. Omit entry in a_n_list to set to default of 0 or add missing values",
          sep = ""
        ))
      }
    } else {
      temp <- list(rep(0, length(temp0)))
      names(temp) <- i
      a_n_list <- c(a_n_list, temp)
    }
  }
  df0 <- data.table()
  for (i in names(df)) {
    if (i %in% events) {
      if (i == events[1]) {
        temp <- c()
        for (j in events) {
          temp <- c(temp, unlist(df[, j, with = FALSE], use.names = FALSE))
        }
        df0[, "events"] <- temp
      }
      temp <- c()
      for (j in events) {
        if (i == j) {
          temp <- c(temp, rep(1, nrow(df)))
        } else {
          temp <- c(temp, rep(0, nrow(df)))
        }
      }
      df0[, i] <- temp
    } else {
      temp <- rep(unlist(df[, i, with = FALSE], use.names = FALSE), length(events))
      df0[, i] <- temp
    }
  }
  names <- c()
  term_n <- c()
  tform <- c()
  keep_constant <- c()
  a_n <- c()
  if ("shared" %in% names(name_list)) {
    names <- c(names, unlist(name_list["shared"], use.names = FALSE))
    term_n <- c(term_n, unlist(term_n_list["shared"], use.names = FALSE))
    tform <- c(tform, unlist(tform_list["shared"], use.names = FALSE))
    keep_constant <- c(keep_constant, unlist(keep_constant_list["shared"],
      use.names = FALSE
    ))
    a_n <- c(a_n, unlist(a_n_list["shared"], use.names = FALSE))
  }
  for (i in events) {
    if (i %in% names(name_list)) {
      interactions <- c()
      new_names <- c()
      for (j in unlist(name_list[i], use.names = FALSE)) {
        interactions <- c(interactions, paste(j, "?*?", i, sep = ""))
        new_names <- c(new_names, paste(j, "_", i, sep = ""))
      }
      vals <- interact_them(df0, interactions, new_names)
      df0 <- vals$df
      new_names <- vals$cols
      names <- c(names, new_names)
      term_n <- c(term_n, unlist(term_n_list[i], use.names = FALSE))
      tform <- c(tform, unlist(tform_list[i], use.names = FALSE))
      keep_constant <- c(keep_constant, unlist(keep_constant_list[i],
        use.names = FALSE
      ))
      a_n <- c(a_n, unlist(a_n_list[i], use.names = FALSE))
    }
  }
  return(list("df" = df0, "names" = names, "term_n" = term_n, "tform" = tform, "keep_constant" = keep_constant, "a_n" = a_n))
}

#' Defines Interactions
#'
#' \code{interact_them} uses user provided interactions define interaction terms and update the data.table. assumes interaction is "+" or "*" and applies basic anti-aliasing to avoid duplicates
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
interact_them <- function(df, interactions, new_names, verbose = 0) {
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
  verbose <- Check_Verbose(verbose)
  cols <- c()
  for (i in seq_len(length(interactions))) {
    interac <- interactions[i]
    formula <- unlist(strsplit(interac, "[?]"), use.names = FALSE)
    if (length(formula) != 3) {
      stop(paste(
        "Error: Iteration:", interac, "has incorrect length of",
        length(formula), "but should be 3."
      ))
    }
    newcol <- paste(formula[1], formula[2], formula[3], sep = "")
    if (new_names[i] != "") {
      newcol <- new_names[i]
    }
    col1 <- formula[1]
    col2 <- formula[3]
    if (paste(formula[1], "?", formula[2], "?", formula[3], sep = "") %in% interactions[i + seq_len(length(interactions))]) {
      if (verbose >= 2) {
        warning(paste("Warning: interation ", i, "is duplicated")) # nocov
      }
    } else if (paste(formula[3], "?", formula[2], "?", formula[1], sep = "") %in% interactions[i + seq_len(length(interactions))]) {
      if (verbose >= 2) {
        warning(paste(
          "Warning: the reverse of interation ", i,
          "is duplicated"
        )) # nocov
      }
    } else {
      if (formula[2] == "+") {
        df[, newcol] <- df[, col1, with = FALSE] +
          df[, col2, with = FALSE]
        cols <- c(cols, newcol)
      } else if (formula[2] == "*") {
        df[, newcol] <- df[, col1, with = FALSE] *
          df[, col2, with = FALSE]
        cols <- c(cols, newcol)
      } else {
        stop(paste("Error: Incorrect operation of", formula[2]))
      }
    }
  }
  list("df" = df, "cols" = cols)
}

#' Automatically applies a normalization to either an input or output
#'
#' \code{apply_norm} applies a normalization factor
#'
#' @param df The data.table with columns to be normalized
#' @param norm The normalization option used, currently max or mean
#' @param input boolean if the normalization is being performed on the input values or on an output
#' @param values list of values using during normalization
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns list with the normalized values
apply_norm <- function(df, norm, names, input, values, model_control) {
  if (input) {
    a_n <- values$a_n
    cons_mat <- values$cons_mat
    if (tolower(norm) == "null") {
      # nothing changes
      norm_weight <- rep(1, length(names))
    } else if (tolower(norm) %in% c("max", "mean")) {
      # weight by the maximum value
      norm_weight <- c()
      if (tolower(norm) == "max") {
        for (i in seq_along(names)) {
          val <- summarise(df, max_value = max(abs(get(names[i]))))[[1]]
          if (val == 0.0) {
            warning(paste("Warning: Maximum value for ", names[i], " was 0. Normalization not applied to column.", sep = ""))
            val <- 1.0
          }
          norm_weight <- c(norm_weight, val)
        }
      } else if (tolower(norm) == "mean") {
        for (i in seq_along(names)) {
          val <- summarise(df, mean_value = mean(get(names[i])))[[1]]
          if (val == 0.0) {
            warning(paste("Warning: Average value for ", names[i], " was 0. Normalization not applied to column.", sep = ""))
            val <- 1.0
          }
          norm_weight <- c(norm_weight, val)
        }
      } else {
        stop(paste("Error: norm option ", norm, " wasn't coded yet", sep = ""))
      }
      for (i in seq_along(names)) {
        if (typeof(a_n) != "list") {
          a_n[i] <- a_n[i] * norm_weight[i]
        } else {
          for (j in seq_len(length(a_n))) {
            a_n[[j]][i] <- a_n[[j]][i] * norm_weight[i]
          }
        }
        if (match(names[i], names)[1] == i) {
          df[, names[i]] <- df[, names[i], with = FALSE] / norm_weight[i]
        }
        norm_weight <- c(norm_weight, val)
      }
      if (model_control[["constraint"]] == TRUE) {
        for (i in seq_along(names)) {
          cons_mat[, i] <- cons_mat[, i] / norm_weight[i]
        }
      }
    } else {
      stop(gettextf(
        "Error: Normalization arguement '%s' not valid.",
        norm
      ), domain = NA)
    }
    output <- list(
      "a_n" = a_n,
      "cons_mat" = cons_mat,
      "norm_weight" = norm_weight,
      "df" = df
    )
  } else {
    res <- values$output
    norm_weight <- values$norm_weight
    keep_constant <- res$Parameter_Lists$keep_constant
    if (tolower(norm) == "null") {
      # nothing changes
    } else if (tolower(norm) %in% c("mean", "max")) {
      # weight by the maximum value
      if (model_control$single) {
        for (i in seq_along(names)) {
          res$beta_0[i] <- res$beta_0[i] / norm_weight[i]
        }
      } else {
        for (i in seq_along(names)) {
          if (keep_constant[i] == 0) {
            i_der <- i - sum(head(keep_constant, i))
            res$First_Der[i_der] <- res$First_Der[i_der] * norm_weight[i]
            res$beta_0[i] <- res$beta_0[i] / norm_weight[i]
            res$Standard_Deviation[i] <- res$Standard_Deviation[i] / norm_weight[i]
            for (j in seq_along(names)) {
              if (keep_constant[j] == 0) {
                j_der <- j - sum(head(keep_constant, j))
                res$Second_Der[i_der, j_der] <- res$Second_Der[i_der, j_der] * norm_weight[i] * norm_weight[j]
                res$Covariance[i_der, j_der] <- res$Covariance[i_der, j_der] / norm_weight[i] / norm_weight[j]
              }
            }
          }
        }
        if (model_control[["constraint"]] == TRUE) {
          for (i in seq_along(names)) {
            res$constraint_matrix[, i] <- res$constraint_matrix[, i] * norm_weight[i]
          }
        }
      }
    } else {
      stop(gettextf(
        "Error: Normalization arguement '%s' not valid.",
        norm
      ), domain = NA)
    }
    output <- res
  }
  output
}

#' Checks system OS
#'
#' \code{get_os} checks the system OS, part of configuration script
#'
#' @noRd
#' @return returns a string representation of OS
get_os <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin") {
      os <- "osx" # nocov
    }
  } else { ## mystery machine
    os <- .Platform$OS.type # nocov
    if (grepl("^darwin", R.version$os)) { # nocov
      os <- "osx" # nocov
    }
    if (grepl("linux-gnu", R.version$os)) { # nocov
      os <- "linux" # nocov
    }
  }
  tolower(os)
}

#' Checks default c++ compiler
#'
#' \code{gcc_version} Checks default c++ compiler, part of configuration script
#'
#' @noRd
#' @return returns a string representation of gcc, clang, or c++ output
gcc_version <- function() {
  #  tstart <- Sys.time()
  out <- tryCatch(run("c++", "-v"),
    error = function(cnd) list(stdout = "")
  )
  #  tend <- Sys.time()
  #  print("out call")
  #  print(tend - tstart)
  out0 <- str_match(out$stdout, "gcc version")[1]
  if (!is.na(out0)) {
    out <- "gcc"
  } else {
    out0 <- str_match(out$stdout, "clang version")[1] # nocov
    if (!is.na(out0)) { # nocov
      out <- "clang" # nocov
    } else {
      out <- out$stdout # nocov
    }
  }
  out
}

#' Checks how R was compiled
#'
#' \code{Rcomp_version} Checks how R was compiled, part of configuration script
#'
#' @noRd
#' @return returns a string representation of gcc, clang, or R CMD config CC output
Rcomp_version <- function() {
  out <- rcmd("config", "CC")
  out0 <- str_match(out$stdout, "clang")[1]
  if (!is.na(out0)) {
    out <- "clang" # nocov
  } else {
    out0 <- str_match(out$stdout, "gcc")[1]
    if (!is.na(out0)) {
      out <- "gcc"
    } else {
      out <- out$stdout # nocov
    }
  }
  out
}

#' Checks default R c++ compiler
#'
#' \code{Rcpp_version} checks ~/.R/Makevars script for default compilers set, part of configuration script
#'
#' @noRd
#' @return returns a string representation of gcc, clang, or head ~/.R/Makevars
Rcpp_version <- function() {
  out <- tryCatch(run("head", "~/.R/Makevars", stderr_to_stdout = TRUE),
    error = function(cnd) list(stdout = "")
  )
  out0 <- str_match(out$stdout, "clang")[1]
  if (!is.na(out0)) {
    out <- "clang" # nocov
  } else {
    out0 <- str_match(out$stdout, "gcc")[1]
    if (!is.na(out0)) {
      out <- "gcc" # nocov
    } else {
      out <- out$stdout # nocov
    }
  }
  out
}

#' Checks OS, compilers, and OMP
#'
#' \code{System_Version} checks OS, default R c++ compiler, and if OMP is enabled
#'
#' @return returns a list of results
#' @export
#' @family Output and Information Functions
System_Version <- function() {
  tstart <- Sys.time()
  os <- get_os()
  gcc <- gcc_version()
  Rcomp <- Rcomp_version()
  OMP <- OMP_Check()
  list(
    "Operating System" = os, "Default c++" = gcc, "R Compiler" = Rcomp,
    "OpenMP Enabled" = OMP
  )
}

#' General purpose verbosity check
#'
#' \code{Check_Verbose} checks and assigns verbosity values
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns correct verbose value
#'
Check_Verbose <- function(verbose) {
  if (verbose %in% c(0, 1, 2, 3, 4)) {
    verbose <- as.integer(verbose)
  } else if (verbose %in% c(TRUE, FALSE)) {
    if (verbose) {
      verbose <- 3 # nocov
    } else {
      verbose <- 0 # nocov
    }
  } else {
    stop("Error: verbosity argument not valid")
  }
  verbose
}

#' uses a table, list of categories, and list of event summaries to generate person-count tables
#'
#' \code{Event_Count_Gen} generates event-count tables
#'
#' @inheritParams R_template
#' @param table dataframe with every category/event column needed
#' @param categ list with category columns and methods, methods can be either strings or lists of boundaries
#' @param events list of columns to summarize, supports counts and means and renaming the summary column
#'
#' @return returns a grouped table and a list of category boundaries used
#' @family Data Cleaning Functions
#' @export
#' @examples
#' library(data.table)
#' a <- c(0, 1, 2, 3, 4, 5, 6)
#' b <- c(1, 2, 3, 4, 5, 6, 7)
#' c <- c(0, 1, 0, 0, 0, 1, 0)
#' table <- data.table::data.table(
#'   "a" = a,
#'   "b" = b,
#'   "c" = c
#' )
#' categ <- list(
#'   "a" = "0/3/5]7",
#'   "b" = list(
#'     lower = c(-1, 3, 6),
#'     upper = c(3, 6, 10),
#'     name = c("low", "medium", "high")
#'   )
#' )
#' event <- list(
#'   "c" = "count AS cases",
#'   "a" = "mean", "b" = "mean"
#' )
#' e <- Event_Count_Gen(table, categ, event, T)
#'
Event_Count_Gen <- function(table, categ, events, verbose = FALSE) {
  df <- as_tibble(table)
  `%>%` <- dplyr::`%>%`
  #
  categ_cols <- c()
  categ_bounds <- list()
  for (cat in names(categ)) {
    cat_str <- ""
    if (!cat %in% names(table)) {
      stop(paste("Error: ", cat, " not in table", sep = ""))
    }
    if (length(categ[[cat]]) > 1) { # list of bounds
      temp0 <- categ[[cat]]$lower
      temp1 <- categ[[cat]]$upper
      if ("name" %in% names(categ[[cat]])) { # assign names to the levels
        temp2 <- categ[[cat]]$name
      } else { # number the categories
        temp2 <- seq_len(length(temp0)) # 1:length(temp0)
      }
      num_categ <- length(temp0) # number of categories
      cat_col <- paste(cat, "category", sep = "_") # name of the category
      categ_cols <- c(categ_cols, cat_col) # add to list
      df <- df %>% mutate("{cat_col}" := "Unassigned") # add to tibble
      for (i in 1:num_categ) { # for each category
        L <- as.numeric(temp0[i]) # lower bound
        if (grepl("]", temp1[i], fixed = TRUE)) {
          U <- as.numeric(gsub("]", "", temp1[i])) # assign upper
          a_col_categ <- case_when(df[[cat]] <= U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]]) # rows within the bounds and unassigned are assigned the name
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ") # add to list of intervals
        } else {
          U <- as.numeric(temp1[i]) # assign upper
          if (L == U) { # discrete case
            a_col_categ <- case_when(df[[cat]] == U & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          } else { # interval case
            a <- case_when(df[[cat]] < U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
          }
        }
        df[[cat_col]] <- a_col_categ # update column
      }
      categ_bounds[[cat_col]] <- cat_str # add to list of intervals
    } else { # string of bounds
      cat_str <- ""
      temp <- categ[[cat]]
      temp <- gsub("/", " / ", temp)
      temp <- gsub("]", " ] ", temp)
      temp <- strsplit(temp, "\\s+")[[1]] # adding and splitting by spaces
      num_categ <- (length(temp) - 1) / 2 # get number of categories
      cat_col <- paste(cat, "category", sep = "_")
      categ_cols <- c(categ_cols, cat_col)
      df <- df %>% mutate("{cat_col}" := "Unassigned") # add to tibble
      for (i in 1:num_categ) {
        L <- as.numeric(temp[2 * i - 1]) # get upper and lower bounds
        U <- as.numeric(temp[2 * i + 1])
        if (L == U) { # dicrete value case
          a_col_categ <- case_when(df[[cat]] == U & df[[cat_col]] == "Unassigned" ~ as.character(i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
        } else if (temp[2 * i] == "/") { # strict less than upper bound
          a_col_categ <- case_when(df[[cat]] < U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
        } else { # within and including bounds
          a_col_categ <- case_when(df[[cat]] <= U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
        }
        df[[cat_col]] <- a_col_categ # update column
      }
      categ_bounds[[cat_col]] <- cat_str # add interval boundaries
    }
  }
  df_group <- df %>%
    group_by(across(all_of(categ_cols))) %>%
    summarize("COUNT" = n(), .groups = "drop") # group by columns and summarize by counts
  for (evt in names(events)) { # iterate through events
    if (grepl(" AS ", events[[evt]], fixed = TRUE)) { # get method and updated name
      temp <- gsub(" AS ", " ", events[[evt]])
      temp <- strsplit(temp, "\\s+")[[1]]
      col_name <- temp[2]
      method <- temp[1]
    } else {
      col_name <- evt
      temp <- strsplit(events[[evt]], "\\s+")[[1]]
      method <- temp[[length(temp)]]
    }
    if (method == "count") { # summarize by count
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method == "mean") { # summarize by mean
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := mean(.data[[evt]]), .groups = "drop")
    }
    df_group[[col_name]] <- df_temp[[col_name]] # add to grouped tibble
  }
  # return the grouped list and list of catgegory bounds
  list(df = as.data.table(df_group), bounds = categ_bounds)
}

#' uses a table, list of categories, list of summaries, list of events, and person-year information to generate person-time tables
#'
#' \code{Event_Time_Gen} generates event-time tables
#'
#' @inheritParams R_template
#' @param table dataframe with every category/event column needed
#' @param pyr list with entry and exit lists, containing day/month/year columns in the table
#' @param categ list with category columns and methods, methods can be either strings or lists of boundaries, includes a time category or entry/exit are both required for the pyr list
#' @param summaries list of columns to summarize, supports counts, means, and weighted means by person-year and renaming the summary column
#' @param events list of events or interests, checks if events are within each time interval
#'
#' @return returns a grouped table and a list of category boundaries used
#' @family Data Cleaning Functions
#' @export
#' @examples
#' library(data.table)
#' a <- c(0, 1, 2, 3, 4, 5, 6)
#' b <- c(1, 2, 3, 4, 5, 6, 7)
#' c <- c(0, 1, 0, 0, 0, 1, 0)
#' d <- c(1, 2, 3, 4, 5, 6, 7)
#' e <- c(2, 3, 4, 5, 6, 7, 8)
#' f <- c(
#'   1900, 1900, 1900, 1900,
#'   1900, 1900, 1900
#' )
#' g <- c(1, 2, 3, 4, 5, 6, 7)
#' h <- c(2, 3, 4, 5, 6, 7, 8)
#' i <- c(
#'   1901, 1902, 1903, 1904,
#'   1905, 1906, 1907
#' )
#' table <- data.table::data.table(
#'   "a" = a, "b" = b, "c" = c,
#'   "d" = d, "e" = e, "f" = f,
#'   "g" = g, "h" = h, "i" = i
#' )
#' categ <- list(
#'   "a" = "-1/3/5]7",
#'   "b" = list(
#'     lower = c(-1, 3, 6), upper = c(3, 6, 10),
#'     name = c("low", "medium", "high")
#'   ),
#'   "time AS time" = list(
#'     "day" = c(1, 1, 1, 1, 1),
#'     "month" = c(1, 1, 1, 1, 1),
#'     "year" = c(1899, 1903, 1910)
#'   )
#' )
#' summary <- list(
#'   "c" = "count AS cases",
#'   "a" = "mean",
#'   "b" = "weighted_mean"
#' )
#' events <- list("c")
#' pyr <- list(
#'   entry = list(year = "f", month = "e", day = "d"),
#'   exit = list(year = "i", month = "h", day = "g"),
#'   unit = "years"
#' )
#' e <- Event_Time_Gen(table, pyr, categ, summary, events, T)
#'
Event_Time_Gen <- function(table, pyr, categ, summaries, events, verbose = FALSE) {
  df <- as_tibble(table)
  `%within%` <- lubridate::`%within%`
  `%>%` <- dplyr::`%>%`
  #
  categ_cols <- c()
  categ_bounds <- list()
  for (cat in names(categ)) { # for each category
    cat_str <- ""
    if (grepl(" AS ", cat, fixed = TRUE)) { # get the column and name
      temp <- strsplit(gsub(" AS ", " ", cat), "\\s+")[[1]]
      cat_col <- temp[2]
      cat_df <- temp[1]
    } else {
      cat_df <- cat
      cat_col <- paste(cat, "category", sep = "_")
    }
    if (cat_col %in% names(df)) { # check that the category doesn't already exist in the original dataframe
      stop(paste("Error: ", cat_col, " already exists, use ' AS ' to rename if needed", sep = ""))
    }
    if (length(categ[[cat]]) > 1) { # boundary as lists
      if ("lower" %in% names(categ[[cat]])) { # lower and upper boundary intervals
        if (!cat_df %in% names(table)) {
          stop(paste("Error: ", cat_df, " not in table", sep = ""))
        }
        temp0 <- categ[[cat]]$lower
        temp1 <- categ[[cat]]$upper
        if ("name" %in% names(categ[[cat]])) { # check for names for each level
          temp2 <- categ[[cat]]$name
        } else {
          temp2 <- seq_len(length(temp0)) # 1:length(temp0)
        }
        num_categ <- length(temp0)
        categ_cols <- c(categ_cols, cat_col)
        df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
        for (i in 1:num_categ) { # for each level
          L <- as.numeric(temp0[i])
          if (grepl("]", temp1[i], fixed = TRUE)) { # check for including the upper limit
            U <- as.numeric(gsub("]", "", temp1[i])) # get upper limit
            a_col_categ <- case_when(df[[cat_df]] <= U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]]) # assign the level to unassigned rows
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ") # add boundary information to list of intervals
          } else {
            U <- as.numeric(temp1[i]) # get upper limit
            if (L == U) { # discrete case
              a_col_categ <- case_when(df[[cat_df]] == U & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
              cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
            } else { # interval case
              a_col_categ <- case_when(df[[cat_df]] < U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
              cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
            }
          }
          df[[cat_col]] <- a_col_categ # update tibble
        }
      } else { # calender time scale
        if ("day" %in% names(categ[[cat]])) { # determine the day, month, year data for the category
          day <- categ[[cat]]$day # nocov
          if ("month" %in% names(categ[[cat]])) { # nocov
            month <- categ[[cat]]$month # nocov
          } else { # nocov
            month <- rep(1, length(day)) # nocov
          }
          if ("year" %in% names(categ[[cat]])) { # nocov
            year <- categ[[cat]]$year # nocov
          } else {
            year <- rep(1900, length(day)) # nocov
          }
        } else if ("month" %in% names(categ[[cat]])) { # nocov
          month <- categ[[cat]]$month # nocov
          day <- rep(1, length(month)) # nocov
          if ("year" %in% names(categ[[cat]])) { # nocov
            year <- categ[[cat]]$year # nocov
          } else { # nocov
            year <- rep(1900, length(month)) # nocov
          }
        } else if ("year" %in% names(categ[[cat]])) { # nocov
          year <- categ[[cat]]$year # nocov
          day <- rep(1, length(year)) # nocov
          month <- rep(1, length(month)) # nocov
        } else {
          stop("Error: calender category missing 'day', 'month', and 'year'") # all three are required
        }
        num_categ <- length(year) - 1
        categ_cols <- c(categ_cols, cat_col)
        df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
        if ("PYR" %in% names(df)) {
          stop("Error: either multiple time categorizations used or 'PYR' column already exists")
        }
        df <- df %>% mutate("PYR" := 0)
        interval <- "unassigned"
        if ("exit" %in% names(pyr)) { # format the exit as date column
          pyr_exit <- pyr$exit
          if ("year" %in% names(pyr_exit)) { # nocov
            if ("month" %in% names(pyr_exit)) { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(year = table[[pyr_exit$year]], month = table[[pyr_exit$month]], day = table[[pyr_exit$day]]) # nocov
              } else { # nocov
                exit <- make_date(year = table[[pyr_exit$year]], month = table[[pyr_exit$month]]) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(year = table[[pyr_exit$year]], day = table[[pyr_exit$day]]) # nocov
              } else { # nocov
                exit <- make_date(year = table[[pyr_exit$year]]) # nocov
              }
            }
          } else {
            if ("month" %in% names(pyr_exit)) { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(month = table[[pyr_exit$month]], day = table[[pyr_exit$day]]) # nocov
              } else { # nocov
                exit <- make_date(month = table[[pyr_exit$month]]) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(day = table[[pyr_exit$day]]) # nocov
              } else {
                stop("Error: person-year exit missing day, month, and year")
              }
            }
          }
          if ("entry" %in% names(pyr)) { # format the entry as date column
            pyr_entry <- pyr$entry
            interval <- "interval"
            if ("year" %in% names(pyr_entry)) { # nocov
              if ("month" %in% names(pyr_entry)) { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
                } else { # nocov
                  entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]]) # nocov
                }
              } else { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(year = table[[pyr_entry$year]], day = table[[pyr_entry$day]]) # nocov
                } else { # nocov
                  entry <- make_date(year = table[[pyr_entry$year]]) # nocov
                }
              }
            } else { # nocov
              if ("month" %in% names(pyr_entry)) { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
                } else { # nocov
                  entry <- make_date(month = table[[pyr_entry$month]]) # nocov
                }
              } else { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(day = table[[pyr_entry$day]]) # nocov
                } else {
                  stop("Error: person-year entry missing day, month, and year")
                }
              }
            }
          } else {
            pyr_entry <- list()
            interval <- "left trunc" # if there is an exit and no entry then the data is left truncated
          }
        } else if ("entry" %in% names(pyr)) {
          pyr_entry <- pyr$entry # format the entry as date column
          if ("year" %in% names(pyr_entry)) { # nocov
            if ("month" %in% names(pyr_entry)) { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
              } else { # nocov
                entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]]) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(year = table[[pyr_entry$year]], day = table[[pyr_entry$day]]) # nocov
              } else { # nocov
                entry <- make_date(year = table[[pyr_entry$year]]) # nocov
              }
            }
          } else { # nocov
            if ("month" %in% names(pyr_entry)) { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
              } else { # nocov
                entry <- make_date(month = table[[pyr_entry$month]]) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(day = table[[pyr_entry$day]]) # nocov
              } else {
                stop("Error: person-year entry missing day, month, and year")
              }
            }
          }
          pyr_exit <- list()
          interval <- "right trunc" # entry and no exit is right truncated
        } else {
          stop("Error: Date columns given for category, but person-year data not in date format") # pyr does not have exit or entry values
        }
        df_added <- tibble()
        pyr_unit <- "years" # default person-years to years
        if ("unit" %in% names(pyr)) {
          pyr_unit <- pyr$unit
        }
        for (time_i in 1:num_categ) { # for every time interval
          istart <- make_date(year = year[time_i], month = month[time_i], day = day[time_i]) # interval start
          iend <- make_date(year = year[time_i + 1], month = month[time_i + 1], day = day[time_i + 1]) # interval end
          cat_str <- paste(cat_str, paste("[", istart, " to ", iend, "]", sep = ""), sep = " ") # prepare the interval info
          categ_interval <- interval(istart, iend) # define as date interval
          c_categ <- list()
          if (interval == "left trunc") {
            # only exit
            a_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # exit within interval, set to i
            for (evt in events) {
              c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # exit within interval set to event column value
            }
            b_categ <- case_when(a_categ == as.character(time_i) ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = 0) # for every row with exit in interval, track the duration from interval start
            risk_interval <- interval(iend, exit) # interval from interval end to row end
            a_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # track rows which end after the interval, set to i
            b_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(categ_interval), pyr_unit), .default = b_categ) # rows which end after interval are at risk the full interval
          } else if (interval == "right trunc") {
            # only entry
            a_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") #  start during the interval, set to i
            for (evt in events) {
              c_categ[[evt]] <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # start during interval, set to event value
            }
            b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(interval(entry, iend)), pyr_unit), .default = 0) # every row which starts during the interval, tracks duration from entry to interval end
            risk_interval <- interval(entry, istart)
            a_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # tracks which rows start before the interval
            b_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(categ_interval), pyr_unit), .default = b_categ) # rows which start before the interval are at risk the full interval
          } else {
            # both entry and exit
            risk_interval <- interval(entry, exit)
            a_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # interval fully contained, set to i
            a_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval ends during row interval, set to i
            a_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval starts during row interval, set to i
            a_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # row interval fully contained in category interval, set to i
            for (evt in events) {
              c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # row ends during category interval, set to event value
            }
            b_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = 0) # rows which end during the category interval, track category interval start to row end
            b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(interval(entry, iend)), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
            b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(categ_interval), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
            b_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(risk_interval), pyr_unit), .default = b_categ) # row interval fully in category interval, track full row interval
          }
          index_kept <- seq_len(nrow(df))
          index_kept <- index_kept[a_categ == time_i] # indexes which contain the category interval to some level
          b_categ <- b_categ[a_categ == time_i] # durations for kept indexes
          row_kept <- slice(df, index_kept) # dataframe at kept rows
          row_kept <- row_kept %>% mutate("{cat_col}" := as.character(time_i)) # time category value
          row_kept <- row_kept %>% mutate("PYR" := b_categ) # duration value
          for (evt in events) {
            d_categ <- c_categ[[evt]]
            d_categ <- d_categ[a_categ == time_i]
            row_kept <- row_kept %>% mutate("{evt}" := d_categ) # event values
          }
          df_added <- bind_rows(df_added, row_kept) # new updates dataset
        }
        df <- df_added
      }
      categ_bounds[[cat_col]] <- cat_str # update the list of category boundaries
    } else { # boundary as string, not a time category
      if (!cat_df %in% names(table)) {
        stop(paste("Error: ", cat_df, " not in table", sep = ""))
      }
      cat_str <- ""
      temp <- categ[[cat]]
      temp <- gsub("/", " / ", temp)
      temp <- gsub("]", " ] ", temp)
      temp <- strsplit(temp, "\\s+")[[1]] # seperate values and delimiters
      num_categ <- (length(temp) - 1) / 2
      categ_cols <- c(categ_cols, cat_col)
      df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize column
      for (time_i in 1:num_categ) {
        L <- as.numeric(temp[2 * time_i - 1])
        U <- as.numeric(temp[2 * time_i + 1]) # get lower and upper bounds
        if (L == U) { # discrete case
          a_categ <- case_when(df[[cat_df]] == U & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
        } else if (temp[2 * time_i] == "/") { # strictly below upper bound
          a_categ <- case_when(df[[cat_df]] < U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
        } else { # including both bounds
          a_categ <- case_when(df[[cat_df]] <= U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = df[[cat_col]])
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
        }
        df[[cat_col]] <- a_categ # update tibble
      }
      categ_bounds[[cat_col]] <- cat_str
    }
  }
  if ("PYR" %in% names(df)) {
    # good, we need a person-years measure
  } else {
    if ("exit" %in% names(pyr)) { # format the exit as date column
      pyr_exit <- pyr$exit
      if ("year" %in% names(pyr_exit)) { # nocov
        if ("month" %in% names(pyr_exit)) { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(year = table[[pyr_exit$year]], month = table[[pyr_exit$month]], day = table[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(year = table[[pyr_exit$year]], month = table[[pyr_exit$month]]) # nocov
          }
        } else { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(year = table[[pyr_exit$year]], day = table[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(year = table[[pyr_exit$year]]) # nocov
          }
        }
      } else {
        if ("month" %in% names(pyr_exit)) { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(month = table[[pyr_exit$month]], day = table[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(month = table[[pyr_exit$month]]) # nocov
          }
        } else { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(day = table[[pyr_exit$day]]) # nocov
          } else {
            stop("Error: person-year exit missing day, month, and year")
          }
        }
      }
      if ("entry" %in% names(pyr)) { # format the entry as date column
        pyr_entry <- pyr$entry
        interval <- "interval"
        if ("year" %in% names(pyr_entry)) { # nocov
          if ("month" %in% names(pyr_entry)) { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(year = table[[pyr_entry$year]], month = table[[pyr_entry$month]]) # nocov
            }
          } else { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(year = table[[pyr_entry$year]], day = table[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(year = table[[pyr_entry$year]]) # nocov
            }
          }
        } else { # nocov
          if ("month" %in% names(pyr_entry)) { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(month = table[[pyr_entry$month]], day = table[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(month = table[[pyr_entry$month]]) # nocov
            }
          } else { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(day = table[[pyr_entry$day]]) # nocov
            } else {
              stop("Error: person-year entry missing day, month, and year")
            }
          }
        }
      } else {
        stop("Error: Both entry and exit needed to calculate person-years, without a time category for reference")
      }
      risk_interval <- interval(entry, exit)
      pyr_unit <- "years" # default person-years to years
      if ("unit" %in% names(pyr)) {
        pyr_unit <- pyr$unit
      }
      df <- df %>% mutate("PYR" := as.numeric(as.duration(risk_interval), pyr_unit))
    } else {
      stop("Error: Both entry and exit needed to calculate person-years, without a time category for reference")
    }
  }
  # check that required columns are not already in use
  if ("PYR" %in% names(summaries)) { # storing person-year durations
    stop("Error: 'PYR' listed as a event column, either remove or rename with ' AS '")
  }
  if ("AT_RISK" %in% names(summaries)) { # storing number at risk
    stop("Error: 'AT_RISK' listed as a event column, either remove or rename with ' AS '")
  }
  df_group <- df %>%
    group_by(across(all_of(categ_cols))) %>%
    summarize("AT_RISK" = n(), "PYR" = sum(.data[["PYR"]]), .groups = "drop") # group by categories and define the durations and counts
  for (evt in names(summaries)) { # for each event summary
    if (grepl(" AS ", summaries[[evt]], fixed = TRUE)) { # get the method and column name
      temp <- gsub(" AS ", " ", summaries[[evt]])
      temp <- strsplit(temp, "\\s+")[[1]]
      col_name <- temp[2]
      method <- temp[1]
    } else {
      col_name <- evt
      temp <- strsplit(summaries[[evt]], "\\s+")[[1]]
      method <- temp[[length(temp)]]
    }
    if (method == "count") { # sum of event across each category combination
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method == "mean") { # mean value of event across each category combination
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := mean(.data[[evt]]), .groups = "drop")
    } else if (method == "weighted_mean") { # mean value weighted by person-years
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := weighted.mean(.data[[evt]], .data[["PYR"]]), .groups = "drop")
    }
    df_group[[col_name]] <- df_temp[[col_name]]
  }
  #
  list(df = as.data.table(df_group), bounds = categ_bounds)
}

#' Prints a cox regression output clearly
#'
#' \code{print.coxres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.coxres <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a poisson regression output clearly
#'
#' \code{print.poisres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.poisres <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a case-control regression output clearly
#'
#' \code{print.caseconres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class caseconres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.caseconres <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a logistic regression output clearly
#'
#' \code{print.logitres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class logitres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.logitres <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a cox likelihood boundary regression output clearly
#'
#' \code{print.coxresbound} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxresbound
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.coxresbound <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a poisson likelihood boundary regression output clearly
#'
#' \code{print.poisresbound} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisresbound
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @export
#' @family Output and Information Functions
print.poisresbound <- function(x, ...) {
  exargs <- list(...)
  digits <- 2
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_Output(x, digits)
}

#' Prints a regression output clearly
#'
#' \code{Interpret_Output} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return return nothing, prints the results to console
Interpret_Output <- function(out_list, digits = 2) {
  # make sure the output isn't an error
  passed <- out_list$Status
  message("|-------------------------------------------------------------------|")
  if (!is.na(passed)) {
    if ("Likelihood_Goal" %in% names(out_list)) {
      # likelihood boundary output
      model <- out_list$model
      modelcontrol <- out_list$modelcontrol
      para_number <- modelcontrol$para_number
      #
      name <- model$names[para_number]
      tform <- model$tform[para_number]
      term_n <- model$term_n[para_number]
      beta_0 <- out_list$beta_0[para_number]
      #
      limits <- out_list$Parameter_Limits
      neg <- out_list$Negative_Risk_Limit_Hit
      conv <- out_list$Limit_Converged
      lik_bound <- out_list$Likelihood_Boundary
      lik_goal <- out_list$Likelihood_Goal
      message("Likelihood Boundary Results")
      message(paste("Solving for the boundary of element: ", para_number, "\nApplied to column: '", name, "'\nSubterm: ", tform, "\nTerm number: ", term_n, sep = ""))
      if (neg[1]) {
        message("Lower limit was not found")
      } else {
        if (conv[1]) {
          message(paste("Lower limit converged to at ", round(limits[1], digits), " at a score of ", round(lik_bound[1], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
        } else {
          message(paste("Lower limit reached ", round(limits[1], digits), " at a score of ", round(lik_bound[1], digits), " with of goal of ", round(lik_goal, digits), " but did not converge", sep = ""))
        }
      }
      message(paste("Central estimate was ", round(beta_0, digits), sep = ""))
      if (neg[2]) {
        message("Upper limit was not found")
      } else {
        if (conv[2]) {
          message(paste("Upper limit converged to at ", round(limits[2], digits), " at a score of ", round(lik_bound[2], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
        } else {
          message(paste("Upper limit reached ", round(limits[2], digits), " at a score of ", round(lik_bound[2], digits), " with of goal of ", round(lik_goal, digits), " but did not converge", sep = ""))
        }
      }
    } else {
      # Check if its a multidose problem
      if (out_list$Survival_Type == "Cox_Multidose") {
        message("Currently the multiple realization code is not setup for printing results, due to the potentially large number of realizations")
      } else if (out_list$Survival_Type == "CaseControl") {
        # case control output
        # get the model details
        names <- out_list$Parameter_Lists$names
        tforms <- out_list$Parameter_Lists$tforms
        term_n <- out_list$Parameter_Lists$term_n
        beta_0 <- out_list$beta_0
        strata_odds <- out_list$StrataOdds
        if ("Standard_Deviation" %in% names(out_list)) {
          stdev <- out_list$Standard_Deviation
          pval <- 2 * pnorm(-abs(beta_0 / stdev))
          res_table <- data.table(
            "Covariate" = names,
            "Subterm" = tforms,
            "Term Number" = term_n,
            "Central Estimate" = beta_0,
            "Standard Error" = stdev,
            "2-tail p-value" = pval
          )
        } else {
          res_table <- data.table(
            "Covariate" = names,
            "Subterm" = tforms,
            "Term Number" = term_n,
            "Central Estimate" = beta_0
          )
        }
        message("Final Results")
        print(res_table)
        deviance <- out_list$Deviation
        iteration <- out_list$Control_List$Iteration
        step_max <- out_list$Control_List$`Maximum Step`
        deriv_max <- out_list$Control_List$`Derivative Limiting`
        converged <- out_list$Converged
        #
        freepara <- out_list$FreeParameters
        freestrata <- out_list$FreeSets
        #
        message("\nMatched Case-Control Model Used")
        message(paste("Deviance: ", round(deviance, digits), sep = ""))
        message(paste(freestrata, " out of ", length(strata_odds), " matched sets used Unconditional Likelihood", sep = ""))
        if (!is.null(converged)) {
          message(paste("Iterations run: ", iteration, "\nmaximum step size: ", formatC(step_max, format = "e", digits = digits), ", maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
          if (converged) {
            message("Analysis converged")
          } else {
            message("Analysis did not converge, check convergence criteria or run further")
          }
        }
      } else {
        # get the model details
        names <- out_list$Parameter_Lists$names
        tforms <- out_list$Parameter_Lists$tforms
        term_n <- out_list$Parameter_Lists$term_n
        beta_0 <- out_list$beta_0
        if ("Standard_Deviation" %in% names(out_list)) {
          stdev <- out_list$Standard_Deviation
          pval <- 2 * pnorm(-abs(beta_0 / stdev))
          res_table <- data.table(
            "Covariate" = names,
            "Subterm" = tforms,
            "Term Number" = term_n,
            "Central Estimate" = beta_0,
            "Standard Error" = stdev,
            "2-tail p-value" = pval
          )
        } else {
          res_table <- data.table(
            "Covariate" = names,
            "Subterm" = tforms,
            "Term Number" = term_n,
            "Central Estimate" = beta_0
          )
        }
        message("Final Results")
        print(res_table)
        # get the model results
        LogLik <- out_list$LogLik
        AIC <- out_list$AIC
        BIC <- out_list$BIC
        deviation <- out_list$Deviation
        iteration <- out_list$Control_List$Iteration
        step_max <- out_list$Control_List$`Maximum Step`
        deriv_max <- out_list$Control_List$`Derivative Limiting`
        converged <- out_list$Converged
        if (is(out_list, "coxres")) {
          # cox model
          message("\nCox Model Used")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  AIC: ", round(AIC, digits), sep = ""))
        } else if (is(out_list, "poisres")) {
          # poisson model
          message("\nPoisson Model Used")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  Deviation: ", round(deviation, digits), ",  AIC: ", round(AIC, digits), ",  BIC: ", round(BIC, digits), sep = ""))
        } else if (is(out_list, "logitres")) {
          # logistic model
          message("\nLogisitic Model Used")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  Deviation: ", round(deviation, digits), ",  AIC: ", round(AIC, digits), ",  BIC: ", round(BIC, digits), sep = ""))
        } else {
          message("\nUnknown Model Used")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  AIC: ", round(AIC, digits), sep = ""))
        }
        if (!is.null(converged)) {
          message(paste("Iterations run: ", iteration, "\nmaximum step size: ", formatC(step_max, format = "e", digits = digits), ", maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
          if (converged) {
            message("Analysis converged")
          } else {
            message("Analysis did not converge, check convergence criteria or run further")
          }
        }
      }
    }
  } else {
    message(paste("Regression Failed"))
  }
  if ("RunTime" %in% names(out_list)) {
    run_time_sec <- as.numeric(out_list$RunTime, units = "secs")
    run_time_min <- as.numeric(out_list$RunTime, units = "mins")
    run_time_hour <- as.numeric(out_list$RunTime, units = "hours")
    if (run_time_sec < 60) {
      message(paste("Run finished in ", round(run_time_sec, digits), " seconds", sep = ""))
    } else if (run_time_min < 60) {
      message(paste("Run finished in ", round(run_time_min, digits), " minutes", sep = ""))
    } else {
      message(paste("Run finished in ", round(run_time_hour, digits), " hours", sep = ""))
    }
    #    message(paste("Run finished in ", out_list$RunTime))
  }
  message("|-------------------------------------------------------------------|")
}
