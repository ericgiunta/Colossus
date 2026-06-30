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
    temp_count <- temp_count + str_count(i, stringr::fixed("(")) - str_count(i, stringr::fixed(")")) # check to see if the current selections have the same number
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
  if (startsWith(string, "c(")) {
    sub_str <- substr(string, 3, nchar(string) - 1)
    args <- nested_split(sub_str) # get the entries of the vector
    args <- gsub("\"", "", args, fixed = TRUE) # remove string literals
    args <- unlist(lapply(args, parse_literal_string), use.names = FALSE) # make sure every entry is processed to a final type
    return(args)
  } else if (startsWith(string, "list(")) {
    sub_str <- substr(string, 6, nchar(string) - 1)
    args <- nested_split(sub_str) # make sure every entry is processed
    #
    factor_list <- list()
    for (i in seq_along(args)) {
      para_cur <- args[i]
      para_break <- lapply(strsplit(para_cur, "", fixed = TRUE), function(x) which(x == "="))[[1]] # split into the name and value
      if (length(para_break) == 0) {
        # no name, just add to list
        factor_list[[i]] <- parse_literal_string(para_cur) # process the item to the final value
      } else {
        item_name <- substr(para_cur, 1, para_break - 1)
        item_value <- substr(para_cur, para_break + 1, nchar(para_cur))
        item_name <- gsub("\"", "", item_name, fixed = TRUE)
        item_value <- gsub("\"", "", item_value, fixed = TRUE)
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
    # needs to detect numbers
    if (all(vapply(string, function(x) grepl("^[\\-]{0,1}[0-9]*\\.{0,1}[0-9]*$", x), logical(1))) || all(vapply(string, function(x) grepl("^[\\-]{0,1}[0-9]+e[\\-]{0,1}[0-9]+$", x), logical(1)))) {
      options(warn = 0) # checks for an integer, decimal, decimal places or scientific notation
      return(as.numeric(string))
    }
    options(warn = 0)
  }
  string
}

#' checks the input constraint matrix, to verify that it is valid
#'
#' \code{check_constraints} checks the input constraint matrix, applies fixes, and updates the model_control
#'
#' @noRd
#' @family Data Cleaning Functions
#' @return returns a list with the updated values
check_constraints <- function(a_n, model_control, cons_mat, cons_vec = c(0), verbose = 2) {
  # check that the constraint matrix is valid
  # It should be a numeric matrix
  if (is.null(cons_mat)) {
    stop("Error: Constraint matrix was empty.")
  } else if (suppressWarnings(any(is.character(cons_mat)))) {
    stop("Error: Constraint matrix was not numeric")
  } else if (is.Date(cons_mat)) {
    stop("Error: Constraint matrix was a date")
  } else if (!is.null(levels(cons_mat))) {
    stop("Error: Constraint matrix was a factor")
  } else if ((suppressWarnings(any(is.na(as.matrix(as.numeric(cons_mat)))))) || (is.list(cons_mat))) {
    stop("Error: Constraint matrix was not a numeric matrix.")
  } else {
    if (!is.matrix(cons_mat)) {
      row_count <- 1
    } else {
      row_count <- nrow(cons_mat)
    }
    cons_mat <- matrix(as.numeric(unlist(cons_mat)), nrow = row_count)
  }
  if (nrow(cons_mat) < 1) {
    stop("Error: Constraint matrix was empty.")
  }
  if (ncol(cons_mat) < length(a_n)) {
    if (verbose >= 2) {
      warning("Warning: The constraint matrix was missing columns.")
    }
    cons_mat <- cbind(cons_mat, matrix(rep(0, nrow(cons_mat) * (length(a_n) - ncol(cons_mat))), nrow = nrow(cons_mat)))
  } else if (ncol(cons_mat) > length(a_n)) {
    stop("Error: Constraint matrix had too many columns.")
  }
  if (!missing(cons_vec)) {
    if (is.null(cons_vec)) {
      stop("Error: Constraint vector was empty.")
    } else if (suppressWarnings(any(is.character(cons_vec)))) {
      stop("Error: Constraint vector was not numeric")
    } else if (is.Date(cons_vec)) {
      stop("Error: Constraint vector was a date")
    } else if (!is.null(levels(cons_vec))) {
      stop("Error: Constraint vector was a factor")
    } else if ((suppressWarnings(any(is.na(as.numeric(cons_vec))))) || (is.list(cons_vec))) {
      stop("Error: Constraint vector was not a numeric vector.")
    } else {
      cons_vec <- as.numeric(cons_vec)
    }
    if (length(cons_vec) < nrow(cons_mat)) {
      if (verbose >= 2) {
        warning("Warning: The constraint solution vector was too short. Zeros added for missing values.")
      }
      cons_vec <- c(cons_vec, rep(0.0, nrow(cons_mat) - length(cons_vec)))
    } else if (length(cons_vec) > nrow(cons_mat)) {
      stop("Error: The Constraint vector had too many entries")
    }
  } else {
    # Values may be set to zero by default, check each row
    for (i in 1:nrow(cons_mat)) {
      if (sum(cons_mat[i, ] != 0) == 1) {
        # There is only one non-zero entry in this row, and the solution vector defaults to zero
        stop("Error: Constraint solution would set parameter to zero by default. It is preferable to manually set cons_vec to zero if this is intended.")
      }
    }
    cons_vec <- c(rep(0.0, nrow(cons_mat)))
  }
  #
  row_remove <- c()
  for (i in 1:nrow(cons_mat)) {
    if (all(cons_mat[i, ] == rep(0, ncol(cons_mat)))) {
      row_remove <- c(row_remove, 0)
    } else {
      row_remove <- c(row_remove, NA)
    }
  }
  if (!all(is.na(row_remove))) {
    if (verbose >= 2) {
      warning("Warning: Atleast one row of the constraint matrix was all zero.")
    }
    cons_mat <- matrix(cons_mat[is.na(row_remove), ], nrow = nrow(cons_mat))
    cons_vec <- cons_vec[is.na(row_remove)]
  }
  # If we have more than one row, then they could be linearly dependent
  if (nrow(cons_mat) > 1) {
    # We want to simplify the matrix
    # Start by removing the linearly dependent rows and identifying incorrect rows
    total_mat <- cbind(cons_mat, unlist(cons_vec, use.names = FALSE))
    row_removed <- c(0)
    for (i in 2:nrow(total_mat)) {
      new_full_rank <- qr(total_mat[1:i, ])$rank # rank for the full solution matrix
      new_base_rank <- qr(cons_mat[1:i, ])$rank # rank for the constraint values
      if (new_full_rank == i - sum(row_removed)) {
        # the two full values were linearly independent
        # The original constraint values could still be linearly dependent
        if (new_base_rank == i - sum(row_removed)) {
          # it is all linearly independent
          row_removed <- c(row_removed, 0)
        } else {
          # The constraint matrix had linearly dependent values, but the vector portions were not equal
          row_removed <- c(row_removed, 1)
          stop("Error: Two linearly dependent constraint matrix rows had different solutions.")
        }
      } else {
        # the row was linearly dependent
        row_removed <- c(row_removed, 1)
      }
    }
    # Remove the offending rows
    total_mat <- total_mat[(row_removed == 0), ]
    # Convert to reduced row echelon format to simplify future calculations
    total_mat <- rref(total_mat)
    # Pull out the new matrix and vector
    cons_mat <- total_mat[, 1:(ncol(total_mat) - 1)]
    cons_vec <- total_mat[, ncol(total_mat)]
  }
  #
  if (nrow(cons_mat) * ncol(cons_mat) > 0) {
    model_control$constraint <- TRUE
  } else {
    if (verbose >= 2) {
      warning("Warning: Constraint matrix was given, but was not valid.")
    }
    model_control$constraint <- FALSE
  }
  #
  res <- list("control" = model_control, "mat" = cons_mat, "vec" = cons_vec)
}

#' Calculates and applies the interaction between a list of factor columns
#'
#' \code{Make_Interaction_Strata} iterates through a list of factors and finds all of the valid interactions
#'
#' @noRd
#' @param col_list vector of strata names
#' @param keep_base boolean if the baseline factor values should be kept
#' @param filter_df boolean if the dataframe should be filtered. If false, a new category is created for any 0 case strata
#' @family Data Cleaning Functions
#' @return returns a list with the data and new columns
Make_Interaction_Strata <- function(df, event0, col_list, control = list(verbose = TRUE), keep_base = TRUE, filter_df = TRUE) {
  vals <- col_list
  og_name <- names(df)
  combs <- c()
  if (length(vals) == 1) {
    # factor
    df$comb_strata <- as.integer(factor(df[[vals]])) - 1
    combs <- unique(df$comb_strata)
  } else {
    # there are multiple to combine
    # get the levels for each element
    if (filter_df) {
      for (term_i in seq_along(vals)) {
        factor_col <- vals[term_i]
        if (is.null(levels(df[[factor_col]]))) {
          df[[factor_col]] <- factor(df[[factor_col]]) # Only convert to factor if needed
        }
        i_levels <- levels(df[[factor_col]]) # get the levels of the factor
        if (!keep_base) {
          level_ref <- levels(df[[factor_col]])[1] # We can remove the baseline from the returned combinations
          i_levels <- i_levels[i_levels != level_ref]
        }
        # We only care about combinations that have events, so we can filter out some combinations early if the base level has no events
        for (col in i_levels) {
          temp <- sum(df[get(factor_col) == col, ][[event0]]) # get number of events
          if (temp == 0) { # if none then we remove that data and the level column
            if (control$verbose >= 2) {
              # nocov start
              warning(paste("Warning: no events for strata group:", col,
                sep = " "
              ))
              # nocov end
            }
            df <- df[get(factor_col) != col, ] # remove data
          }
        }
      }
    }
    # We want to combine the strata values together
    df$comb_strata <- ""
    for (term_i in seq_along(vals)) {
      factor_col <- vals[term_i]
      if (term_i == 1) {
        df$comb_strata <- df[[factor_col]]
      } else {
        df$comb_strata <- paste(df$comb_strata, df[[factor_col]], sep = ":")
      }
    }
    df$comb_strata <- as.integer(factor(df$comb_strata)) # converts to integer levels
    combs <- unique(df$comb_strata)
    if (filter_df) {
      df_end <- df[get(event0) == 1, ] # get the event data
      combs <- unique(df_end$comb_strata) # check for the strata with events
      comb_tot <- unique(df$comb)
      comb_remove <- comb_tot[!comb_tot %in% combs] # get the strata that were not in the event data
      for (comb in comb_remove) {
        df <- df["comb_strata" != comb, ] # remove data
      }
    }
  }
  list(data = df, combs = "comb_strata", levels = combs)
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
#'   UserID = c(112, 114, 213, 214, 115, 116, 117),
#'   Starting_Age = c(18, 20, 18, 19, 21, 20, 18),
#'   Ending_Age = c(30, 45, NA, 47, 36, NA, 55),
#'   Cancer_Status = c(0, 0, 1, 0, 1, 0, 0)
#' )
#' df <- Replace_Missing(df, c("Starting_Age", "Ending_Age"), 70)
Replace_Missing <- function(df, name_list, msv, verbose = FALSE) {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  verbose <- Check_Verbose(verbose)
  if (is.na(msv)) {
    stop("Error: The missing-value replacement is also NA")
  }
  for (j in name_list) {
    #
    if (j %in% names(df)) {
      # fine
    } else {
      stop("Error: ", j, " missing from column names")
    }
    if (sum(is.na(df[[j]]))) {
      data.table::set(df, which(is.na(df[[j]])), j, msv)
      # nocov start
      if (verbose >= 3) {
        message(paste("Note: Column ", j, " had replaced values",
          sep = ""
        ))
      }
      # nocov end
    }
  }
  df
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
    verbose = 0, lr = 0.75, maxiter = 20,
    halfmax = 5, epsilon = 1e-4, ll_epsilon = 1e-4,
    deriv_epsilon = 1e-4, step_max = 1.0,
    change_all = TRUE, thres_step_max = 100.0,
    ties = "breslow",
    ncores = as.numeric(detectCores())
  )
  #
  names(control) <- tolower(names(control)) # set the names to lowercase
  names(control) <- lapply(names(control), function(x) tryCatch(match.arg(x, choices = names(control_def)), error = function(error_message) x)) # match against expected values
  control <- control[!duplicated(names(control))] # filter down
  # nocov start
  if ((identical(Sys.getenv("TESTTHAT"), "true")) || (identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))) {
    control_def$ncores <- min(c(2, as.numeric(detectCores()))) # reduces cores if testing is running
  }
  if (Sys.getenv("ColossusOMP") == "") {
    syscheck <- System_Version()
    OpenMP <- syscheck[["OpenMP Enabled"]]
    if (!OpenMP) {
      Sys.setenv(ColossusOMP = "FALSE")
    } else {
      Sys.setenv(ColossusOMP = "TRUE")
      Sys.setenv(ColossusGCC = "TRUE")
      if (control["verbose"] > 1) {
        if (system.file(package = "processx") == "") {
          warning("Warning: processx is missing, which is used to check default compiler. Parallelization is disabled by default.")
        }
        if (system.file(package = "callr") == "") {
          warning("Warning: callr is missing, which is used to check default compiler. Parallelization is disabled by default.")
        }
      }
      os <- syscheck[["Operating System"]]
      if (os == "linux") {
        cpp_compiler <- syscheck[["Default c++"]]
        if (cpp_compiler != "") {
          if (cpp_compiler == "package_missing") {
            # just going to assume it will not work
            Sys.setenv(ColossusGCC = "FALSE") # nocov
          } else if (cpp_compiler == "gcc") {
            R_compiler <- syscheck[["R Compiler"]]
            if (R_compiler != "gcc") {
              Sys.setenv(ColossusGCC = "FALSE") # nocov
            }
          } else if (cpp_compiler == "clang") {
            Sys.setenv(ColossusGCC = "FALSE") # nocov
          }
        } else {
          R_compiler <- syscheck[["R Compiler"]]
          if (R_compiler != "gcc") {
            Sys.setenv(ColossusGCC = "FALSE") # nocov
          }
        }
      }
    }
  }
  # nocov end
  if ("verbose" %in% names(control)) {
    control$verbose <- Check_Verbose(control$verbose)
  } else {
    control["verbose"] <- control_def["verbose"]
  }
  # nocov start
  if (Sys.getenv("ColossusOMP") == "FALSE") {
    if (control["verbose"] > 1) {
      warning("Warning: OpenMP not detected, cores set to 1")
    }
    control$ncores <- 1 # nocov
  } else if ((Sys.getenv("R_COLOSSUS_NOT_CRAN") == "") && (Sys.getenv("ColossusGCC") == "FALSE")) {
    control$ncores <- 1 # nocov
    if (control["verbose"] > 1) {
      # In the past, a fedora machine running clang had issues with parallel code. Now linux machine using clang only use 1 core.
      warning("Warning: linux machine not using gcc, cores set to 1. Set R_COLOSSUS_NOT_CRAN environemnt variable to skip check")
    }
  }
  # nocov end
  # Take advantage of code to check control values
  control_args <- intersect(names(control), names(formals(ColossusControl)))
  control_def <- do.call(ColossusControl, control[control_args])
  control <- c(control[!(names(control) %in% names(control_def))], control_def)
  control
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
  control_def_names <- c(
    "single", "basic", "null", "cr", "linear_err",
    "gradient", "constraint", "strata", "surv",
    "schoenfeld", "risk",
    "risk_subset", "log_bound", "pearson", "deviance",
    "mcml", "observed_info", "time_risk",
    "logit_odds", "logit_ident", "logit_loglink",
    "logit_probit"
  )
  name_full_list <- c(control_def_names, "qchi", "alpha", "para_number", "maxstep", "manual", "search_mult", "step_size", "momentum", "adadelta", "adam", "momentum_decay", "learning_decay", "epsilon_decay", "constraint", "penalty_weight", "penalty_method")
  names(control) <- tolower(names(control)) # set the names to lowercase
  names(control) <- lapply(names(control), function(x) tryCatch(match.arg(x, choices = name_full_list), error = function(error_message) x)) # match against expected values
  control <- control[!duplicated(names(control))] # filter down
  for (nm in control_def_names) {
    if (nm %in% names(control)) {
      if (((length(control[[nm]]) > 1) || is.list(control[[nm]]))) {
        stop(paste0("Error: ", nm, " was not a single value."))
      }
      if (suppressWarnings(is.na(as.logical(control[nm])))) {
        stop(paste0("Error: ", nm, " had a non-logical value"))
      } else {
        control[nm] <- as.logical(control[nm])
      }
    } else {
      control[nm] <- FALSE
    }
  }
  link_vec <- c(control$logit_odds, control$logit_ident, control$logit_loglink, control$logit_probit)
  if (sum(link_vec) == 0) {
    control["logit_odds"] <- TRUE
  } else if (sum(link_vec) > 1) {
    stop("Error: Multiple link functions used, only select one link function")
  }
  if (control["null"] == TRUE) {
    control["single"] <- TRUE
  }
  if ("step_size" %in% names(control)) {
    if (suppressWarnings(any(is.na(as.numeric(control["step_size"]))))) {
      stop("Error: model control step size had a non-numeric value")
    } else {
      control["step_size"] <- as.numeric(control["step_size"])
    }
  } else {
    control["step_size"] <- 0.5
  }
  if ("unique_values" %in% names(control)) {
    if (suppressWarnings(any(is.na(as.numeric(control["unique_values"]))))) {
      stop("Error: unique values had a non-numeric value")
    } else {
      control["unique_values"] <- as.integer(control["unique_values"])
    }
  } else {
    control["unique_values"] <- 2
  }
  if ("gmix_theta" %in% names(control)) {
    if (suppressWarnings(!is(control$gmix_theta, "numeric"))) {
      stop("Error: gmix theta had a non-numeric value")
    }
    if (length(control$gmix_theta) > 1) {
      stop("Error: gmix theta had multiple values.")
    }
    if (abs(control$gmix_theta - 0.5) > 0.5) {
      warning("Warning: gmix theta was outside of [0,1.0]. Not the intended use, errors may occur.")
    }
  } else {
    control["gmix_theta"] <- 0.5
  }
  if ("gmix_term" %in% names(control)) {
    if (suppressWarnings(!is(control$gmix_term, "numeric"))) {
      if (is(control$gmix_term, "character")) {
        if (all(vapply(control$gmix_term, grepl, x = "er", logical(1)))) {
          control$gmix_term <- as.numeric(control$gmix_term == "e")
        } else {
          # Remaining entry was wrong
          stop("Error: Gmix term had an invalid option. Please only use 'e/r' or '0/1'.")
        }
      } else {
        stop("Error: Gmix term had an invalid option. Please only use 'e/r' or '0/1'.")
      }
    } else {
      control$gmix_term <- as.integer(control$gmix_term)
    }
    if (any(abs(control[["gmix_term"]] - 0.5) > 0.5)) {
      stop("Error: gmix term was outside of [0,1.0]. Not the intended use, errors will occur.")
    }
  } else {
    control["gmix_term"] <- c(0)
  }
  if ("conditional_threshold" %in% names(control)) {
    if (suppressWarnings(any(is.na(as.numeric(control["conditional_threshold"]))))) {
      stop("Error: conditional threshold had a non-numeric value")
    } else {
      control["conditional_threshold"] <- as.integer(control["conditional_threshold"])
    }
    control["conditional_threshold"] <- max(c(control[["conditional_threshold"]], 0))
  } else {
    control["conditional_threshold"] <- 50
  }
  if (control[["log_bound"]]) {
    if ("qchi" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["qchi"]))))) {
        stop("Error: chi squared value had a non-numeric value")
      } else {
        control["qchi"] <- as.numeric(control["qchi"])
      }
    } else {
      if ("alpha" %in% names(control)) {
        if (suppressWarnings(any(is.na(as.numeric(control["alpha"]))))) {
          stop("Error: alpha value had a non-numeric value")
        } else {
          control["alpha"] <- as.numeric(control["alpha"])
        }
        control["qchi"] <- qchisq(1 - control[["alpha"]], df = 1) / 2
      } else {
        control["alpha"] <- 0.05
        control["qchi"] <- qchisq(1 - control[["alpha"]], df = 1) / 2
      }
    }
    if ("para_num" %in% names(control)) {
      warning("Warning: para_num detected in model_control, did you mean para_number?") # nocov
      control["para_number"] <- control["para_num"]
    }
    if ("para_number" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["para_number"]))))) {
        stop("Error: parameter number was non-numeric")
      } else {
        control["para_number"] <- as.integer(control["para_number"])
      }
    } else {
      control["para_number"] <- 1
    }
    if ("maxstep" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["maxstep"]))))) {
        stop("Error: maximum step size was non-numeric")
      } else {
        control["maxstep"] <- as.integer(control["maxstep"])
      }
    } else {
      control["maxstep"] <- 10
    }
    if ("manual" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.logical(control["manual"]))))) {
        stop("Error: manual search was non-logical")
      } else {
        control["manual"] <- as.logical(control["manual"])
      }
    } else {
      control["manual"] <- FALSE
    }
    if ("search_mult" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["search_mult"]))))) {
        stop("Error: solution search multiplier was non-numeric")
      } else {
        control["search_mult"] <- as.numeric(control["search_mult"])
      }
    } else {
      control["search_mult"] <- 1.0
    }
    if ("step_size" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["step_size"]))))) {
        stop("Error: step size was non-numeric")
      } else {
        control["step_size"] <- as.numeric(control["step_size"])
      }
    } else {
      control["step_size"] <- 0.5
    }
    if ("bisect" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.logical(control["bisect"]))))) {
        stop("Error: bisection search was non-logical")
      } else {
        control["bisect"] <- as.logical(control["bisect"])
      }
    } else {
      control["bisect"] <- FALSE
    }
  }
  if (control[["gradient"]]) {
    control_def_names <- c(
      "momentum", "adadelta", "adam"
    )
    for (nm in control_def_names) {
      if (nm %in% names(control)) {
        if (suppressWarnings(any(is.na(as.logical(control[nm]))))) {
          stop(paste0("Error: ", nm, " was non-logical"))
        } else {
          control[nm] <- as.logical(control[nm])
        }
      } else {
        control[nm] <- FALSE
      }
    }
    # add different parameters
    if ("momentum_decay" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["momentum_decay"]))))) {
        stop("Error: momentum decay was non-numeric")
      } else {
        control["momentum_decay"] <- as.numeric(control["momentum_decay"])
      }
    } else {
      control["momentum_decay"] <- 0.9
    }
    if ("learning_decay" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["learning_decay"]))))) {
        stop("Error: learning decay was non-numeric")
      } else {
        control["learning_decay"] <- as.numeric(control["learning_decay"])
      }
    } else {
      control["learning_decay"] <- 0.999
    }
    if ("epsilon_decay" %in% names(control)) {
      if (suppressWarnings(any(is.na(as.numeric(control["epsilon_decay"]))))) {
        stop("Error: epsilon decay was non-numeric")
      } else {
        control["epsilon_decay"] <- as.numeric(control["epsilon_decay"])
      }
    } else {
      control["epsilon_decay"] <- 1e-4
    }
    if (control["constraint"] == TRUE) {
      # needs to add constraint penalty
      if ("penalty_weight" %in% names(control)) {
        if (suppressWarnings(any(is.na(as.numeric(control["penalty_weight"]))))) {
          stop("Error: penalty weight was non-numeric")
        } else {
          control["penalty_weight"] <- as.numeric(control["penalty_weight"])
        }
      } else {
        control["penalty_weight"] <- 1.0
      }
      if ("penalty_method" %in% names(control)) {
        if (suppressWarnings(any(is.na(as.character(control["penalty_method"]))))) {
          stop("Error: penalty method was not a string")
        } else {
          control["penalty_method"] <- as.character(control["penalty_method"])
        }
      } else {
        control["penalty_method"] <- "sqr_error"
      }
    }
  }
  control
}

#' Automatically checks the number of starting guesses
#'
#' \code{Check_Iters} checks the number of iterations and number of guesses, and corrects
#'
#' @noRd
#' @inheritParams R_template
#' @family Data Cleaning Functions
#' @return returns a list with the corrected control list and a_n
Check_Iters <- function(control, a_n) {
  if ("maxiters" %in% names(control)) {
    if (length(control$maxiters) == length(a_n) + 1) {
      # all good, it matches
    } else {
      if (control$verbose >= 3) {
        # nocov start
        message(paste("Note: Initial starts:", length(a_n),
          ", Number of iterations provided:",
          length(control$maxiters),
          ". Colossus requires one more iteration counts than number of guesses (for best guess)",
          sep = " "
        ))
        # nocov end
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
  list(control = control, a_n = a_n)
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
#' tforms <- list(cov_0 = "quad", cov_1 = "exp")
#' paras <- list(cov_0 = c(1, 3.45), cov_1 = c(1.2, 4.5, 0.1))
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
  full_paras
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
  b1
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
#' df <- data.table::data.table(a = a, b = b, c = c)
#' col_list <- c("c")
#' val <- factorize(df, col_list)
#' df <- val$df
#' new_col <- val$cols
#'
factorize <- function(df, col_list, verbose = 0) {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  col_list <- vapply(col_list, function(x) tryCatch(match.arg(x, choices = names(df)), error = function(error_message) x), USE.NAMES = FALSE, FUN.VALUE = "character")
  verbose <- Check_Verbose(verbose)
  cols <- c()
  col0 <- names(df)
  tnum <- c()
  for (i in seq_along(col_list)) {
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
  list(df = df, cols = cols)
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
#' a <- c(0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6)
#' b <- c(1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7)
#' c <- c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' d <- c(3, 4, 5, 6, 7, 8, 9, 1, 2, 1, 1, 2, 1, 2)
#' e <- c(1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2)
#' df <- data.table(a = a, b = b, c = c, d = d, e = e)
#' keep_constant <- c(0)
#' a_n <- c(-0.1, 0.1, 0.1, 0.2)
#' control <- list(ncores = 1, maxiter = 10, verbose = 0)
#' model <- Cox(a, b, c) ~ plinear(d * d, 0) + loglinear(factor(e))
#' alternative_model <- CoxRun(model, df,
#'   control = control,
#'   a_n = a_n, keep_constant = c(0, 1, 0)
#' )
#' null_model <- CoxRun(Cox(a, b, c) ~ null(), df, control = control)
#' score <- Likelihood_Ratio_Test(alternative_model, null_model)
#'
Likelihood_Ratio_Test <- function(alternative_model, null_model) {
  alt_is_null <- alternative_model$modelcontrol$null
  null_is_null <- null_model$modelcontrol$null
  if (("LogLik" %in% names(alternative_model)) && ("LogLik" %in% names(null_model))) {
    #
    alt_is_null <- alternative_model$modelcontrol$null
    null_is_null <- null_model$modelcontrol$null
    if (alt_is_null) {
      stop("Error: Alternative model shouldn't be null.")
    } else {
      alt_count <- length(alternative_model$beta_0) - sum(alternative_model$model$keep_constant)
    }
    if (alt_is_null) {
      null_count <- 0
    } else {
      null_count <- length(null_model$beta_0) - sum(null_model$model$keep_constant)
    }
    #
    freedom <- alt_count - null_count
    val <- 2 * (alternative_model$LogLik - null_model$LogLik)
    pval <- pchisq(val, freedom, lower.tail = FALSE)
    return(list(Difference = val, `p value` = pval))
  } else {
    stop("Error: models input did not contain LogLik values")
  }
  NULL # nocov
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
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  verbose <- Check_Verbose(verbose)
  cols <- vapply(cols, function(x) tryCatch(match.arg(x, choices = names(df)), error = function(error_message) x), USE.NAMES = FALSE, FUN.VALUE = "character")
  if (length(cols) > 1) {
    features_pair <- combn(cols, 2, simplify = FALSE) # list all column pairs
    terms_pair <- combn(term_n, 2, simplify = FALSE) # list all term pairs
    toRemove <- c() # init a vector to store duplicates
    for (pair_n in seq_along(features_pair)) {
      # put the pairs for testing into temp objects
      pair <- unlist(features_pair[pair_n])
      term <- unlist(terms_pair[pair_n])
      f1 <- pair[1]
      f2 <- pair[2]
      if (!(f1 %in% names(df))) {
        stop("Error: ", f1, " not in data.table") # nocov
      }
      if (!(f2 %in% names(df))) {
        stop("Error: ", f2, " not in data.table") # nocov
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
              # nocov start
              warning(paste("Warning: ", f1, " and ", f2,
                " are equal",
                sep = ""
              ))
              # nocov end
            }
            toRemove <- c(toRemove, f2) # build the list of duplicates
          }
          if (min(df[[f2]]) == max(df[[f2]])) {
            if (min(df[[f2]]) == 0) {
              if (verbose >= 2) {
                # nocov start
                warning(paste("Warning: ", f2,
                  " is equal to zero, removed.",
                  sep = ""
                ))
                # nocov end
              }
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
      stop("Error: ", f1, " not in data.table")
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
  c() # nocov
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
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  # We need to make sure that the variable is a string, not length 0, not empty string, and not null
  name_vec <- c("starting age", "ending age", "event")
  for (i in 1:3) {
    name <- name_vec[i]
    if (!is(ce[[i]], "character")) {
      stop(paste0("Error: The ", name, " column must be a string")) # nocov
    }
    if (length(ce[[i]]) == 0) {
      stop(paste0("Error: The ", name, " column must not be empty")) # nocov
    }
    if (length(ce[[i]]) > 1) {
      stop(paste0("Error: The ", name, " column had multiple values")) # nocov
    }
    if ((ce[[i]] == "") || (is.null(ce[[i]])) || (is.na(ce[[i]]))) {
      stop(paste0("Error: The ", name, " column must not be empty")) # nocov
    }
  }
  if (ce[1] == ce[2]) {
    stop("Error: The starting and ending interval times were set to the same column, they must be different") # nocov
  }
  #
  verbose <- Check_Verbose(verbose)
  if (ce[1] %in% c("%trunc%", "right_trunc")) {
    if (ce[2] %in% c("%trunc%", "left_trunc")) {
      stop("Error: Both endpoints are truncated, not acceptable")
    }
    tname <- ce[2]
    if (!is.numeric(df[[tname]])) {
      stop("Error: Age column was not numeric: ", tname)
    }
    tmin <- min(df[, get(tname)]) - 1
    if (!("right_trunc" %in% names(df))) {
      df[, ":="(right_trunc = tmin)]
    }
    ce[1] <- "right_trunc"
  } else if (ce[2] %in% c("%trunc%", "left_trunc")) {
    tname <- ce[1]
    if (!is.numeric(df[[tname]])) {
      stop("Error: Age column was not numeric: ", tname)
    }
    tmax <- max(df[, get(tname)]) + 1
    if (!("left_trunc" %in% names(df))) {
      df[, ":="(left_trunc = tmax)]
    }
    ce[2] <- "left_trunc"
  }
  list(df = df, ce = ce)
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
#' df <- data.table::data.table(a = a, b = b, c = c)
#' time1 <- "%trunc%"
#' time2 <- "a"
#' event <- "c"
#' control <- list(
#'   lr = 0.75, maxiter = -1, halfmax = 5, epsilon = 1e-9,
#'   deriv_epsilon = 1e-9, step_max = 1.0,
#'   thres_step_max = 100.0,
#'   verbose = FALSE, ties = "breslow", double_step = 1
#' )
#' grt_f <- function(df, time_col) {
#'   return((df[, "b"] * df[, get(time_col)])[[1]])
#' }
#' func_form <- c("lin")
#' df_new <- gen_time_dep(
#'   df, time1, time2, event, TRUE, 0.01, c("grt"), c(),
#'   c(grt_f), paste("test", "_new.csv", sep = ""), func_form, 1
#' )
#' file.remove("test_new.csv")
#'
gen_time_dep <- function(df, time1, time2, event0, iscox, dt, new_names, dep_cols, func_form, fname, tform, nthreads = as.numeric(detectCores())) {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  # ------------------------------------------------------------------------------ #
  # Make data.table use the set number of threads too
  if ((identical(Sys.getenv("TESTTHAT"), "true")) || (identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))) {
    nthreads <- min(c(2, nthreads))
  }
  thread_0 <- setDTthreads(nthreads) # save the old number and set the new number
  # ------------------------------------------------------------------------------ #
  dfn <- names(df)
  ce <- c(time1, time2, event0)
  t_check <- Check_Trunc(df, ce)
  df <- t_check$df
  ce <- t_check$ce
  time1 <- ce[1]
  time2 <- ce[2]
  dfn_same <- dfn[!(dfn %in% dep_cols)]
  dfn_dep <- c()
  if (length(new_names) != length(func_form)) {
    stop("Error: new_names vector should be the same size as the list of functions applied")
  }
  if (length(new_names) != length(tform)) {
    stop("Error: new_names vector should be the same size as the list of interpolation method used")
  }
  func_form <- c(func_form)
  for (i in seq_along(new_names)) {
    name0 <- paste(new_names[i], 0, sep = "_")
    name1 <- paste(new_names[i], 1, sep = "_")
    func <- func_form[i]
    df[, name0] <- lapply(func, function(f) f(df, time1))
    df[, name1] <- lapply(func, function(f) f(df, time2))
    dfn_dep <- c(dfn_dep, name0, name1)
  }
  for (i in seq_along(tform)) {
    temp <- tform[i]
    if (temp != "lin") {
      a <- substr(temp, 1, 5)
      if (a != "step?") {
        stop("Error: Interpolation method not recognized: ", temp)
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
  # Revert data.table core change
  thread_1 <- setDTthreads(thread_0) # revert the old number
  # ------------------------------------------------------------------------------ #
  df_new
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
#' df <- data.table::data.table(m0 = m0, m1 = m1, d0 = d0, d1 = d1, y0 = y0, y1 = y1)
#' df <- Date_Shift(df, c("m0", "d0", "y0"), c("m1", "d1", "y1"), "date_since")
#'
Date_Shift <- function(df, dcol0, dcol1, col_name, units = "days") {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
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
  df[, def_cols, with = FALSE]
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
#'   m0 = m0, m1 = m1,
#'   d0 = d0, d1 = d1,
#'   y0 = y0, y1 = y1
#' )
#' tref <- strptime("3-22-1997", format = "%m-%d-%Y", tz = "UTC")
#' df <- Time_Since(df, c("m1", "d1", "y1"), tref, "date_since")
#'
Time_Since <- function(df, dcol0, tref, col_name, units = "days") {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
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
  df[, def_cols, with = FALSE]
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
#' df <- data.table(t0 = a, t1 = b, e0 = c, e1 = d, fac = e)
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
#' name_list <- list(shared = names_shared, e0 = names_e0, e1 = names_e1)
#' term_n_list <- list(shared = term_n_shared, e0 = term_n_e0, e1 = term_n_e1)
#' tform_list <- list(shared = tform_shared, e0 = tform_e0, e1 = tform_e1)
#' keep_constant_list <- list(
#'   shared = keep_constant_shared,
#'   e0 = keep_constant_e0, e1 = keep_constant_e1
#' )
#' a_n_list <- list(shared = a_n_shared, e0 = a_n_e0, e1 = a_n_e1)
#' val <- Joint_Multiple_Events(
#'   df, events, name_list, term_n_list,
#'   tform_list, keep_constant_list, a_n_list
#' )
#'
Joint_Multiple_Events <- function(df, events, name_list, term_n_list = list(), tform_list = list(), keep_constant_list = list(), a_n_list = list()) {
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  # filling missing values
  for (i in names(name_list)) {
    temp0 <- unlist(name_list[i], use.names = FALSE)
    if (i %in% names(term_n_list)) {
      temp1 <- unlist(term_n_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(
          "Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in term_n_list has ",
          length(temp1),
          " items. Omit entry in term_n_list to set to default of term 0 or add missing values"
        )
      }
    } else {
      temp <- list(rep(0, length(temp0)))
      names(temp) <- i
      term_n_list <- c(term_n_list, temp)
    }
    if (i %in% names(tform_list)) {
      temp1 <- unlist(tform_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(
          "Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in tform_list has ",
          length(temp1),
          " items. Omit entry in tform_list to set to default of 'loglin' or add missing values"
        )
      }
    } else {
      temp <- list(rep("loglin", length(temp0)))
      names(temp) <- i
      tform_list <- c(tform_list, temp)
    }
    if (i %in% names(keep_constant_list)) {
      temp1 <- unlist(keep_constant_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(
          "Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in keep_constant_list has ",
          length(temp1),
          " items. Omit entry in tform_list to set to default of 0 or add missing values"
        )
      }
    } else {
      temp <- list(rep(0, length(temp0)))
      names(temp) <- i
      keep_constant_list <- c(keep_constant_list, temp)
    }
    if (i %in% names(a_n_list)) {
      temp1 <- unlist(a_n_list[i], use.names = FALSE)
      if (length(temp0) != length(temp1)) {
        stop(
          "Error: item ", i, " in name_list has ",
          length(temp0),
          " items, but same item in a_n_list has ",
          length(temp1),
          " items. Omit entry in a_n_list to set to default of 0 or add missing values"
        )
      }
    } else {
      temp <- list(rep(0.1, length(temp0)))
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
  list(df = df0, names = names, term_n = term_n, tform = tform, keep_constant = keep_constant, a_n = a_n)
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
  # nocov start
  if (class(df)[[1]] != "data.table") {
    tryCatch(
      {
        setDT(df)
      },
      error = function(error_message) {
        df <- data.table(df)
      }
    )
  }
  # nocov end
  verbose <- Check_Verbose(verbose)
  cols <- c()
  for (i in seq_along(interactions)) {
    interac <- interactions[i]
    formula <- unlist(strsplit(interac, "?", fixed = TRUE), use.names = FALSE)
    if (length(formula) != 3) {
      stop(
        "Error: Interaction:", interac, "has incorrect length of",
        length(formula), "but should be 3."
      )
    }
    newcol <- paste(formula[1], formula[2], formula[3], sep = "")
    if (new_names[i] != "") {
      newcol <- new_names[i]
    }
    col1 <- formula[1]
    col2 <- formula[3]
    if (paste(formula[1], "?", formula[2], "?", formula[3], sep = "") %in% interactions[i + seq_along(interactions)]) {
      if (verbose >= 2) {
        warning(paste("Warning: interation ", i, "is duplicated")) # nocov
      }
    } else if (paste(formula[3], "?", formula[2], "?", formula[1], sep = "") %in% interactions[i + seq_along(interactions)]) {
      if (verbose >= 2) {
        # nocov start
        warning(paste(
          "Warning: the reverse of interation ", i,
          "is duplicated"
        ))
        # nocov end
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
        stop("Error: Incorrect operation of", formula[2])
      }
    }
  }
  list(df = df, cols = cols)
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
#' @noRd
#' @family Data Cleaning Functions
#' @return returns list with the normalized values
apply_norm <- function(df, norm, names, input, values, model_control) {
  if ((!is(norm, "character")) || (length(norm) > 1)) {
    stop("Error: normalization factor was not an individual string.")
  }
  if (input) {
    a_n <- values$a_n
    cons_mat <- values$cons_mat
    tforms <- values$tform
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
            warning(paste("Warning: Maximum value for ", names[i], " was 0. Normalization not applied to column.", sep = "")) # nocov
            val <- 1.0
          } else if (tforms[i] == "step_slope") {
            # Forcing to 1, no need to normalize this one
            val <- 1.0
          }
          norm_weight <- c(norm_weight, val)
        }
      } else if (tolower(norm) == "mean") {
        for (i in seq_along(names)) {
          val <- summarise(df, mean_value = mean(get(names[i])))[[1]]
          if (val == 0.0) {
            warning(paste("Warning: Average value for ", names[i], " was 0. Normalization not applied to column.", sep = "")) # nocov
            val <- 1.0
          } else if (tforms[i] == "step_slope") {
            # Forcing to 1, no need to normalize this one
            val <- 1.0
          }
          norm_weight <- c(norm_weight, val)
        }
      } else {
        stop(gettextf(
          "Error: Normalization arguement '%s' not valid.",
          norm
        ), domain = NA)
      }
      for (i in seq_along(names)) {
        if (typeof(a_n) != "list") {
          if (grepl("_int", tforms[i], fixed = TRUE)) {
            a_n[i] <- a_n[i] / norm_weight[i]
          } else {
            a_n[i] <- a_n[i] * norm_weight[i]
          }
        } else {
          for (j in seq_along(a_n)) {
            if (grepl("_int", tforms[i], fixed = TRUE)) {
              a_n[[j]][i] <- a_n[[j]][i] / norm_weight[i]
            } else {
              a_n[[j]][i] <- a_n[[j]][i] * norm_weight[i]
            }
          }
        }
        if (match(names[i], names)[1] == i) {
          df[, names[i]] <- df[, names[i], with = FALSE] / norm_weight[i]
        }
      }
      if (model_control[["constraint"]] == TRUE) {
        for (i in seq_along(names)) {
          if (grepl("_int", tforms[i], fixed = TRUE)) {
            cons_mat[, i] <- cons_mat[, i] * norm_weight[i]
          } else {
            cons_mat[, i] <- cons_mat[, i] / norm_weight[i]
          }
        }
      }
    } else {
      stop(gettextf(
        "Error: Normalization arguement '%s' not valid.",
        norm
      ), domain = NA)
    }
    output <- list(
      a_n = a_n,
      cons_mat = cons_mat,
      norm_weight = norm_weight,
      df = df
    )
  } else {
    res <- values$output
    norm_weight <- values$norm_weight
    keep_constant <- res$Parameter_Lists$keep_constant
    tforms <- values$tform
    if (any(norm_weight != 1.0)) {
      if (tolower(norm) == "null") {
        # nothing changes
      } else if (tolower(norm) %in% c("mean", "max")) {
        # weight by the maximum value
        if (model_control$single) {
          for (i in seq_along(names)) {
            if (grepl("_int", tforms[i], fixed = TRUE)) {
              res$beta_0[i] <- res$beta_0[i] * norm_weight[i]
            } else {
              res$beta_0[i] <- res$beta_0[i] / norm_weight[i]
            }
          }
        } else {
          for (i in seq_along(names)) {
            if (keep_constant[i] == 0) {
              i_der <- i - sum(head(keep_constant, i))
              if (grepl("_int", tforms[i], fixed = TRUE)) {
                res$First_Der[i_der] <- res$First_Der[i_der] / norm_weight[i]
                res$beta_0[i] <- res$beta_0[i] * norm_weight[i]
                res$Standard_Error[i] <- res$Standard_Error[i] * norm_weight[i]
              } else {
                res$First_Der[i_der] <- res$First_Der[i_der] * norm_weight[i]
                res$beta_0[i] <- res$beta_0[i] / norm_weight[i]
                res$Standard_Error[i] <- res$Standard_Error[i] / norm_weight[i]
              }
              for (j in seq_along(names)) {
                if (keep_constant[j] == 0) {
                  j_der <- j - sum(head(keep_constant, j))
                  if (grepl("_int", tforms[i], fixed = TRUE)) {
                    if (grepl("_int", tforms[j], fixed = TRUE)) {
                      res$Second_Der[i_der, j_der] <- res$Second_Der[i_der, j_der] / norm_weight[i] / norm_weight[j]
                      res$Covariance[i_der, j_der] <- res$Covariance[i_der, j_der] * norm_weight[i] * norm_weight[j]
                    } else {
                      res$Second_Der[i_der, j_der] <- res$Second_Der[i_der, j_der] / norm_weight[i] * norm_weight[j]
                      res$Covariance[i_der, j_der] <- res$Covariance[i_der, j_der] * norm_weight[i] / norm_weight[j]
                    }
                  } else {
                    if (grepl("_int", tforms[j], fixed = TRUE)) {
                      res$Second_Der[i_der, j_der] <- res$Second_Der[i_der, j_der] * norm_weight[i] / norm_weight[j]
                      res$Covariance[i_der, j_der] <- res$Covariance[i_der, j_der] / norm_weight[i] * norm_weight[j]
                    } else {
                      res$Second_Der[i_der, j_der] <- res$Second_Der[i_der, j_der] * norm_weight[i] * norm_weight[j]
                      res$Covariance[i_der, j_der] <- res$Covariance[i_der, j_der] / norm_weight[i] / norm_weight[j]
                    }
                  }
                }
              }
            } else {
              if (grepl("_int", tforms[i], fixed = TRUE)) {
                res$beta_0[i] <- res$beta_0[i] * norm_weight[i]
              } else {
                res$beta_0[i] <- res$beta_0[i] / norm_weight[i]
              }
            }
          }
          if (model_control[["constraint"]] == TRUE) {
            for (i in seq_along(names)) {
              if (grepl("_int", tforms[i], fixed = TRUE)) {
                res$constraint_matrix[, i] <- res$constraint_matrix[, i] / norm_weight[i]
              } else {
                res$constraint_matrix[, i] <- res$constraint_matrix[, i] * norm_weight[i]
              }
            }
          }
        }
      } else {
        stop(gettextf(
          "Error: Normalization arguement '%s' not valid.",
          norm
        ), domain = NA)
      }
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
    if (grepl("^darwin", R.version$os, fixed = TRUE)) { # nocov
      os <- "osx" # nocov
    }
    if (grepl("linux-gnu", R.version$os, fixed = TRUE)) { # nocov
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
  if (system.file(package = "processx") == "") {
    out <- "package_missing"
  } else {
    out <- tryCatch(run("c++", "-v"),
      error = function(cnd) list(stdout = "")
    )
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
  if (system.file(package = "callr") == "") {
    out <- "package_missing"
  } else {
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
  if (system.file(package = "processx") == "") {
    out <- "package_missing"
  } else {
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
  #
  omp_check <- FALSE
  gcc_check <- FALSE
  if (!OMP) {
    omp_check <- FALSE
  } else {
    omp_check <- TRUE
    gcc_check <- TRUE
    if (os == "linux") {
      if (gcc != "") {
        if (gcc == "package_missing") {
          # just going to assume it will not work
          gcc_check <- FALSE
        } else if (cpp_compiler == "gcc") {
          if (Rcomp != "gcc") {
            gcc_check <- FALSE
          }
        } else if (cpp_compiler == "clang") {
          gcc_check <- FALSE
        }
      } else {
        if (Rcomp != "gcc") {
          gcc_check <- FALSE
        }
      }
    }
  }
  #
  omp_allowed <- TRUE
  if (!omp_check) {
    omp_allowed <- FALSE
  } else if (!gcc_check) {
    omp_allowed <- FALSE
  }
  #
  list(
    `Operating System` = os, `Default c++` = gcc, `R Compiler` = Rcomp,
    `OpenMP Enabled` = OMP, `Multicore Allowed` = omp_allowed
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
#'   a = a,
#'   b = b,
#'   c = c
#' )
#' categ <- list(
#'   a = "0/3/5]7",
#'   b = list(
#'     lower = c(-1, 3, 6),
#'     upper = c(3, 6, 10),
#'     name = c("low", "medium", "high")
#'   )
#' )
#' event <- list(
#'   c = "count AS cases",
#'   a = "mean", b = "mean"
#' )
#' e <- Event_Count_Gen(table, categ, event, T)
#'
Event_Count_Gen <- function(table, categ = list(), events = list(), verbose = FALSE) {
  df <- as_tibble(table)
  `%>%` <- dplyr::`%>%`
  #
  categ_cols <- c()
  categ_bounds <- list()
  names(categ) <- lapply(names(categ), function(x) tryCatch(match.arg(x, choices = names(table)), error = function(error_message) x)) # match against columns in the table
  for (cat in names(categ)) {
    cat_str <- ""
    if (!cat %in% names(table)) {
      stop("Error: ", cat, " not in table")
    }
    if (length(categ[[cat]]) > 1) { # list of bounds
      names(categ[[cat]]) <- tolower(names(categ[[cat]])) # set the names to lowercase
      temp0 <- categ[[cat]]$lower
      temp1 <- categ[[cat]]$upper
      if ("name" %in% names(categ[[cat]])) { # assign names to the levels
        temp2 <- categ[[cat]]$name
      } else { # number the categories
        temp2 <- seq_along(temp0)
      }
      num_categ <- length(temp0) # number of categories
      cat_col <- paste(cat, "category", sep = "_") # name of the category
      categ_cols <- c(categ_cols, cat_col) # add to list
      df <- df %>% mutate("{cat_col}" := "Unassigned") # add to tibble
      for (i in 1:num_categ) { # for each category
        L <- as.numeric(temp0[i]) # lower bound
        if (grepl("]", temp1[i], fixed = TRUE)) {
          U <- as.numeric(gsub("]", "", temp1[i], fixed = TRUE)) # assign upper
          a_col_categ <- case_when(df[[cat]] <= U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]]) # rows within the bounds and unassigned are assigned the name
          cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ") # add to list of intervals
        } else {
          U <- as.numeric(temp1[i]) # assign upper
          if (L > U) {
            stop(paste("Error"))
          }
          if (L == U) { # discrete case
            a_col_categ <- case_when(df[[cat]] == U & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          } else { # interval case
            a_col_categ <- case_when(df[[cat]] < U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ as.character(temp2[i]), .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
          }
        }
        df[[cat_col]] <- a_col_categ # update column
      }
      categ_bounds[[cat_col]] <- cat_str # add to list of intervals
    } else {
      if (tolower(categ[[cat]]) == "factor") {
        # It is a factor
        cat_col <- paste(cat, "category", sep = "_")
        categ_cols <- c(categ_cols, cat_col)
        df[[cat_col]] <- as.factor(df[[cat]])
        cat_str <- paste(levels(df[[cat_col]]), collapse = "/")
        categ_bounds[[cat_col]] <- cat_str # add to list of intervals
      } else {
        # string of bounds
        cat_str <- ""
        temp <- categ[[cat]]
        temp <- gsub("/", " / ", temp, fixed = TRUE)
        temp <- gsub("]", " ] ", temp, fixed = TRUE)
        temp <- strsplit(temp, "\\s+")[[1]] # adding and splitting by spaces
        match_index <- which(temp %in% c("/", "]"))
        cat_col <- paste(cat, "category", sep = "_")
        categ_cols <- c(categ_cols, cat_col)
        df <- df %>% mutate("{cat_col}" := "Unassigned") # add to tibble
        match_i <- 1
        for (time_i in match_index) {
          # We want to figure out if there is a label
          L <- as.numeric(temp[time_i - 1])
          if (time_i == match_index[length(match_index)]) {
            # We are in the last entry
            if (length(temp) - time_i == 2) {
              # There is a label
              entry_label <- temp[time_i + 1]
              U <- as.numeric(temp[time_i + 2])
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              entry_label <- paste(L, U, sep = " - ")
            }
          } else {
            # Not at last entry
            if (match_index[match_i + 1] - time_i == 3) {
              # There is a label
              entry_label <- temp[time_i + 1]
              U <- as.numeric(temp[time_i + 2])
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              entry_label <- paste(L, U, sep = " - ")
            }
          }
          if (L == U) { # dicrete value case
            a_col_categ <- case_when(df[[cat]] == U & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          } else if (temp[time_i] == "/") { # strict less than upper bound
            a_col_categ <- case_when(df[[cat]] < U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
          } else { # within and including bounds
            a_col_categ <- case_when(df[[cat]] <= U & df[[cat]] >= L & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          }
          df[[cat_col]] <- a_col_categ # update column
          match_i <- match_i + 1
        }
        categ_bounds[[cat_col]] <- cat_str # add interval boundaries
      }
    }
  }
  df_group <- df %>%
    group_by(across(all_of(categ_cols))) %>%
    summarize(COUNT = n(), .groups = "drop") # group by columns and summarize by counts
  names(events) <- lapply(names(events), function(x) tryCatch(match.arg(x, choices = names(table)), error = function(error_message) x)) # match against columns in the table
  for (evt in names(events)) { # iterate through events
    if (grepl(" AS ", events[[evt]], fixed = TRUE)) { # get method and updated name
      temp <- gsub(" AS ", " ", events[[evt]], fixed = TRUE)
      temp <- strsplit(temp, "\\s+")[[1]]
      col_name <- temp[2]
      method <- temp[1]
    } else {
      col_name <- evt
      temp <- strsplit(events[[evt]], "\\s+")[[1]]
      method <- temp[[length(temp)]]
    }
    #
    method <- tolower(method)
    method <- lapply(method, function(x) tryCatch(match.arg(x, choices = c("count", "sum", "mean", "rsum", "rmean")), error = function(error_message) x))[[1]] # match against expected values
    #
    if (method %in% c("count", "sum", "rsum")) { # summarize by count
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method %in% c("mean", "rmean")) { # summarize by mean
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := mean(.data[[evt]]), .groups = "drop")
    } else {
      stop(paste("Error: ", method, " is not a recognized aggregation method", sep = ""))
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
#'   a = a, b = b, c = c,
#'   d = d, e = e, f = f,
#'   g = g, h = h, i = i
#' )
#' categ <- list(
#'   a = "-1/3/5]7"
#' )
#' summary <- list(
#'   c = "count AS cases"
#' )
#' events <- list("c")
#' pyr <- list(
#'   entry = list(year = "f"),
#'   exit = list(year = "i"),
#'   unit = "years"
#' )
#' e <- Event_Time_Gen(table, pyr, categ, summary, events)
#'
Event_Time_Gen <- function(table, pyr = list(), categ = list(), summaries = list(), events = c(), verbose = FALSE) {
  df <- as_tibble(table)
  `%within%` <- lubridate::`%within%`
  `%>%` <- dplyr::`%>%`
  #
  table_names <- names(table)
  # Setting default date values
  year_default <- 1900
  month_default <- 1
  day_default <- 1
  # Checking for errors or valid data
  if (length(events) == 0) {
    stop("Error: no events were given")
  }
  event_cols <- c()
  event_names <- c()
  evt_list <- list()
  for (evt in events) {
    if (grepl(" AS ", evt, fixed = TRUE)) { # get the column and name
      temp <- strsplit(gsub(" AS ", " ", evt, fixed = TRUE), "\\s+")[[1]]
      evt_col <- temp[2]
      evt_df <- temp[1]
    } else {
      evt_df <- evt
      evt_col <- evt
    }
    event_cols <- c(event_cols, evt_df)
    event_names <- c(event_names, evt_col)
    evt_list[evt_df] <- paste("count AS ", evt_col, sep = "")
  }
  summaries <- c(evt_list, summaries[!(names(summaries) %in% names(evt_list))])
  names(summaries) <- lapply(names(summaries), function(x) tryCatch(match.arg(x, choices = table_names), error = function(error_message) x)) # match against columns in the table
  #
  names(pyr) <- tolower(names(pyr)) # set the names to lowercase
  #
  categ_cols <- c()
  categ_bounds <- list()
  for (cat in names(categ)) { # for each category
    cat_str <- ""
    if (grepl(" AS ", cat, fixed = TRUE)) { # get the column and name
      temp <- strsplit(gsub(" AS ", " ", cat, fixed = TRUE), "\\s+")[[1]]
      cat_col <- temp[2]
      cat_df <- temp[1]
    } else {
      cat_df <- cat
      cat_col <- paste(cat, "category", sep = "_")
    }
    cat_df <- lapply(cat_df, function(x) tryCatch(match.arg(x, choices = table_names), error = function(error_message) x))[[1]] # match against columns in the table
    if (cat_col %in% names(df)) { # check that the category doesn't already exist in the original dataframe
      stop("Error: ", cat_col, " already exists, use ' AS ' to rename if needed")
    }
    if (!is.null(names(categ[[cat]]))) { # boundary as lists
      names(categ[[cat]]) <- tolower(names(categ[[cat]])) # set the names to lowercase
      if ("lower" %in% names(categ[[cat]])) { # lower and upper boundary intervals
        if (!cat_df %in% table_names) {
          stop("Error: ", cat_df, " not in table")
        }
        temp0 <- categ[[cat]]$lower
        temp1 <- categ[[cat]]$upper
        if ("name" %in% names(categ[[cat]])) { # check for names for each level
          temp2 <- categ[[cat]]$name
        } else {
          temp2 <- seq_along(temp0)
        }
        num_categ <- length(temp0)
        categ_cols <- c(categ_cols, cat_col)
        df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
        for (i in 1:num_categ) { # for each level
          L <- as.numeric(temp0[i])
          if (grepl("]", temp1[i], fixed = TRUE)) { # check for including the upper limit
            U <- as.numeric(gsub("]", "", temp1[i], fixed = TRUE)) # get upper limit
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
            month <- rep(month_default, length(day)) # nocov
          }
          if ("year" %in% names(categ[[cat]])) { # nocov
            year <- categ[[cat]]$year # nocov
          } else {
            year <- rep(year_default, length(day)) # nocov
          }
        } else if ("month" %in% names(categ[[cat]])) { # nocov
          month <- categ[[cat]]$month # nocov
          day <- rep(day_default, length(month)) # nocov
          if ("year" %in% names(categ[[cat]])) { # nocov
            year <- categ[[cat]]$year # nocov
          } else { # nocov
            year <- rep(year_default, length(month)) # nocov
          }
        } else if ("year" %in% names(categ[[cat]])) { # nocov
          year <- categ[[cat]]$year # nocov
          day <- rep(day_default, length(year)) # nocov
          month <- rep(month_default, length(year)) # nocov
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
        entry <- make_date(day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df)))
        exit <- make_date(day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df)))
        if ("exit" %in% names(pyr)) { # format the exit as date column
          pyr_exit <- pyr$exit
          names(pyr_exit) <- tolower(names(pyr_exit)) # set the names to lowercase
          if ("year" %in% names(pyr_exit)) { # nocov
            if ("month" %in% names(pyr_exit)) { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
              } else { # nocov
                exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = rep(day_default, nrow(df))) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(year = df[[pyr_exit$year]], day = df[[pyr_exit$day]], month = rep(month_default, nrow(df))) # nocov
              } else { # nocov
                exit <- make_date(year = df[[pyr_exit$year]], day = rep(day_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
              }
            }
          } else {
            if ("month" %in% names(pyr_exit)) { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(month = df[[pyr_exit$month]], day = df[[pyr_exit$day]], year = rep(year_default, nrow(df))) # nocov
              } else { # nocov
                exit <- make_date(month = df[[pyr_exit$month]], year = rep(year_default, nrow(df)), day = rep(day_default, nrow(df))) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_exit)) { # nocov
                exit <- make_date(day = df[[pyr_exit$day]], year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
              } else {
                stop("Error: person-year exit missing day, month, and year")
              }
            }
          }
          if ("entry" %in% names(pyr)) { # format the entry as date column
            pyr_entry <- pyr$entry
            names(pyr_entry) <- tolower(names(pyr_entry)) # set the names to lowercase
            interval <- "interval"
            if ("year" %in% names(pyr_entry)) { # nocov
              if ("month" %in% names(pyr_entry)) { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
                } else { # nocov
                  entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = rep(day_default, nrow(df))) # nocov
                }
              } else { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(year = df[[pyr_entry$year]], day = df[[pyr_entry$day]], month = rep(month_default, nrow(df))) # nocov
                } else { # nocov
                  entry <- make_date(year = df[[pyr_entry$year]], day = rep(day_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
                }
              }
            } else { # nocov
              if ("month" %in% names(pyr_entry)) { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(month = df[[pyr_entry$month]], day = df[[pyr_entry$day]], year = rep(year_default, nrow(df))) # nocov
                } else { # nocov
                  entry <- make_date(month = df[[pyr_entry$month]], day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df))) # nocov
                }
              } else { # nocov
                if ("day" %in% names(pyr_entry)) { # nocov
                  entry <- make_date(day = df[[pyr_entry$day]], month = rep(month_default, nrow(df)), year = rep(year_default, nrow(df))) # nocov
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
          names(pyr_entry) <- tolower(names(pyr_entry)) # set the names to lowercase
          if ("year" %in% names(pyr_entry)) { # nocov
            if ("month" %in% names(pyr_entry)) { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
              } else { # nocov
                entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = rep(day_default, nrow(df))) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(year = df[[pyr_entry$year]], day = df[[pyr_entry$day]], month = rep(month_default, nrow(df))) # nocov
              } else { # nocov
                entry <- make_date(year = df[[pyr_entry$year]], day = rep(day_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
              }
            }
          } else { # nocov
            if ("month" %in% names(pyr_entry)) { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(month = df[[pyr_entry$month]], day = df[[pyr_entry$day]], year = rep(year_default, nrow(df))) # nocov
              } else { # nocov
                entry <- make_date(month = df[[pyr_entry$month]], day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df))) # nocov
              }
            } else { # nocov
              if ("day" %in% names(pyr_entry)) { # nocov
                entry <- make_date(day = df[[pyr_entry$day]], month = rep(month_default, nrow(df)), year = rep(year_default, nrow(df))) # nocov
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
        #
        df[["entry"]] <- entry
        df[["exit"]] <- exit
        df[["interval_dur"]] <- as.duration(interval(entry, exit))
        df <- df |> filter(interval_dur >= days(0))
        entry <- df$entry
        exit <- df$exit
#        kept <- interval_dur >= days(0)
#        df <- df[kept]
        #
        df_added <- tibble()
        pyr_unit <- "years" # default person-years to years
        if ("unit" %in% names(pyr)) {
          pyr_unit <- pyr$unit
        }
        for (time_i in 1:num_categ) { # for every time interval
          istart <- make_date(year = year[time_i], month = month[time_i], day = day[time_i]) # interval start
          iend <- make_date(year = year[time_i + 1], month = month[time_i + 1], day = day[time_i + 1]) # interval end
          # We don't want the upper limit to be inclusive, so we call roll it back 1 day
          # Start by checking if the category is one day
          bin_dur <- as.duration(interval(istart, iend))
          # The lowest unit of time is days, so we want to adjust any interval longer than 0 days down one
          if ((bin_dur > days(1)) && (time_i < num_categ)) {
            iend <- iend - days(1)
          }
          # Now the interval is only fully inclusive if the width is zero
          # The end of the categorical interval is counted, so if it is within, then we add the extra day
          cat_str <- paste(cat_str, paste("[", istart, " to ", iend, "]", sep = ""), sep = " ") # prepare the interval info
          categ_interval <- interval(istart, iend) # define as date interval
          c_categ <- list()
          if (interval == "left trunc") {
            # only exit
            a_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # exit within interval, set to i
            for (evt in event_cols) {
              c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # exit within interval set to event column value
            }
            b_categ <- case_when(a_categ == as.character(time_i) ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = 0) # for every row with exit in interval, track the duration from interval start to row end
            risk_interval <- interval(iend, exit) # interval from interval end to row end
            a_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # track rows which end after the interval, set to i
            b_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(categ_interval) + ddays(1), pyr_unit), .default = b_categ) # rows which end after interval are at risk the full interval
          } else if (interval == "right trunc") {
            # only entry
            a_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") #  start during the interval, set to i
            for (evt in event_cols) {
              c_categ[[evt]] <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # start during interval, set to event value
            }
            b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(interval(entry, iend)) + ddays(1), pyr_unit), .default = 0) # every row which starts during the interval, tracks duration from entry to interval end
            risk_interval <- interval(entry, istart)
            a_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # tracks which rows start before the interval
            b_categ <- case_when(as.numeric(as.duration(risk_interval)) > 0 & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(categ_interval) + ddays(1), pyr_unit), .default = b_categ) # rows which start before the interval are at risk the full interval
          } else {
            # both entry and exit
            risk_interval <- interval(entry, exit)
            a_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # category fully contained, set to i
            a_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval ends during row interval, set to i
            a_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval starts during row interval, set to i
            a_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # row interval fully contained in category interval, set to i
            for (evt in event_cols) {
              c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # row ends during category interval, set to event value
            }
            b_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(risk_interval), pyr_unit), .default = -1) # row interval fully in category interval, track full row interval
            b_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = b_categ) # rows which end during the category interval, track category interval start to row end
            b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(entry, iend)) + ddays(1), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
            b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(categ_interval) + ddays(1), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
            b_categ <- case_when( b_categ == -1 ~ 0.0, .default = b_categ) # set every unused interval to 0
          }
          index_kept <- seq_len(nrow(df))
          index_kept <- index_kept[a_categ == time_i] # indexes which contain the category interval to some level
          b_categ <- b_categ[a_categ == time_i] # durations for kept indexes
          row_kept <- slice(df, index_kept) # dataframe at kept rows
          row_kept <- row_kept %>% mutate("{cat_col}" := as.character(time_i)) # time category value
          row_kept <- row_kept %>% mutate("PYR" := b_categ) # duration value
          for (evt in event_cols) {
            d_categ <- c_categ[[evt]]
            d_categ <- d_categ[a_categ == time_i]
            row_kept <- row_kept %>% mutate("{evt}" := d_categ) # event values
          }
          # We don't want to keep rows with negative durations
          row_kept <- row_kept %>% filter("PYR" >= 0)
          df_added <- bind_rows(df_added, row_kept) # new updates dataset
        }
        df <- df_added
      }
      categ_bounds[[cat_col]] <- cat_str # update the list of category boundaries
    } else { # boundary as string, not a time category
      if (!cat_df %in% table_names) {
        stop("Error: ", cat_df, " not in table")
      }
      if (tolower(categ[[cat]]) == "factor") {
        # It is a factor
        cat_col <- paste(cat, "category", sep = "_")
        categ_cols <- c(categ_cols, cat_col)
        df[[cat_col]] <- as.factor(df[[cat]])
        cat_str <- paste(levels(df[[cat_col]]), collapse = "/")
        categ_bounds[[cat_col]] <- cat_str # add to list of intervals
      } else {
        cat_str <- ""
        temp <- categ[[cat]]
        temp <- gsub("/", " / ", temp, fixed = TRUE)
        temp <- gsub("]", " ] ", temp, fixed = TRUE)
        temp <- strsplit(temp, "\\s+")[[1]] # seperate values and delimiters
        match_index <- which(temp %in% c("/", "]"))
        if (length(match_index) < 1) {
          stop(paste("Error: Category ", cat, " did not have categories.", sep = ""))
        }
        categ_cols <- c(categ_cols, cat_col)
        df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize column
        match_i <- 1
        for (time_i in match_index) {
          # We want to figure out if there is a label
          L <- as.numeric(temp[time_i - 1])
          if (time_i == match_index[length(match_index)]) {
            # We are in the last entry
            if (length(temp) - time_i == 2) {
              # There is a label
              entry_label <- temp[time_i + 1]
              U <- as.numeric(temp[time_i + 2])
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              entry_label <- paste(L, U, sep = " - ")
            }
          } else {
            # Not at last entry
            if (match_index[match_i + 1] - time_i == 3) {
              # There is a label
              entry_label <- temp[time_i + 1]
              U <- as.numeric(temp[time_i + 2])
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              entry_label <- paste(L, U, sep = " - ")
            }
          }
          if (L == U) { # discrete case
            a_categ <- case_when(df[[cat_df]] == U & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          } else if (temp[time_i] == "/") { # strictly below upper bound
            a_categ <- case_when(df[[cat_df]] < U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, ")", sep = ""), sep = " ")
          } else { # including both bounds
            a_categ <- case_when(df[[cat_df]] <= U & df[[cat_df]] >= L & df[[cat_col]] == "Unassigned" ~ entry_label, .default = df[[cat_col]])
            cat_str <- paste(cat_str, paste("[", L, ", ", U, "]", sep = ""), sep = " ")
          }
          df[[cat_col]] <- a_categ # update tibble
          match_i <- match_i + 1
        }
        categ_bounds[[cat_col]] <- cat_str
      }
    }
  }
  if ("PYR" %in% names(df)) {
    # good, we need a person-years measure
  } else {
    if ("exit" %in% names(pyr)) { # format the exit as date column
      pyr_exit <- pyr$exit
      names(pyr_exit) <- tolower(names(pyr_exit)) # set the names to lowercase
      if ("year" %in% names(pyr_exit)) { # nocov
        if ("month" %in% names(pyr_exit)) { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]]) # nocov
          }
        } else { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(year = df[[pyr_exit$year]], day = df[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(year = df[[pyr_exit$year]]) # nocov
          }
        }
      } else {
        if ("month" %in% names(pyr_exit)) { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
          } else { # nocov
            exit <- make_date(month = df[[pyr_exit$month]]) # nocov
          }
        } else { # nocov
          if ("day" %in% names(pyr_exit)) { # nocov
            exit <- make_date(day = df[[pyr_exit$day]]) # nocov
          } else {
            stop("Error: person-year exit missing day, month, and year")
          }
        }
      }
      if ("entry" %in% names(pyr)) { # format the entry as date column
        pyr_entry <- pyr$entry
        names(pyr_entry) <- tolower(names(pyr_entry)) # set the names to lowercase
        interval <- "interval"
        if ("year" %in% names(pyr_entry)) { # nocov
          if ("month" %in% names(pyr_entry)) { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]]) # nocov
            }
          } else { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(year = df[[pyr_entry$year]], day = df[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(year = df[[pyr_entry$year]]) # nocov
            }
          }
        } else { # nocov
          if ("month" %in% names(pyr_entry)) { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
            } else { # nocov
              entry <- make_date(month = df[[pyr_entry$month]]) # nocov
            }
          } else { # nocov
            if ("day" %in% names(pyr_entry)) { # nocov
              entry <- make_date(day = df[[pyr_entry$day]]) # nocov
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
    summarize(AT_RISK = n(), PYR = sum(.data[["PYR"]]), .groups = "drop") # group by categories and define the durations and counts
  for (evt in names(summaries)) { # for each event summary
    if (grepl(" AS ", summaries[[evt]], fixed = TRUE)) { # get the method and column name
      temp <- gsub(" AS ", " ", summaries[[evt]], fixed = TRUE)
      temp <- strsplit(temp, "\\s+")[[1]]
      col_name <- temp[2]
      method <- temp[1]
    } else {
      col_name <- evt
      temp <- strsplit(summaries[[evt]], "\\s+")[[1]]
      method <- temp[[length(temp)]]
    }
    method <- tolower(method)
    method <- lapply(method, function(x) tryCatch(match.arg(x, choices = c("count", "sum", "mean", "rsum", "rmean", "weighted_mean")), error = function(error_message) x))[[1]] # match against expected values
    if (method %in% c("count", "sum", "rsum")) { # sum of event across each category combination
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method %in% c("mean", "rmean")) { # mean value of event across each category combination
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

#' Prints regression output cleanly
#'
#' \code{general_print} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @noRd
#' @return return nothing, prints the results to console
#' @family Output and Information Functions
general_print <- function(x, ...) {
  exargs <- list(...)
  digits <- 3
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

#' Prints fma regression output cleanly
#'
#' \code{general_fma_print} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @noRd
#' @return return nothing, prints the results to console
#' @family Output and Information Functions
general_fma_print <- function(x, ...) {
  exargs <- list(...)
  digits <- 3
  if ("digits" %in% names(exargs)) {
    digits <- exargs$digits
  } else if (length(exargs) == 1) {
    if (is.numeric(exargs[[1]])) {
      if (as.integer(exargs[[1]]) == exargs[[1]]) {
        digits <- exargs[[1]]
      }
    }
  }
  Interpret_FMA_Output(x, digits)
}

#' Prints a cox regression output clearly
#'
#' \code{print.coxres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.coxres <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a poisson regression output clearly
#'
#' \code{print.poisres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.poisres <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a case-control regression output clearly
#'
#' \code{print.caseconres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class caseconres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.caseconres <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a logistic regression output clearly
#'
#' \code{print.logitres} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class logitres
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.logitres <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a cox likelihood boundary regression output clearly
#'
#' \code{print.coxresbound} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxresbound
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.coxresbound <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a poisson likelihood boundary regression output clearly
#'
#' \code{print.poisresbound} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisresbound
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.poisresbound <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a logistic likelihood boundary regression output clearly
#'
#' \code{print.logitresbound} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class logitresbound
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.logitresbound <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a cox MCML regression output clearly
#'
#' \code{print.coxresmcml} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxresmcml
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.coxresmcml <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a poisson MCML regression output clearly
#'
#' \code{print.poisresmcml} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisresmcml
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.poisresmcml <- function(x, ...) {
  general_print(x, ...)
}

#' Prints a cox FMA regression output clearly
#'
#' \code{print.coxresfma} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class coxresmcml
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.coxresfma <- function(x, ...) {
  general_fma_print(x, ...)
}

#' Prints a poisson FMA regression output clearly
#'
#' \code{print.poisresfma} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @param x result object from a regression, class poisresmcml
#' @param ... can include the number of digits, named digit, or an unnamed integer entry assumed to be digits
#'
#' @return return nothing, prints the results to console
#' @noRd
#' @export
#' @family Output and Information Functions
print.poisresfma <- function(x, ...) {
  general_fma_print(x, ...)
}

#' Prints a regression output clearly
#'
#' \code{Interpret_Output} uses the list output from a regression, prints off a table of results and summarizes the score and convergence.
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return return nothing, prints the results to console
Interpret_Output <- function(out_list, digits = 3) {
  # make sure the output isn't an error
  passed <- out_list$Status
  message("|", paste(rep("-", options()$width), collapse = ""), "|")
  if (!is.na(passed)) {
    if ("Likelihood_Goal" %in% names(out_list)) {
      # likelihood boundary output
      model <- out_list$model
      modelcontrol <- out_list$modelcontrol
      para_number <- modelcontrol$para_number
      strata <- model$strata
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
      if (is(out_list, "coxresbound")) {
        message("Proportional Hazards Model")
      } else if (is(out_list, "poisresbound")) {
        message("Poisson Model")
      }
      if (all(strata != "NONE")) {
        message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
      }
      message(paste("Solving for the boundary of element: ", para_number, "\nApplied to column: '", name, "'\nSubterm: ", tform, "\nTerm number: ", term_n, sep = ""))
      if (neg[1]) {
        message(paste("Lower limit was not found, last step was at ", format(limits[1], digits = digits), " at a score of ", round(lik_bound[1], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
      } else {
        if (conv[1]) {
          message(paste("Lower limit converged to at ", format(limits[1], digits = digits), " at a score of ", round(lik_bound[1], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
        } else {
          message(paste("Lower limit reached ", format(limits[1], digits = digits), " at a score of ", round(lik_bound[1], digits), " with of goal of ", round(lik_goal, digits), " but did not converge", sep = ""))
        }
      }
      message(paste("Central estimate was ", format(beta_0, digits = digits), sep = ""))
      if (neg[2]) {
        message(paste("Upper limit was not found, last step was at ", format(limits[2], digits = digits), " at a score of ", round(lik_bound[2], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
      } else {
        if (conv[2]) {
          message(paste("Upper limit converged to at ", format(limits[2], digits = digits), " at a score of ", round(lik_bound[2], digits), " with of goal of ", round(lik_goal, digits), sep = ""))
        } else {
          message(paste("Upper limit reached ", format(limits[2], digits = digits), " at a score of ", round(lik_bound[2], digits), " with of goal of ", round(lik_goal, digits), " but did not converge", sep = ""))
        }
      }
    } else {
      # Check if its a multidose problem
      #      if (out_list$Survival_Type == "Cox_Multidose") {
      #        message("Currently the multiple realization code is not setup for printing results, due to the potentially large number of realizations")
      if (is(out_list, "caseconres")) {
        # case control output
        # get the model details
        null_model <- out_list$modelcontrol$null
        strata_odds <- out_list$StrataOdds
        KeptRecords <- out_list$UsedRecords
        RemovedRecords <- out_list$RejectedRecords
        if (!null_model) {
          names <- out_list$Parameter_Lists$names
          tforms <- out_list$Parameter_Lists$tforms
          term_n <- out_list$Parameter_Lists$term_n
          beta_0 <- out_list$beta_0
          keep_constant <- out_list$Parameter_Lists$keep_constant == 1
          if ("Standard_Error" %in% names(out_list)) {
            stdev <- out_list$Standard_Error
            pval <- 2 * pnorm(-abs(beta_0 / stdev))
            CI_low <- as.numeric(format(beta_0 - 1.96*stdev, digits = digits))
            CI_high <- as.numeric(format(beta_0 + 1.96*stdev, digits = digits))
            CI <- paste0("(",CI_low," - ", CI_high,")")
            res_table <- data.table(
              Covariate = names,
              Subterm = tforms,
              `Term Number` = term_n,
              Constant = keep_constant,
              `Central Estimate` = as.numeric(format(beta_0, digits = digits)),
              `Standard Error` = as.numeric(format(stdev, digits = digits)),
              `95% Confidence Interval` = CI,
              `2-tail p-value` = as.numeric(format(pval, digits = digits))
            )
          } else {
            res_table <- data.table(
              Covariate = names,
              Subterm = tforms,
              `Term Number` = term_n,
              Constant = keep_constant,
              `Central Estimate` = as.numeric(format(beta_0, digits = digits))
            )
          }
          if (!any(keep_constant)) {
            res_table <- res_table[, names(res_table)[names(res_table) != "Constant"], with = FALSE]
          }
          if (min(term_n) == max(term_n)) {
            res_table <- res_table[, names(res_table)[names(res_table) != "Term Number"], with = FALSE]
          }
        }
        deviance <- out_list$Deviance
        iteration <- out_list$Control_List$Iteration
        step_max <- out_list$Control_List$`Maximum Step`
        deriv_max <- out_list$Control_List$`Derivative Limiting`
        delta_ll <- out_list$Control_List$Delta_LogLik
        converged <- out_list$Converged
        #
        iter_lim <- out_list$control$maxiter
        ll_lim <- out_list$control$ll_epsilon
        step_lim <- out_list$control$epsilon
        deriv_lim <- out_list$control$deriv_epsilon
        #
        freepara <- out_list$FreeParameters
        freestrata <- out_list$FreeSets
        strata <- out_list$model$strata
        time_model <- out_list$modelcontrol$time_risk
        #
        modelform <- out_list$model$modelform
        form_type <- case_when(
          modelform == "M" ~ "Multiplicative Model Used: T0*T1*T2*...",
          modelform == "ME" ~ "Multiplicative-Excess Model Used: T0*(1+T1)*(1+T2)*...",
          modelform == "A" ~ "Additive Model Used: T0+T1+T2+...",
          modelform == "PA" ~ "Product-Additive Model Used: T0*(T1+T2+...)",
          modelform == "PAE" ~ "Product-Additive-Excess Model Used: T0*(1+T1+T2+...)",
          modelform == "GMIX" ~ "Geometric-Mixture Model Used: T0 *((1+T1)*(1+T2)*...)^(t)*(1+T1+T2+...)^(1-t)",
          .default = "Unknown"
        )
        #
        message("Final Results")
        if (null_model) {
          message("Null model used")
        } else {
          print(res_table)
        }
        #
        message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
        message("\nMatched Case-Control Model Used")
        if ((!null_model) && (min(term_n) != max(term_n))) {
          message(form_type)
        }
        if (all(strata != "NONE")) {
          if (time_model) {
            message("Model stratified by ", paste(shQuote(strata), " and time at risk", collapse = ", "))
          } else {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
          }
        } else if (time_model) {
          message("Model stratified by time at risk")
        } else {
          message("No risk grouping applied")
        }
        message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
        message(paste("Deviance: ", round(deviance, digits), sep = ""))
        message(paste(freestrata, " out of ", length(strata_odds), " matched sets used Unconditional Likelihood", sep = ""))
        if (!is.null(converged)) {
          if (iteration == 0) {
            message(paste("Iterations run: ", iteration, "\nmaximum step size: None taken, maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
          } else {
            message(paste("Iterations run: ", iteration, "\nmaximum step size: ", formatC(step_max, format = "e", digits = digits), ", maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
            if (delta_ll > 0) {
              message(paste("Last iteration improved the log-likelihood by: ", formatC(delta_ll, format = "e", digits = digits), sep = ""))
            } else {
              message("Log-likelihood was not improved in last iteration")
            }
          }
          if (converged) {
            message("Analysis converged")
          } else {
            if (iteration >= iter_lim) {
              message("Analysis did not converge, iteration limit was hit. Regression may converge with additional iterations ('maxiter').")
            } else if (step_max <= step_lim) {
              message("Analysis did not converge, step size limit was hit. Regression may converge if limit is reduced ('epsilon').")
            } else {
              message("Analysis did not converge.")
            }
          }
          neg_lim <- out_list$Control_List$`Ended on Negative Limit`
          if (neg_lim) {
            message("Warning: The regression ended after hitting a negative risk.")
          }
        }
        message("Records Used: ", KeptRecords, ", Records Removed: ", RemovedRecords)
      } else {
        # get the model details
        KeptRecords <- out_list$UsedRecords
        RemovedRecords <- out_list$RejectedRecords
        null_model <- out_list$modelcontrol$null
        if (!null_model) {
          ##
          names <- out_list$Parameter_Lists$names
          tforms <- out_list$Parameter_Lists$tforms
          term_n <- out_list$Parameter_Lists$term_n
          beta_0 <- out_list$beta_0
          keep_constant <- out_list$Parameter_Lists$keep_constant == 1
          if ("Standard_Error" %in% names(out_list)) {
            stdev <- out_list$Standard_Error
            pval <- 2 * pnorm(-abs(beta_0 / stdev))
            CI_low <- beta_0 - 1.96*stdev
            CI_high <- beta_0 + 1.96*stdev
            CI <- paste0("(",CI_low," - ", CI_high,")")
            res_table <- data.table(
              Covariate = names,
              Subterm = tforms,
              `Term Number` = term_n,
              Constant = keep_constant,
              `Central Estimate` = as.numeric(format(beta_0, digits = digits)),
              `Standard Error` = as.numeric(format(stdev, digits = digits)),
              `95% Confidence Interval` = CI,
              `2-tail p-value` = as.numeric(format(pval, digits = digits))
            )
          } else {
            res_table <- data.table(
              Covariate = names,
              Subterm = tforms,
              `Term Number` = term_n,
              Constant = keep_constant,
              `Central Estimate` = as.numeric(format(beta_0, digits = digits))
            )
          }
          if (!any(keep_constant)) {
            res_table <- res_table[, names(res_table)[names(res_table) != "Constant"], with = FALSE]
          }
          if (min(term_n) == max(term_n)) {
            res_table <- res_table[, names(res_table)[names(res_table) != "Term Number"], with = FALSE]
          }
          #
          modelform <- out_list$model$modelform
          form_type <- case_when(
            modelform == "M" ~ "Multiplicative Model Used: T0*T1*T2*...",
            modelform == "ME" ~ "Multiplicative-Excess Model Used: T0*(1+T1)*(1+T2)*...",
            modelform == "A" ~ "Additive Model Used: T0+T1+T2+...",
            modelform == "PA" ~ "Product-Additive Model Used: T0*(T1+T2+...)",
            modelform == "PAE" ~ "Product-Additive-Excess Model Used: T0*(1+T1+T2+...)",
            modelform == "GMIX" ~ "Geometric-Mixture Model Used: T0 *((1+T1)*(1+T2)*...)^(t)*(1+T1+T2+...)^(1-t)",
            .default = "Unknown"
          )
          #
        }
        message("Final Results")
        if (null_model) {
          message("Null model used")
        } else {
          print(res_table)
        }
        # get the model results
        LogLik <- out_list$LogLik
        AIC <- out_list$AIC
        BIC <- out_list$BIC
        deviation <- out_list$Deviance
        iteration <- out_list$Control_List$Iteration
        step_max <- out_list$Control_List$`Maximum Step`
        deriv_max <- out_list$Control_List$`Derivative Limiting`
        delta_ll <- out_list$Control_List$Delta_LogLik
        strata <- out_list$model$strata
        strata_level <- out_list$strata_levels
        cens_weight <- out_list$model$weight
        converged <- out_list$Converged
        #
        iter_lim <- out_list$control$maxiter
        ll_lim <- out_list$control$ll_epsilon
        step_lim <- out_list$control$epsilon
        deriv_lim <- out_list$control$deriv_epsilon
        #
        message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
        if (is(out_list, "coxres")) {
          if (cens_weight == "NONE") {
            # cox model
            message("\nCox Model Used")
          } else {
            # fine-gray model
            message(paste("\nFine-Gray Model Used, weighted by ", cens_weight, sep = ""))
          }
          #
          tstart <- out_list$model$start_age
          tend <- out_list$model$end_age
          event <- out_list$model$event
          if (tstart == "right_trunc") {
            message("Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
          } else if (tend == "left_trunc") {
            message("Entry Age Column was: '", tstart, "', Outcome Column was: '", event, "'")
          } else {
            message("Entry Age Column was: '", tstart, "', Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
          }
          if (cens_weight != "NONE") {
            message("Survival Weighting Column was :'", cens_weight, "'")
          }
          #
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
          }
          risk_groups <- out_list$RiskGroups
          message("Risk Groups Used: ", risk_groups)
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  AIC: ", round(AIC, digits), sep = ""))
        } else if (is(out_list, "coxresmcml")) {
          if (cens_weight == "NONE") {
            # cox model
            message("\nCox Model Used")
          } else {
            # fine-gray model
            message(paste("\nFine-Gray Model Used, weighted by ", cens_weight, sep = ""))
          }
          #
          tstart <- out_list$model$start_age
          tend <- out_list$model$end_age
          event <- out_list$model$event
          if (tstart == "right_trunc") {
            message("Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
          } else if (tend == "left_trunc") {
            message("Entry Age Column was: '", tstart, "', Outcome Column was: '", event, "'")
          } else {
            message("Entry Age Column was: '", tstart, "', Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
          }
          if (cens_weight != "NONE") {
            message("Survival Weighting Column was :'", cens_weight, "'")
          }
          #
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
          }
          risk_groups <- out_list$RiskGroups
          message("Risk Groups Used: ", risk_groups)
          realizations <- out_list$realizations
          message("Realizations Used: ", realizations)
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  AIC: ", round(AIC, digits), sep = ""))
        } else if (is(out_list, "poisres")) {
          # poisson model
          message("\nPoisson Model Used")
          pyr_col <- out_list$model$person_year
          evt_col <- out_list$model$event
          message("Person-year Column: '", pyr_col, "'")
          message("Event Column: '", evt_col, "'")
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
            message("Strata split into ", strata_level, " distinct levels", sep = "")
          }
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  Deviance: ", round(deviation, digits), ",  AIC: ", round(AIC, digits), ",  BIC: ", round(BIC, digits), sep = ""))
        } else if (is(out_list, "poisresmcml")) {
          # poisson model
          message("\nPoisson Model Used")
          pyr_col <- out_list$model$person_year
          evt_col <- out_list$model$event
          message("Person-year Column: '", pyr_col, "'")
          message("Event Column: '", evt_col, "'")
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
            message("Strata split into ", strata_level, " distinct levels", sep = "")
          }
          realizations <- out_list$realizations
          message("Realizations Used: ", realizations)
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  Deviance: ", round(deviation, digits), ",  AIC: ", round(AIC, digits), ",  BIC: ", round(BIC, digits), sep = ""))
        } else if (is(out_list, "logitres")) {
          # logistic model
          message("\nLogisitic Model Used")
          odds <- out_list$modelcontrol$logit_odds
          link <- "Unknown"
          if (odds) {
            link <- "Odds Ratio"
          } else {
            ident <- out_list$modelcontrol$logit_ident
            if (ident) {
              link <- "Identity"
            } else {
              loglink <- out_list$modelcontrol$logit_loglink
              if (loglink) {
                link <- "Complementary Log"
              } else {
                probit <- out_list$modelcontrol$logit_probit
                if (probit) {
                  link <- "Probability Unit (probit)"
                }
              }
            }
          }
          message(link, " Linking Function Used")
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
          }
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  Deviance: ", round(deviation, digits), ",  AIC: ", round(AIC, digits), ",  BIC: ", round(BIC, digits), sep = ""))
        } else {
          message("\nUnknown Model Used")
          if ((!null_model) && (min(term_n) != max(term_n))) {
            message(form_type)
          }
          if (all(strata != "NONE")) {
            message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
          }
          message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
          message(paste("-2*Log-Likelihood: ", round(-2 * LogLik, digits), ",  AIC: ", round(AIC, digits), sep = ""))
        }
        if (!is.null(converged)) {
          if (iteration == 0) {
            message(paste("Iterations run: ", iteration, "\nmaximum step size: None taken, maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
          } else {
            message(paste("Iterations run: ", iteration, "\nmaximum step size: ", formatC(step_max, format = "e", digits = digits), ", maximum first derivative: ", formatC(deriv_max, format = "e", digits = digits), sep = ""))
            if (delta_ll > 0) {
              message(paste("Last iteration improved the log-likelihood by: ", formatC(delta_ll, format = "e", digits = digits), sep = ""))
            } else {
              message("Log-likelihood was not improved in last iteration")
            }
          }
          if (converged) {
            message("Analysis converged")
          } else {
            if (iteration >= iter_lim) {
              message("Analysis did not converge, iteration limit was hit. Regression may converge with additional iterations ('maxiter').")
            } else if (step_max <= step_lim) {
              message("Analysis did not converge, step size limit was hit. Regression may converge if limit is reduced ('epsilon').")
            } else {
              message("Analysis did not converge.")
            }
          }
          neg_lim <- out_list$Control_List$`Ended on Negative Limit`
          if (neg_lim) {
            message("Warning: The last iteration encountered a negative risk.")
          }
        }
        message("Records Used: ", KeptRecords, ", Records Removed: ", RemovedRecords)
      }
    }
  } else {
    message(paste("Regression Failed"))
  }
  if ("RunTime" %in% names(out_list)) {
    run_time_sec <- as.numeric(out_list$RunTime, units = "secs")
    run_time_min <- as.numeric(out_list$RunTime, units = "mins")
    run_time_hour <- as.numeric(out_list$RunTime, units = "hours")
    # nocov start
    if (run_time_sec < 60) {
      message(paste("Run finished in ", round(run_time_sec, digits), " seconds", sep = ""))
    } else if (run_time_min < 60) {
      message(paste("Run finished in ", round(run_time_min, digits), " minutes", sep = ""))
    } else {
      message(paste("Run finished in ", round(run_time_hour, digits), " hours", sep = ""))
    }
    # nocov end
  }
  message("|", paste(rep("-", options()$width), collapse = ""), "|")
}

#' Prints a FMA regression output clearly
#'
#' \code{Interpret_FMA_Output} uses the list output from a FMA regression, summarizes the model and some results.
#'
#' @inheritParams R_template
#'
#' @noRd
#' @return return nothing, prints the results to console
Interpret_FMA_Output <- function(out_list, digits = 3) {
  # make sure the output isn't an error
  passed <- out_list$Status
  message("|", paste(rep("-", options()$width), collapse = ""), "|")

  # get the model details
  KeptRecords <- out_list$UsedRecords
  RemovedRecords <- out_list$RejectedRecords
  ##
  names <- out_list$Parameter_Lists$names
  tforms <- out_list$Parameter_Lists$tforms
  term_n <- out_list$Parameter_Lists$term_n
  keep_constant <- out_list$Parameter_Lists$keep_constant == 1
  res_table <- data.table(
    Covariate = names,
    Subterm = tforms,
    `Term Number` = term_n,
    Constant = keep_constant
  )
  if (!any(keep_constant)) {
    res_table <- res_table[, names(res_table)[names(res_table) != "Constant"], with = FALSE]
  }
  if (min(term_n) == max(term_n)) {
    res_table <- res_table[, names(res_table)[names(res_table) != "Term Number"], with = FALSE]
  }
  #
  modelform <- out_list$model$modelform
  form_type <- case_when(
    modelform == "M" ~ "Multiplicative Model Used: T0*T1*T2*...",
    modelform == "ME" ~ "Multiplicative-Excess Model Used: T0*(1+T1)*(1+T2)*...",
    modelform == "A" ~ "Additive Model Used: T0+T1+T2+...",
    modelform == "PA" ~ "Product-Additive Model Used: T0*(T1+T2+...)",
    modelform == "PAE" ~ "Product-Additive-Excess Model Used: T0*(1+T1+T2+...)",
    modelform == "GMIX" ~ "Geometric-Mixture Model Used: T0 *((1+T1)*(1+T2)*...)^(t)*(1+T1+T2+...)^(1-t)",
    .default = "Unknown"
  )
  #
  message("Final Serial Analysis Results")
  print(res_table)
  # get the model results
  LogLik <- out_list$LogLik
  AIC <- out_list$AIC
  BIC <- out_list$BIC
  deviation <- out_list$Deviance
  strata <- out_list$model$strata
  strata_level <- out_list$strata_levels
  cens_weight <- out_list$model$weight
  #
  message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
  if (is(out_list, "coxresfma")) {
    if (cens_weight == "NONE") {
      # cox model
      message("\nCox Model Used")
    } else {
      # fine-gray model
      message(paste("\nFine-Gray Model Used, weighted by ", cens_weight, sep = ""))
    }
    #
    tstart <- out_list$model$start_age
    tend <- out_list$model$end_age
    event <- out_list$model$event
    if (tstart == "right_trunc") {
      message("Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
    } else if (tend == "left_trunc") {
      message("Entry Age Column was: '", tstart, "', Outcome Column was: '", event, "'")
    } else {
      message("Entry Age Column was: '", tstart, "', Survival Age Column was: '", tend, "', Outcome Column was: '", event, "'")
    }
    if (cens_weight != "NONE") {
      message("Survival Weighting Column was :'", cens_weight, "'")
    }
    #
    if (min(term_n) != max(term_n)) {
      message(form_type)
    }
    if (all(strata != "NONE")) {
      message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
    }
    risk_groups <- out_list$RiskGroups
    message("Risk Groups Used: ", risk_groups)
    realizations <- out_list$realizations
    message("Realizations Used: ", realizations)
    message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
  } else if (is(out_list, "poisresfma")) {
    # poisson model
    message("\nPoisson Model Used")
    pyr_col <- out_list$model$person_year
    evt_col <- out_list$model$event
    message("Person-year Column: '", pyr_col, "'")
    message("Event Column: '", evt_col, "'")
    if (min(term_n) != max(term_n)) {
      message(form_type)
    }
    if (all(strata != "NONE")) {
      message("Model stratified by ", paste(shQuote(strata), collapse = ", "))
      message("Strata split into ", strata_level, " distinct levels", sep = "")
    }
    realizations <- out_list$realizations
    message("Realizations Used: ", realizations)
    message("|", paste(rep("-", as.integer(options()$width / 2)), collapse = " "), "|")
  } else {
    stop("\nUnknown Model Used")
  }
  if ("RunTime" %in% names(out_list)) {
    run_time_sec <- as.numeric(out_list$RunTime, units = "secs")
    run_time_min <- as.numeric(out_list$RunTime, units = "mins")
    run_time_hour <- as.numeric(out_list$RunTime, units = "hours")
    # nocov start
    if (run_time_sec < 60) {
      message(paste("Run finished in ", round(run_time_sec, digits), " seconds", sep = ""))
    } else if (run_time_min < 60) {
      message(paste("Run finished in ", round(run_time_min, digits), " minutes", sep = ""))
    } else {
      message(paste("Run finished in ", round(run_time_hour, digits), " hours", sep = ""))
    }
    # nocov end
  }
  message("|", paste(rep("-", options()$width), collapse = ""), "|")
  if (all(passed != "PASSED")) {
    message("Error: Every realization failed.")
  } else if (any(passed != "PASSED")) {
    message("Warning: Atleast one realization failed.")
  }
  message("Further results can be extracted using 'output$LogLik', 'output$Parameters', 'output$Standard_Error', etc.")
  message("|", paste(rep("-", options()$width), collapse = ""), "|")
}
