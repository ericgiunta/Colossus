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
#' @family Table Generation Functions
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
  table_names <- names(table)
  categ_cols <- c()
  categ_bounds <- list()
  #
  categ_res <- Category_Process(df, table_names, categ, categ_cols, categ_bounds)
  df <- categ_res$df
  categ_cols <- categ_res$categ_cols
  categ_bounds <- categ_res$categ_bounds
  #
  df_group <- generate_summaries(df, events, c(), c(), categ_cols, fcount = FALSE, lcount = FALSE, time_table = FALSE)
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
#' @param time_scale list with the time scale information, either a calendar category or an age category
#' @param categ list with category columns and methods, methods can be either strings or lists of boundaries
#' @param summaries list of columns to summarize, supports counts, means, and weighted means by person-year and renaming the summary column
#' @param events list of events or interests, checks if events are within each time interval
#' @param fcount boolean if the `first at risk` count should be returned in the data. Returns the number of observations starting in each combination of categories.
#' @param lcount boolean if the `last at risk` count should be returned in the data. Returns the number of observations ending in each combination of categories.
#' @param studyid id used to determine distinct subjects used for first and last at risk intervals
#'
#' @return returns a grouped table and a list of category boundaries used
#' @family Table Generation Functions
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
#' calendar_categ <- list(
#'   type = "calendar",
#'   day = c(1, 1, 1),
#'   month = c(1, 1, 1),
#'   year = c(1901, 1904, 1908)
#' )
#' time_scale <- list("time" = calendar_categ)
#' summary <- list(
#'   c = "count AS cases"
#' )
#' events <- list("c")
#' pyr <- list(
#'   entry = list(year = "f"),
#'   exit = list(year = "i"),
#'   unit = "years"
#' )
#' e <- Event_Time_Gen(table, pyr, time_scale, categ, summary, events)
#'
Event_Time_Gen <- function(table, pyr = list(), time_scale = list(), categ = list(), summaries = list(), events = c(), fcount = FALSE, lcount = FALSE, studyid = "studyID", verbose = FALSE) {
  df <- as_tibble(table)
  # setting commands and known columns to dummy values to avoid `no visible binding` warnings
  `%within%` <- lubridate::`%within%`
  `%>%` <- dplyr::`%>%`
  `%m+%` <- lubridate::`%m+%`
  interval_dur <- ""
  entry_age <- ""
  #
  if (nrow(df) < 1) {
    stop("Error: The table passed had no data.")
  }
  # check for any restricted names, which will be overwritten
  restricted_names <- c("F_AT_RISK", "L_AT_RISK", "day_entry", "month_entry", "year_entry", "day_exit", "month_exit", "year_exit", "birth", "entry", "exit", "entry_age", "exit_age", "interval_dur", "_arbitrary_age", "PYR")
  for (r_name in restricted_names) {
    if (r_name %in% names(df)) {
      warning(paste0("Warning: `", r_name, "` is a restricted name and may be overwritten."))
    }
  }
  #
  table_names <- names(table)
  # Setting default date values
  year_default <- 1900
  month_default <- 1
  day_default <- 1
  # We want to determine which type of timescale we are using
  # If there are day/month/year info then we are in calendar time-scale
  # We want to prepare extra columns for any filled in data
  names(pyr) <- tolower(names(pyr))
  if (!any(names(pyr) %in% c("entry", "exit"))) {
    stop("Error: person-year columns did not contain entry or exit values.")
  }
  if (missing(studyid)) {
    studyid <- "studyID"
    df[["studyID"]] <- seq_len(nrow(df))
  } else {
    if (!(studyid %in% names(df))) {
      stop(paste0("Error: Study id was set to `", studyid, "`, but it was not in the data."))
    }
  }
  # We start by assuming no time-scale
  time_scale_version <- "unknown"
  if ("exit" %in% names(pyr)) { # format the exit as date column
    if (inherits(pyr$exit, "list")) { # the entry is a list
      names(pyr$exit) <- tolower(names(pyr$exit)) # set the names to lowercase
      pyr$exit$trunc <- FALSE
      if (!any(names(pyr$exit) %in% c("day", "month", "year"))) {
        stop("Error: Interval exit list did not contain day, month, or year.")
      } else {
        time_scale_version <- "calendar"
        if (!("day" %in% names(pyr$exit))) {
          pyr$exit$day <- "day_exit"
          df$day_exit <- day_default
        }
        if (!("month" %in% names(pyr$exit))) {
          pyr$exit$month <- "month_exit"
          df$month_exit <- month_default
        }
        if (!("year" %in% names(pyr$exit))) {
          pyr$exit$year <- "year_exit"
          df$year_exit <- year_default
        }
      }
    } else {
      time_scale_version <- "user"
      pyr$exit <- list("duration" = pyr$exit, "trunc" = FALSE)
    }
  }
  # We cannot determine what it should be, check the entry for more information
  if ("entry" %in% names(pyr)) {
    if (inherits(pyr$entry, "list")) { # the entry is a list
      names(pyr$entry) <- tolower(names(pyr$entry)) # set the names to lowercase
      pyr$entry$trunc <- FALSE
      if (!any(names(pyr$entry) %in% c("day", "month", "year"))) {
        if (time_scale_version == "calendar") {
          stop("Error: Person-year exit used calendar time-scale, but the entry did not.")
        }
        stop("Error: Interval entry list did not contain day, month, or year.")
      } else {
        time_scale_version <- "calendar"
        if (!("day" %in% names(pyr$entry))) {
          pyr$entry$day <- "day_entry"
          df$day_entry <- day_default
        }
        if (!("month" %in% names(pyr$entry))) {
          pyr$entry$month <- "month_entry"
          df$month_entry <- month_default
        }
        if (!("year" %in% names(pyr$entry))) {
          pyr$entry$year <- "year_entry"
          df$year_entry <- year_default
        }
      }
    } else {
      time_scale_version <- "user"
      pyr$entry <- list("duration" = pyr$entry, "trunc" = FALSE)
    }
  }
  # Now we check if there was missing data
  if (!(time_scale_version %in% c("user", "calendar"))) {
    stop("Error: Time scale was not set.")
  } else if (time_scale_version == "calendar") {
    if (!("entry" %in% names(pyr))) {
      pyr$entry <- list(
        "day" = "day_entry",
        "month" = "month_entry",
        "year" = "year_entry"
      )
      pyr$entry$trunc <- TRUE
      df$day_entry <- day_default
      df$month_entry <- month_default
      df$year_entry <- year_default
    }
    if (!("exit" %in% names(pyr))) {
      pyr$exit <- list(
        "day" = "day_exit",
        "month" = "month_exit",
        "year" = "year_exit"
      )
      pyr$exit$trunc <- TRUE
      df$day_exit <- day_default
      df$month_exit <- month_default
      df$year_exit <- year_default
    }
  } else {
    if (!("entry" %in% names(pyr))) {
      pyr$entry <- list("duration" = "entry_age", "trunc" = TRUE)
      df[["entry_age"]] <- 0
    }
    if (!("exit" %in% names(pyr))) {
      pyr$exit <- list("duration" = "exit_age", "trunc" = TRUE)
      df[["exit_age"]] <- 0
    }
  }
  if ((pyr$entry$trunc == pyr$exit$trunc) && (pyr$entry$trunc == TRUE)) {
    stop("Error: Neither entry nor exit dates/durations were given. Atleast one needs to be provided for person-time tables.")
  }
  # Start the first at risk category
  # without any time categories, every value will contibute
  df[["F_AT_RISK"]] <- 1
  df[["L_AT_RISK"]] <- 1
  # Checking for errors or valid data
  if (length(events) == 0) {
    stop("Error: no events were given")
  }
  time_cols <- c()
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
  # Process the Calendar time scale
  if (time_scale_version == "calendar") {
    time_res <- Calendar_Process(df, table_names, pyr, time_scale, event_cols, studyid, verbose)
    df <- time_res$df
    categ_cols <- time_res$categ_cols
    categ_bounds <- time_res$categ_bounds
  } else {
    time_res <- UserScale_Process(df, table_names, pyr, time_scale, event_cols, studyid, verbose)
    df <- time_res$df
    categ_cols <- time_res$categ_cols
    categ_bounds <- time_res$categ_bounds
  }
  # Now the non-time categories
  categ_res <- Category_Process(df, table_names, categ, categ_cols, categ_bounds)
  df <- categ_res$df
  categ_cols <- categ_res$categ_cols
  categ_bounds <- categ_res$categ_bounds
  # check that required columns are not already in use
  if ("PYR" %in% names(summaries)) { # storing person-year durations
    stop("Error: 'PYR' listed as a event column, either remove or rename with ' AS '")
  }
  if ("AT_RISK" %in% names(summaries)) { # storing number at risk
    stop("Error: 'AT_RISK' listed as a event column, either remove or rename with ' AS '")
  }
  df_group <- generate_summaries(df, summaries, event_cols, event_names, categ_cols, fcount, lcount, time_table = TRUE)
  #
  list(df = as.data.table(df_group), bounds = categ_bounds) # , ungrouped = as.data.table(df))
}

#' uses a table and calendar time scale categories to split a table into risk intervals
#'
#' \code{Calendar_Process} generates table split by risk categories
#'
#' @inheritParams R_template
#' @param df dataframe with every category/event column needed
#' @param pyr list with entry and exit lists, containing day/month/year columns in the table
#' @param time_scale list with the time scale information, either a calendar category or an age category
#' @param event_cols list of event columns
#' @param studyid id used to determine distinct subjects used for first and last at risk intervals
#'
#' @return returns calendar based table generation
#' @family Table Generation Functions
#' @noRd
Calendar_Process <- function(df, table_names, pyr = list(), time_scale = list(), event_cols = c(), studyid = "studyID", verbose = FALSE) {
  # setting commands and known columns to dummy values to avoid `no visible binding` warnings
  `%within%` <- lubridate::`%within%`
  `%>%` <- dplyr::`%>%`
  `%m+%` <- lubridate::`%m+%`
  interval_dur <- ""
  entry_age <- ""
  #
  year_default <- 1900
  month_default <- 1
  day_default <- 1
  #
  pyr_entry <- pyr$entry
  pyr_exit <- pyr$exit
  y_entry <- pyr_entry$year
  m_entry <- pyr_entry$month
  d_entry <- pyr_entry$day
  y_exit <- pyr_exit$year
  m_exit <- pyr_exit$month
  d_exit <- pyr_exit$day
  #
  categ_cols <- c()
  categ_bounds <- list()
  for (cat in names(time_scale)) { # for each time category
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
    names(time_scale[[cat]]) <- tolower(names(time_scale[[cat]]))
    if (!("type" %in% names(time_scale[[cat]]))) { # if the type is not listed
      type <- "none"
      if (all(names(time_scale[[cat]]) %in% c("day", "month", "year"))) {
        type <- "calendar" # calendar scale only uses day/month/year
      } else if (any(names(time_scale[[cat]]) %in% c("day", "month", "year"))) {
        type <- "age" # age scale only use day/month/year for the birth and other values for the actual categories
      } else {
        stop("Error: Unused time category options. Currently only calendar and person-age are supported")
      }
      time_scale[[cat]]$type <- type
    } else if (!(time_scale[[cat]]$type %in% c("calendar", "age"))) {
      type <- "none"
      if (all(names(time_scale[[cat]]) %in% c("day", "month", "year"))) {
        type <- "calendar" # calendar scale only uses day/month/year
      } else if (any(names(time_scale[[cat]]) %in% c("day", "month", "year"))) {
        type <- "age" # age scale only use day/month/year for the birth and other values for the actual categories
      } else {
        stop("Error: Unused time category options. Currently only calendar and person-age are supported")
      }
      time_scale[[cat]]$type <- type
    }
    pyr_unit <- "years" # default person-years to years
    if ("unit" %in% names(time_scale[[cat]])) {
      pyr_unit <- time_scale[[cat]]$unit
    }
    if (time_scale[[cat]]$type == "calendar") {
      if ("day" %in% names(time_scale[[cat]])) { # determine the day, month, year data for the category
        day_categ <- time_scale[[cat]]$day # nocov
        if ("month" %in% names(time_scale[[cat]])) { # nocov
          month_categ <- time_scale[[cat]]$month # nocov
        } else { # nocov
          month_categ <- rep(month_default, length(day_categ)) # nocov
        }
        if ("year" %in% names(time_scale[[cat]])) { # nocov
          year_categ <- time_scale[[cat]]$year # nocov
        } else {
          year_categ <- rep(year_default, length(day_categ)) # nocov
        }
      } else if ("month" %in% names(time_scale[[cat]])) { # nocov
        month_categ <- time_scale[[cat]]$month # nocov
        day_categ <- rep(day_default, length(month_categ)) # nocov
        if ("year" %in% names(time_scale[[cat]])) { # nocov
          year_categ <- time_scale[[cat]]$year # nocov
        } else { # nocov
          year_categ <- rep(year_default, length(month_categ)) # nocov
        }
      } else if ("year" %in% names(time_scale[[cat]])) { # nocov
        year_categ <- time_scale[[cat]]$year # nocov
        day_categ <- rep(day_default, length(year_categ)) # nocov
        month_categ <- rep(month_default, length(year_categ)) # nocov
      } else {
        stop("Error: calendar category missing 'day', 'month', and 'year'") # all three are required
      }
      num_categ <- length(year_categ) - 1
      categ_cols <- c(categ_cols, cat_col)
      df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
      #
      interval <- "interval"
      pyr_entry <- pyr$entry
      pyr_exit <- pyr$exit
      entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
      exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
      if (pyr_entry$trunc) {
        entry <- exit - years(1)
      }
      if (pyr_exit$trunc) {
        exit <- entry + years(1)
      }
      if (any(is.na(entry))) {
        stop("Error: Entry date failed. Check dates for any impossible values.")
      }
      if (any(is.na(exit))) {
        stop("Error: Exit date failed. Check dates for any impossible values.")
      }
      #
      y_entry <- pyr_entry$year
      m_entry <- pyr_entry$month
      d_entry <- pyr_entry$day
      y_exit <- pyr_exit$year
      m_exit <- pyr_exit$month
      d_exit <- pyr_exit$day
      #
      df[["entry"]] <- entry
      df[["exit"]] <- exit
      df[["interval_dur"]] <- as.duration(interval(entry, exit))
      df <- df |> filter(interval_dur >= days(0))
      if (nrow(df) < 1) {
        stop("Error: The table passed had no data at risk after filtering for negative durations.")
      }
      entry <- df$entry
      exit <- df$exit
      df_added <- tibble()
      for (time_i in 1:num_categ) { # for every time interval
        istart <- make_date(year = year_categ[time_i], month = month_categ[time_i], day = day_categ[time_i]) # interval start
        iend <- make_date(year = year_categ[time_i + 1], month = month_categ[time_i + 1], day = day_categ[time_i + 1]) # interval end
        if (any(is.na(istart))) {
          stop(paste0("Error: Calender entry date failed (d/m/y: ", day_categ[time_i], "/", month_categ[time_i], "/", year_categ[time_i], "). Check for impossible value."))
        }
        if (any(is.na(iend))) {
          stop(paste0("Error: Calender exit date failed (d/m/y: ", day_categ[time_i + 1], "/", month_categ[time_i + 1], "/", year_categ[time_i + 1], "). Check for impossible value."))
        }
        # We don't want the upper limit to be inclusive, so we call roll it back 1 day
        # Start by checking if the category is one day
        bin_dur <- as.duration(interval(istart, iend))
        # The lowest unit of time is days, so we want to adjust any interval longer than 0 days down one, excluding the last time interval which includes the upper limit
        if ((bin_dur > days(1)) && (time_i < num_categ)) {
          iend <- iend - days(1)
        }
        # We also want to be sure the bin is increasing in time
        if (bin_dur < days(0)) {
          stop(paste0("Error: Calender category starting at (d/m/y: ", day_categ[time_i], "/", month_categ[time_i], "/", year_categ[time_i], ") was not increasing with time."))
        }
        if (pyr_entry$trunc) {
          entry <- istart - years(1)
        }
        if (pyr_exit$trunc) {
          exit <- iend + years(1)
        }
        # Now the interval is only fully inclusive if the width is zero
        # The end of the categorical interval is counted, so if it is within, then we add the extra day
        cat_str <- paste(cat_str, paste("[", istart, " to ", iend, "]", sep = ""), sep = " ") # prepare the interval info
        categ_interval <- interval(istart, iend) # define as date interval
        c_categ <- list()
        # both entry and exit
        risk_interval <- interval(entry, exit)
        # When we are not in final category, we correct the interval end time
        day_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::day(istart), .default = 0)
        month_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::month(istart), .default = 0)
        year_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::year(istart), .default = 0)
        #
        day_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::day(iend %m+% lubridate::days(1)), .default = 0)
        month_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::month(iend %m+% lubridate::days(1)), .default = 0)
        year_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::year(iend %m+% lubridate::days(1)), .default = 0)
        if (time_i == num_categ) {
          end_dif <- as.duration(interval(iend, exit)) # check if the risk and category intervals end at the same time
          #
          day_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::day(iend), .default = 0)
          month_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::month(iend), .default = 0)
          year_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::year(iend), .default = 0)
        }
        a_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # category fully contained, set to i
        #
        day_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(entry), .default = day_enter_categ)
        month_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(entry), .default = month_enter_categ)
        year_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(entry), .default = year_enter_categ)
        if ((bin_dur > days(1)) && (time_i < num_categ)) {
          day_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend %m+% lubridate::days(1)), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend %m+% lubridate::days(1)), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend %m+% lubridate::days(1)), .default = year_exit_categ)
        } else {
          end_dif <- as.duration(interval(iend, exit)) # check if the risk and category intervals end at the same time
          day_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend %m+% lubridate::days(1)), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend %m+% lubridate::days(1)), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend %m+% lubridate::days(1)), .default = year_exit_categ)
          #
          day_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend), .default = year_exit_categ)
        }
        a_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval ends during row interval, set to i
        #
        day_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(istart), .default = day_enter_categ)
        day_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(exit), .default = day_exit_categ)
        month_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(istart), .default = month_enter_categ)
        month_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(exit), .default = month_exit_categ)
        year_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(istart), .default = year_enter_categ)
        year_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(exit), .default = year_exit_categ)
        a_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval starts during row interval, set to i
        #
        day_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(entry), .default = day_enter_categ)
        day_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(exit), .default = day_exit_categ)
        month_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(entry), .default = month_enter_categ)
        month_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(exit), .default = month_exit_categ)
        year_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(entry), .default = year_enter_categ)
        year_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(exit), .default = year_exit_categ)
        a_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # row interval fully contained in category interval, set to i
        #
        for (evt in event_cols) {
          c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # row ends during category interval, set to event value
        }
        b_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(risk_interval), pyr_unit), .default = -1) # row interval fully in category interval, track full row interval
        b_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = b_categ) # rows which end during the category interval, track category interval start to row end
        if ((bin_dur > days(1)) && (time_i < num_categ)) {
          b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(entry, iend)) + ddays(1), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
          b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(categ_interval) + ddays(1), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
        } else {
          b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(entry, iend)), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
          b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(categ_interval), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
        }
        b_categ <- case_when(b_categ == -1 ~ 0.0, .default = b_categ) # set every unused interval to 0
        #
        index_kept <- seq_len(nrow(df))
        index_kept <- index_kept[a_categ == time_i] # indexes which contain the category interval to some level
        b_categ <- b_categ[a_categ == time_i]
        index_kept <- index_kept[b_categ >= 0]
        #
        day_enter_categ <- day_enter_categ[a_categ == time_i]
        day_exit_categ <- day_exit_categ[a_categ == time_i]
        month_enter_categ <- month_enter_categ[a_categ == time_i]
        month_exit_categ <- month_exit_categ[a_categ == time_i]
        year_enter_categ <- year_enter_categ[a_categ == time_i]
        year_exit_categ <- year_exit_categ[a_categ == time_i]
        #
        day_enter_categ <- day_enter_categ[b_categ >= 0]
        day_exit_categ <- day_exit_categ[b_categ >= 0]
        month_enter_categ <- month_enter_categ[b_categ >= 0]
        month_exit_categ <- month_exit_categ[b_categ >= 0]
        year_enter_categ <- year_enter_categ[b_categ >= 0]
        year_exit_categ <- year_exit_categ[b_categ >= 0]
        #
        row_kept <- slice(df, index_kept) # dataframe at kept rows
        row_kept <- row_kept %>% mutate("{cat_col}" := as.character(time_i)) # time category value
        #
        row_kept <- row_kept %>% mutate("{y_entry}" := year_enter_categ)
        row_kept <- row_kept %>% mutate("{m_entry}" := month_enter_categ)
        row_kept <- row_kept %>% mutate("{d_entry}" := day_enter_categ)
        row_kept <- row_kept %>% mutate("{y_exit}" := year_exit_categ)
        row_kept <- row_kept %>% mutate("{m_exit}" := month_exit_categ)
        row_kept <- row_kept %>% mutate("{d_exit}" := day_exit_categ)
        #
        for (evt in event_cols) {
          d_categ <- c_categ[[evt]]
          d_categ <- d_categ[a_categ == time_i]
          d_categ <- d_categ[b_categ >= 0]
          row_kept <- row_kept %>% mutate("{evt}" := d_categ) # event values
        }
        df_added <- bind_rows(df_added, row_kept) # new updates dataset
      }
      categ_bounds[[cat_col]] <- cat_str
      df <- df_added
    } else { # age based time scale
      # Create the age columns
      # Get the birth date values
      if (cat_df %in% names(df)) { # check that the category doesn't already exist in the original dataframe
        stop("Error: ", cat_df, " already exists, use a different name to avoid overwritting the existing column.")
      }
      birth <- make_date(day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df)))
      birth_list <- time_scale[[cat]]
      if ("year" %in% names(birth_list)) { # nocov
        if ("month" %in% names(birth_list)) { # nocov
          if (any(abs(df[[birth_list$month]] - 6.5) > 5.5)) {
            stop("Error: Atleast one month was not valid")
          }
          if ("day" %in% names(birth_list)) { # nocov

            birth <- make_date(year = df[[birth_list$year]], month = df[[birth_list$month]], day = df[[birth_list$day]]) # nocov
          } else { # nocov
            birth <- make_date(year = df[[birth_list$year]], month = df[[birth_list$month]], day = rep(day_default, nrow(df))) # nocov
          }
        } else { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(year = df[[birth_list$year]], day = df[[birth_list$day]], month = rep(month_default, nrow(df))) # nocov
          } else { # nocov
            birth <- make_date(year = df[[birth_list$year]], day = rep(day_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
          }
        }
      } else {
        if ("month" %in% names(birth_list)) { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(month = df[[birth_list$month]], day = df[[birth_list$day]], year = rep(year_default, nrow(df))) # nocov
          } else { # nocov
            birth <- make_date(month = df[[birth_list$month]], year = rep(year_default, nrow(df)), day = rep(day_default, nrow(df))) # nocov
          }
        } else { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(day = df[[birth_list$day]], year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
          } else {
            stop("Error: birth date missing day, month, and year")
          }
        }
      }
      if (any(is.na(birth))) {
        stop("Error: Birth date failed. Check dates for any impossible values.")
      }
      # Get the entry/exit years
      interval <- "interval"
      pyr_entry <- pyr$entry
      pyr_exit <- pyr$exit
      entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
      exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
      if (pyr_entry$trunc) {
        entry <- birth
      }
      if (pyr_exit$trunc) {
        exit <- entry + years(1)
      }
      if (any(is.na(entry))) {
        stop("Error: Entry date failed. Check dates for any impossible values.")
      }
      if (any(is.na(exit))) {
        stop("Error: Exit date failed. Check dates for any impossible values.")
      }
      #
      df[["entry"]] <- entry
      df[["exit"]] <- exit
      df[["birth"]] <- birth
      df[["interval_dur"]] <- as.numeric(as.duration(interval(entry, exit)), pyr_unit)
      df[["entry_age"]] <- as.numeric(as.duration(interval(birth, entry)), pyr_unit)
      df[["exit_age"]] <- as.numeric(as.duration(interval(birth, exit)), pyr_unit)
      df <- df |> filter(interval_dur >= 0)
      if (nrow(df) < 1) {
        stop("Error: The table passed had no data at risk after filtering for negative durations.")
      }
      df <- df |> filter(entry_age >= 0)
      if (nrow(df) < 1) {
        stop("Error: The table passed had no data at risk after filtering for negative `entry age`.")
      }
      birth <- df$birth
      entry <- df$entry
      exit <- df$exit
      # Create the categories and assign
      Ls <- c()
      Us <- c()
      df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
      df <- df %>% mutate("{cat_df}" := "Unassigned") # initialize the new "age" column
      df_added <- tibble()
      if ("lower" %in% names(time_scale[[cat]])) { # lower and upper boundary intervals
        temp0 <- time_scale[[cat]]$lower
        temp1 <- time_scale[[cat]]$upper
        if ("name" %in% names(time_scale[[cat]])) { # check for names for each level
          temp2 <- time_scale[[cat]]$name
        } else {
          temp2 <- seq_along(temp0)
        }
        num_categ <- length(temp0)
        categ_cols <- c(categ_cols, cat_col)
        for (i in 1:num_categ) { # for each level
          L <- as.numeric(temp0[i])
          Ls <- c(Ls, L)
          if (grepl("]", temp1[i], fixed = TRUE)) { # check for including the upper limit
            stop("Error: Time category included upper limit, currently not supported to avoid double-counting events.")
          } else {
            U <- as.numeric(temp1[i]) # get upper limit
            Us <- c(Us, U)
          }
        }
      } else { # boundary as string
        if (!("categories" %in% names(time_scale[[cat]]))) {
          stop("Error: Age categories should be named `categories` in the list.")
        }
        temp <- time_scale[[cat]]$categories
        temp <- gsub("/", " / ", temp, fixed = TRUE)
        temp <- gsub("]", " ] ", temp, fixed = TRUE)
        temp <- strsplit(temp, "\\s+")[[1]] # seperate values and delimiters
        match_index <- which(temp %in% c("/", "]"))
        if (length(match_index) < 1) {
          stop(paste("Error: Category ", cat, " did not have categories.", sep = ""))
        }
        categ_cols <- c(categ_cols, cat_col)
        match_i <- 1
        for (time_i in match_index) {
          # We want to figure out if there is a label
          L <- as.numeric(temp[time_i - 1])
          Ls <- c(Ls, L)
          if (time_i == match_index[length(match_index)]) {
            # We are in the last entry
            if (length(temp) - time_i == 2) {
              # There is a label
              stop("Error: Age categories currently do not support names.")
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              Us <- c(Us, U)
            }
          } else {
            # Not at last entry
            if (match_index[match_i + 1] - time_i == 3) {
              # There is a label
              stop("Error: Age categories currently do not support names.")
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              Us <- c(Us, U)
            }
          }
          match_i <- match_i + 1
        }
      }
      num_categ <- length(Ls)
      for (time_i in 1:num_categ) {
        L <- Ls[time_i]
        U <- Us[time_i]
        if (U < L) {
          stop(paste0("Error: Age category went from", L, " to ", U, " and was not increasing."))
        }
        ## We convert to a calendar based comparison
        # We need to find the entry and exit dates, convert L and U to integer/decimal and add to birth
        if (pyr_unit == "years") {
          # We can start by just adding the years
          L_year <- floor(L)
          L_decimal <- L - L_year
          U_year <- floor(U)
          U_decimal <- U - U_year
          istart <- birth %m+% lubridate::years(L_year)
          iend <- birth %m+% lubridate::years(U_year)
          # Now we want to add in the correct number of days
          istart <- istart %m+% lubridate::days(round(L_decimal * 365.242199))
          iend <- iend %m+% lubridate::days(round(U_decimal * 365.242199))
        } else if (pyr_unit == "months") {
          # We start by adding months
          L_month <- floor(L)
          L_decimal <- L - L_month
          U_month <- floor(U)
          U_decimal <- U - U_month
          istart <- birth %m+% months(L_month)
          iend <- birth %m+% months(U_month)
          # Now we want to add in the correct number of days
          istart <- istart %m+% lubridate::days(round(L_decimal * 365.242199 / 12))
          iend <- iend %m+% lubridate::days(round(U_decimal * 365.242199 / 12))
        } else if (pyr_unit == "days") {
          # The smallest we can use are days, so round them accordingly
          L <- round(L)
          U <- round(U)
          #
          istart <- birth %m+% lubridate::days(L)
          iend <- birth %m+% lubridate::days(U)
        }
        # We don't want the upper limit to be inclusive, so we call roll it back 1 day
        # Start by checking if the category is one day
        bin_dur <- as.duration(interval(istart, iend))
        # The lowest unit of time is days, so we want to adjust any interval longer than 0 days down one, excluding the last time interval which includes the upper limit
        if (time_i < num_categ) {
          iend <- case_when(bin_dur > days(1) ~ iend - days(1), .default = iend)
        }
        if (pyr_exit$trunc) {
          # we do not have an exit time
          exit <- iend + years(1)
        }
        cat_str <- paste(cat_str, paste("[", L, " to ", U, ")", sep = ""), sep = " ") # prepare the interval info
        c_categ <- list()
        categ_interval <- interval(istart, iend) # define as date interval
        # both entry and exit
        risk_interval <- interval(entry, exit)
        # When we are not in final category, we correct the interval end time
        day_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::day(istart), .default = 0)
        month_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::month(istart), .default = 0)
        year_enter_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::year(istart), .default = 0)
        #
        day_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::day(iend %m+% lubridate::days(1)), .default = 0)
        month_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::month(iend %m+% lubridate::days(1)), .default = 0)
        year_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ lubridate::year(iend %m+% lubridate::days(1)), .default = 0)
        if (time_i == num_categ) {
          end_dif <- as.duration(interval(iend, exit)) # check if the risk and category intervals end at the same time
          #
          day_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::day(iend), .default = 0)
          month_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::month(iend), .default = 0)
          year_exit_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & end_dif < days(1) & df[[cat_col]] == "Unassigned" ~ lubridate::year(iend), .default = 0)
        }
        a_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # category fully contained, set to i
        #
        day_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(entry), .default = day_enter_categ)
        month_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(entry), .default = month_enter_categ)
        year_enter_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(entry), .default = year_enter_categ)
        if (time_i < num_categ) {
          day_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend %m+% lubridate::days(1)), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend %m+% lubridate::days(1)), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend %m+% lubridate::days(1)), .default = year_exit_categ)
        } else {
          end_dif <- as.duration(interval(iend, exit)) # check if the risk and category intervals end at the same time
          day_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend %m+% lubridate::days(1)), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend %m+% lubridate::days(1)), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & end_dif >= days(1) & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend %m+% lubridate::days(1)), .default = year_exit_categ)
          #
          day_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(iend), .default = day_exit_categ)
          month_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(iend), .default = month_exit_categ)
          year_exit_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(iend), .default = year_exit_categ)
        }
        a_categ <- case_when(iend %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval ends during row interval, set to i
        #
        day_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(istart), .default = day_enter_categ)
        day_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(exit), .default = day_exit_categ)
        month_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(istart), .default = month_enter_categ)
        month_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(exit), .default = month_exit_categ)
        year_enter_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(istart), .default = year_enter_categ)
        year_exit_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(exit), .default = year_exit_categ)
        a_categ <- case_when(istart %within% risk_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval starts during row interval, set to i
        #
        day_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & day_enter_categ == 0 ~ lubridate::day(entry), .default = day_enter_categ)
        day_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & day_exit_categ == 0 ~ lubridate::day(exit), .default = day_exit_categ)
        month_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & month_enter_categ == 0 ~ lubridate::month(entry), .default = month_enter_categ)
        month_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & month_exit_categ == 0 ~ lubridate::month(exit), .default = month_exit_categ)
        year_enter_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & year_enter_categ == 0 ~ lubridate::year(entry), .default = year_enter_categ)
        year_exit_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" & year_exit_categ == 0 ~ lubridate::year(exit), .default = year_exit_categ)
        a_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # row interval fully contained in category interval, set to i
        #
        for (evt in event_cols) {
          c_categ[[evt]] <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # row ends during category interval, set to event value
        }
        b_categ <- case_when(entry %within% categ_interval & exit %within% categ_interval & df[[cat_col]] == "Unassigned" ~ as.numeric(as.duration(risk_interval), pyr_unit), .default = -1) # row interval fully in category interval, track full row interval
        b_categ <- case_when(exit %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(istart, exit)), pyr_unit), .default = b_categ) # rows which end during the category interval, track category interval start to row end
        if (time_i < num_categ) {
          b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(entry, iend)) + ddays(1), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
          b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(categ_interval) + ddays(1), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
        } else {
          b_categ <- case_when(entry %within% categ_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(interval(entry, iend)), pyr_unit), .default = b_categ) # rows which enter during the category interval, track entry to category interval end
          b_categ <- case_when(istart %within% risk_interval & iend %within% risk_interval & df[[cat_col]] == "Unassigned" & b_categ == -1 ~ as.numeric(as.duration(categ_interval), pyr_unit), .default = b_categ) # category interval fully in row interval, track full category interval
        }
        b_categ <- case_when(b_categ == -1 ~ 0.0, .default = b_categ) # set every unused interval to 0
        #
        index_kept <- seq_len(nrow(df))
        index_kept <- index_kept[a_categ == time_i] # indexes which contain the category interval to some level
        b_categ <- b_categ[a_categ == time_i]
        index_kept <- index_kept[b_categ >= 0]
        #
        day_enter_categ <- day_enter_categ[a_categ == time_i]
        day_exit_categ <- day_exit_categ[a_categ == time_i]
        month_enter_categ <- month_enter_categ[a_categ == time_i]
        month_exit_categ <- month_exit_categ[a_categ == time_i]
        year_enter_categ <- year_enter_categ[a_categ == time_i]
        year_exit_categ <- year_exit_categ[a_categ == time_i]
        #
        day_enter_categ <- day_enter_categ[b_categ >= 0]
        day_exit_categ <- day_exit_categ[b_categ >= 0]
        month_enter_categ <- month_enter_categ[b_categ >= 0]
        month_exit_categ <- month_exit_categ[b_categ >= 0]
        year_enter_categ <- year_enter_categ[b_categ >= 0]
        year_exit_categ <- year_exit_categ[b_categ >= 0]
        #
        row_kept <- slice(df, index_kept) # dataframe at kept rows
        row_kept <- row_kept %>% mutate("{cat_col}" := as.character(time_i)) # time category value
        #
        row_kept <- row_kept %>% mutate("{y_entry}" := year_enter_categ)
        row_kept <- row_kept %>% mutate("{m_entry}" := month_enter_categ)
        row_kept <- row_kept %>% mutate("{d_entry}" := day_enter_categ)
        row_kept <- row_kept %>% mutate("{y_exit}" := year_exit_categ)
        row_kept <- row_kept %>% mutate("{m_exit}" := month_exit_categ)
        row_kept <- row_kept %>% mutate("{d_exit}" := day_exit_categ)
        #

        # {cat_df}
        #
        for (evt in event_cols) {
          d_categ <- c_categ[[evt]]
          d_categ <- d_categ[a_categ == time_i]
          d_categ <- d_categ[b_categ >= 0]
          row_kept <- row_kept %>% mutate("{evt}" := d_categ) # event values
        }
        #
        df_added <- bind_rows(df_added, row_kept) # new updates dataset
      }
      df <- df_added
      # now we calculate the `age` column values
      entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
      exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
      birth <- make_date(day = rep(day_default, nrow(df)), year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df)))
      birth_list <- time_scale[[cat]]
      if ("year" %in% names(birth_list)) { # nocov
        if ("month" %in% names(birth_list)) { # nocov
          if (any(abs(df[[birth_list$month]] - 6.5) > 5.5)) {
            stop("Error: Atleast one month was not valid")
          }
          if ("day" %in% names(birth_list)) { # nocov

            birth <- make_date(year = df[[birth_list$year]], month = df[[birth_list$month]], day = df[[birth_list$day]]) # nocov
          } else { # nocov
            birth <- make_date(year = df[[birth_list$year]], month = df[[birth_list$month]], day = rep(day_default, nrow(df))) # nocov
          }
        } else { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(year = df[[birth_list$year]], day = df[[birth_list$day]], month = rep(month_default, nrow(df))) # nocov
          } else { # nocov
            birth <- make_date(year = df[[birth_list$year]], day = rep(day_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
          }
        }
      } else {
        if ("month" %in% names(birth_list)) { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(month = df[[birth_list$month]], day = df[[birth_list$day]], year = rep(year_default, nrow(df))) # nocov
          } else { # nocov
            birth <- make_date(month = df[[birth_list$month]], year = rep(year_default, nrow(df)), day = rep(day_default, nrow(df))) # nocov
          }
        } else { # nocov
          if ("day" %in% names(birth_list)) { # nocov
            birth <- make_date(day = df[[birth_list$day]], year = rep(year_default, nrow(df)), month = rep(month_default, nrow(df))) # nocov
          } else {
            stop("Error: birth date missing day, month, and year")
          }
        }
      }
      if (pyr_entry$trunc) {
        entry <- birth
        df[[cat_df]] <- as.numeric(as.duration(interval(birth, exit)), pyr_unit)
      } else if (pyr_exit$trunc) {
        df[[cat_df]] <- as.numeric(as.duration(interval(entry, birth)), pyr_unit)
      } else {
        df[[cat_df]] <- (as.numeric(as.duration(interval(birth, exit)), pyr_unit) + as.numeric(as.duration(interval(birth, entry)), pyr_unit)) / 2
      }
      #
      categ_bounds[[cat_col]] <- cat_str
    }
    # after one run through, there should be no truncation
    pyr_entry$trunc <- FALSE
    pyr_exit$trunc <- FALSE
    if (nrow(df) < 1) {
      stop("Error: The table passed had no data at risk after applying time categorization.")
    }
  }
  # Now we find the first and last entry for each ID
  # They can be compared by entry age within each study id, relative to some arbitrary point
  # we don't care if the value is negative, just which dates are earlier or later
  # start by getting the entry and exit dates
  pyr_entry <- pyr$entry
  pyr_exit <- pyr$exit
  entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
  exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
  # now we take a duration from an arbitrary realistic date
  arbitrary_date <- entry[1]
  if ("_arbitrary_age" %in% names(df)) {
    warning("Warning: `_arbitrary_age` is a restricted name and will be overwritten.")
  }
  df[["_arbitrary_age"]] <- as.numeric(as.duration(interval(arbitrary_date, entry)), "years")
  #
  setDT(df)
  setkeyv(df, c(studyid, "_arbitrary_age"))
  # The intervals earlier on should correspond to earlier observations, later on corresponding to later observations
  # We want it to be ordered by studyid first, so that we can fill out a vector sequentially
  id_vals <- unique(df[[studyid]])
  f_categ <- c()
  l_categ <- c()
  for (id_val in id_vals) {
    df_id_len <- nrow(df[df[[studyid]] == id_val, ])
    f_categ <- c(f_categ, 1.0, rep(0.0, df_id_len - 1))
    l_categ <- c(l_categ, rep(0.0, df_id_len - 1), 1.0)
  }
  df[["F_AT_RISK"]] <- f_categ
  df[["L_AT_RISK"]] <- l_categ
  df <- tibble(df)
  # Now we get the person-years
  pyr_entry <- pyr$entry
  pyr_exit <- pyr$exit
  entry <- make_date(year = df[[pyr_entry$year]], month = df[[pyr_entry$month]], day = df[[pyr_entry$day]]) # nocov
  exit <- make_date(year = df[[pyr_exit$year]], month = df[[pyr_exit$month]], day = df[[pyr_exit$day]]) # nocov
  risk_interval <- interval(entry, exit)
  pyr_unit <- "years" # default person-years to years
  if ("unit" %in% names(pyr)) {
    pyr_unit <- pyr$unit
  }
  df <- df %>% mutate("PYR" := as.numeric(as.duration(risk_interval), pyr_unit))
  #
  list(df = df, categ_cols = categ_cols, categ_bounds = categ_bounds)
}

#' uses a table and user time scale categories to split a table into risk intervals
#'
#' \code{UserScale_Process} generates table split by risk categories
#'
#' @inheritParams R_template
#' @param df dataframe with every category/event column needed
#' @param pyr list with entry and exit lists, containing day/month/year columns in the table
#' @param time_scale list with the time scale information, either a calendar category or an age category
#' @param event_cols list of event columns
#' @param studyid id used to determine distinct subjects used for first and last at risk intervals
#'
#' @return returns calendar based table generation
#' @family Table Generation Functions
#' @noRd
UserScale_Process <- function(df, table_names, pyr = list(), time_scale = list(), event_cols = c(), studyid = "studyID", verbose = FALSE) {
  # setting commands and known columns to dummy values to avoid `no visible binding` warnings
  `%within%` <- lubridate::`%within%`
  `%>%` <- dplyr::`%>%`
  `%m+%` <- lubridate::`%m+%`
  interval_dur <- ""
  entry_age <- ""
  #
  pyr_entry <- pyr$entry
  pyr_exit <- pyr$exit
  #
  categ_cols <- c()
  categ_bounds <- list()
  for (cat in names(time_scale)) { # for each time category
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
    names(time_scale[[cat]]) <- tolower(names(time_scale[[cat]]))
    #
    if (any(names(time_scale[[cat]]) %in% c("day", "month", "year"))) {
      stop("Error: Time scale used calender naming, but person-year data used a user time-scale. User defined categories are defined similar to non-time categories (upper/lower limits).")
    }
    # Create the duration column
    if (cat_df %in% names(df)) { # check that the category doesn't already exist in the original dataframe
      stop("Error: ", cat_df, " already exists, use a different name to avoid overwritting the existing column.")
    }
    pyr_entry <- pyr$entry
    pyr_exit <- pyr$exit
    #
    entry_duration <- pyr_entry$duration
    exit_duration <- pyr_exit$duration
    # Get the entry and exit durations
    entry <- df[[entry_duration]]
    exit <- df[[exit_duration]]
    if (pyr_exit$trunc) {
      exit <- entry + 1
    }
    if (any(is.na(entry))) {
      stop("Error: Entry duration had NA.")
    }
    if (any(is.na(exit))) {
      stop("Error: Exit duration had NA.")
    }
    #
    df[["entry"]] <- entry
    df[["exit"]] <- exit
    df[["interval_dur"]] <- exit - entry
    df <- df |> filter(interval_dur >= 0)
    if (nrow(df) < 1) {
      stop("Error: The table passed had no data at risk after filtering for negative durations.")
    }
    entry <- df$entry
    exit <- df$exit
    # Create the categories and assign
    Ls <- c()
    Us <- c()
    df <- df %>% mutate("{cat_col}" := "Unassigned") # initialize the tibble
    df <- df %>% mutate("{cat_df}" := "Unassigned") # initialize the new "age" column
    df_added <- tibble()
    if (inherits(time_scale[[cat]], "list")) {
      if ("lower" %in% names(time_scale[[cat]])) { # lower and upper boundary intervals
        temp0 <- time_scale[[cat]]$lower
        temp1 <- time_scale[[cat]]$upper
        if ("name" %in% names(time_scale[[cat]])) { # check for names for each level
          temp2 <- time_scale[[cat]]$name
        } else {
          temp2 <- seq_along(temp0)
        }
        num_categ <- length(temp0)
        categ_cols <- c(categ_cols, cat_col)
        for (i in 1:num_categ) { # for each level
          L <- as.numeric(temp0[i])
          Ls <- c(Ls, L)
          if (grepl("]", temp1[i], fixed = TRUE)) { # check for including the upper limit
            stop("Error: Time category included upper limit, currently not supported to avoid double-counting events.")
          } else {
            U <- as.numeric(temp1[i]) # get upper limit
            Us <- c(Us, U)
          }
        }
      } else { # boundary as string
        if (!("categories" %in% names(time_scale[[cat]]))) {
          stop("Error: Age categories should be named `categories` in the list.")
        }
        temp <- time_scale[[cat]]$categories
        temp <- gsub("/", " / ", temp, fixed = TRUE)
        temp <- gsub("]", " ] ", temp, fixed = TRUE)
        temp <- strsplit(temp, "\\s+")[[1]] # seperate values and delimiters
        match_index <- which(temp %in% c("/", "]"))
        if (length(match_index) < 1) {
          stop(paste("Error: Category ", cat, " did not have categories.", sep = ""))
        }
        categ_cols <- c(categ_cols, cat_col)
        match_i <- 1
        for (time_i in match_index) {
          # We want to figure out if there is a label
          L <- as.numeric(temp[time_i - 1])
          Ls <- c(Ls, L)
          if (time_i == match_index[length(match_index)]) {
            # We are in the last entry
            if (length(temp) - time_i == 2) {
              # There is a label
              stop("Error: Age categories currently do not support names.")
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              Us <- c(Us, U)
            }
          } else {
            # Not at last entry
            if (match_index[match_i + 1] - time_i == 3) {
              # There is a label
              stop("Error: Age categories currently do not support names.")
            } else {
              # No label
              U <- as.numeric(temp[time_i + 1])
              Us <- c(Us, U)
            }
          }
          match_i <- match_i + 1
        }
      }
    } else {
      temp <- time_scale[[cat]]
      temp <- gsub("/", " / ", temp, fixed = TRUE)
      temp <- gsub("]", " ] ", temp, fixed = TRUE)
      temp <- strsplit(temp, "\\s+")[[1]] # seperate values and delimiters
      match_index <- which(temp %in% c("/", "]"))
      if (length(match_index) < 1) {
        stop(paste("Error: Category ", cat, " did not have categories.", sep = ""))
      }
      categ_cols <- c(categ_cols, cat_col)
      match_i <- 1
      for (time_i in match_index) {
        # We want to figure out if there is a label
        L <- as.numeric(temp[time_i - 1])
        Ls <- c(Ls, L)
        if (time_i == match_index[length(match_index)]) {
          # We are in the last entry
          if (length(temp) - time_i == 2) {
            # There is a label
            stop("Error: Age categories currently do not support names.")
          } else {
            # No label
            U <- as.numeric(temp[time_i + 1])
            Us <- c(Us, U)
          }
        } else {
          # Not at last entry
          if (match_index[match_i + 1] - time_i == 3) {
            # There is a label
            stop("Error: Age categories currently do not support names.")
          } else {
            # No label
            U <- as.numeric(temp[time_i + 1])
            Us <- c(Us, U)
          }
        }
        match_i <- match_i + 1
      }
    }
    num_categ <- length(Ls)
    for (time_i in 1:num_categ) {
      L <- Ls[time_i]
      U <- Us[time_i]
      if (U < L) {
        stop(paste0("Error: Age category went from", L, " to ", U, " and was not increasing."))
      }
      if (pyr_exit$trunc) {
        # we do not have an exit time
        exit <- rep(U + 1, length(entry))
      }
      if (pyr_entry$trunc) {
        # we do not have an exit time
        entry <- rep(L - 1, length(exit))
      }
      cat_str <- paste(cat_str, paste("[", L, " to ", U, ")", sep = ""), sep = " ") # prepare the interval info
      c_categ <- list()
      #
      entry_categ <- case_when(entry <= L & U <= exit & df[[cat_col]] == "Unassigned" ~ L, .default = NA)
      exit_categ <- case_when(entry <= L & U <= exit & df[[cat_col]] == "Unassigned" ~ U, .default = NA)
      a_categ <- case_when(entry <= L & U <= exit & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = "0") # category fully contained, set to i
      #
      entry_categ <- case_when(entry <= U & U <= exit & df[[cat_col]] == "Unassigned" & is.na(entry_categ) ~ entry, .default = entry_categ)
      exit_categ <- case_when(entry <= U & U <= exit & df[[cat_col]] == "Unassigned" & is.na(exit_categ) ~ U, .default = exit_categ)
      a_categ <- case_when(entry <= U & U <= exit & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval ends during row interval, set to i
      #
      entry_categ <- case_when(entry <= L & L <= exit & df[[cat_col]] == "Unassigned" & is.na(entry_categ) ~ L, .default = entry_categ)
      exit_categ <- case_when(entry <= L & L <= exit & df[[cat_col]] == "Unassigned" & is.na(exit_categ) ~ exit, .default = exit_categ)
      a_categ <- case_when(entry <= L & L <= exit & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # interval starts during row interval, set to i
      #
      entry_categ <- case_when(L <= entry & exit < U & df[[cat_col]] == "Unassigned" & is.na(entry_categ) ~ entry, .default = entry_categ)
      exit_categ <- case_when(L <= entry & exit < U & df[[cat_col]] == "Unassigned" & is.na(exit_categ) ~ exit, .default = exit_categ)
      a_categ <- case_when(L <= entry & exit < U & df[[cat_col]] == "Unassigned" ~ as.character(time_i), .default = a_categ) # row interval fully contained in category interval, set to i
      #
      for (evt in event_cols) {
        c_categ[[evt]] <- case_when(L <= exit & exit < U & df[[cat_col]] == "Unassigned" ~ df[[evt]], .default = 0) # row ends during category interval, set to event value
      }
      b_categ <- case_when(a_categ == as.character(time_i) ~ exit_categ - entry_categ, .default = 0)
      #
      index_kept <- seq_len(nrow(df))
      index_kept <- index_kept[a_categ == time_i] # indexes which contain the category interval to some level
      b_categ <- b_categ[a_categ == time_i]
      index_kept <- index_kept[b_categ >= 0]
      #
      entry_categ <- entry_categ[a_categ == time_i]
      entry_categ <- entry_categ[b_categ >= 0]
      exit_categ <- exit_categ[a_categ == time_i]
      exit_categ <- exit_categ[b_categ >= 0]
      #
      row_kept <- slice(df, index_kept) # dataframe at kept rows
      #
      row_kept <- row_kept %>% mutate("{cat_col}" := as.character(time_i)) # time category value
      row_kept <- row_kept %>% mutate("{entry_duration}" := entry_categ) # entry duration
      row_kept <- row_kept %>% mutate("{exit_duration}" := exit_categ) # entry duration
      #
      for (evt in event_cols) {
        d_categ <- c_categ[[evt]]
        d_categ <- d_categ[a_categ == time_i]
        d_categ <- d_categ[b_categ >= 0]
        row_kept <- row_kept %>% mutate("{evt}" := d_categ) # event values
      }
      #
      df_added <- bind_rows(df_added, row_kept) # new updates dataset
    }
    df <- df_added
    # now we calculate the `age` column values
    df[[cat_df]] <- (df[[entry_duration]] + df[[exit_duration]]) / 2
    #
    categ_bounds[[cat_col]] <- cat_str
    # after one run through, there should be no truncation
    pyr_entry$trunc <- FALSE
    pyr_exit$trunc <- FALSE
    if (nrow(df) < 1) {
      stop("Error: The table passed had no data at risk after applying time categorization.")
    }
  }
  # Now we find the first and last entry for each ID
  pyr_entry <- pyr$entry
  pyr_exit <- pyr$exit
  #
  setDT(df)
  setkeyv(df, c(studyid, pyr_entry$duration))
  # The intervals earlier on should correspond to earlier observations, later on corresponding to later observations
  # We want it to be ordered by studyid first, so that we can fill out a vector sequentially
  id_vals <- unique(df[[studyid]])
  f_categ <- c()
  l_categ <- c()
  for (id_val in id_vals) {
    df_id_len <- nrow(df[df[[studyid]] == id_val, ])
    f_categ <- c(f_categ, 1.0, rep(0.0, df_id_len - 1))
    l_categ <- c(l_categ, rep(0.0, df_id_len - 1), 1.0)
  }
  df[["F_AT_RISK"]] <- f_categ
  df[["L_AT_RISK"]] <- l_categ
  df <- tibble(df)
  # Now we get the new duration
  df[["PYR"]] <- df[[pyr_exit$duration]] - df[[pyr_entry$duration]]
  #
  list(df = df, categ_cols = categ_cols, categ_bounds = categ_bounds)
}

#' uses a table and summaries to construct a grouped table
#'
#' \code{generate_summaries} generates table grouped with summaries
#'
#' @inheritParams R_template
#'
#' @return returns grouped table
#' @family Table Generation Functions
#' @noRd
generate_summaries <- function(df, summaries, event_cols, event_names, categ_cols, fcount = FALSE, lcount = FALSE, time_table = FALSE) {
  `%>%` <- dplyr::`%>%`
  # Makes a list of named changed events, to allow them to be used with weighting
  evt_name_change <- list()
  for (i in seq_along(event_cols)) {
    evt_name_change[[event_names[i]]] <- event_cols[[i]]
  }
  if (time_table) {
    df_group <- df %>%
      group_by(across(all_of(categ_cols))) %>%
      summarize(AT_RISK = n(), PYR = sum(.data[["PYR"]]), .groups = "drop") # group by categories and define the durations and counts
    if (fcount) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize(F_AT_RISK = sum(.data[["F_AT_RISK"]]), .groups = "drop")
      df_group[["F_AT_RISK"]] <- df_temp[["F_AT_RISK"]]
    }
    if (lcount) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize(L_AT_RISK = sum(.data[["L_AT_RISK"]]), .groups = "drop")
      df_group[["L_AT_RISK"]] <- df_temp[["L_AT_RISK"]]
    }
  } else {
    df_group <- df %>%
      group_by(across(all_of(categ_cols))) %>%
      summarize(COUNT = n(), .groups = "drop") # group by columns and summarize by counts
  }
  for (evt_i in seq_along(names(summaries))) { # for each event summary
    evt <- names(summaries)[evt_i]
    if (!(evt %in% names(df))) {
      stop(paste0("Error: The column `", evt, "` was not present in the data for summary"))
    }
    col_name <- evt
    weight <- "NULL"
    temp <- strsplit(summaries[[evt_i]], "\\s+")[[1]]
    if (length(temp) == 0) {
      # nothing
      stop("Error: Missing summary method")
    } else if (length(temp) == 1) {
      # only method
      method <- tolower(temp)
    } else if (length(temp) == 3) {
      # method and some change
      method <- tolower(temp[1])
      temp[2] <- tolower(temp[2])
      if (!(temp[2] %in% c("by", "as", "weight"))) {
        stop("Error: Only AS, WEIGHT, and BY can be used to change summaries.")
      }
      if (temp[2] %in% c("by", "weight")) {
        weight <- temp[3]
      } else if (temp[2] == "as") {
        col_name <- temp[3]
      }
    } else if (length(temp) == 5) {
      # method and both changes
      method <- tolower(temp[1])
      temp[2] <- tolower(temp[2])
      temp[4] <- tolower(temp[4])
      if (temp[2] == "weight") {
        temp[2] <- "by"
      }
      if (temp[4] == "weight") {
        temp[4] <- "by"
      }
      if (temp[2] == temp[4]) {
        stop(paste0("Error: ", temp[2], " used twice."))
      }
      if (!any(c(temp[2], temp[4]) %in% c("by", "as"))) {
        stop("Error: Only AS, WEIGHT, and BY can be used to change summaries.")
      }
      if (temp[2] == "by") {
        weight <- temp[3]
        col_name <- temp[5]
      } else if (temp[2] == "as") {
        col_name <- temp[3]
        weight <- temp[5]
      }
    } else {
      stop("Error: Summaries should be listed as `method AS new_name BY weight` with optional new name and weight")
    }
    if (time_table) {
      acceptable_methods <- c("count", "sum", "mean", "rsum", "rmean", "weighted_mean", "xsum", "xmean", "weighted_sum", "weighted_xmean", "weighted_xsum")
    } else {
      acceptable_methods <- c("count", "sum", "mean", "rsum", "rmean", "weighted_mean", "weighted_sum")
    }
    method <- lapply(method, function(x) tryCatch(match.arg(x, choices = acceptable_methods), error = function(error_message) x))[[1]] # match against expected values
    if (!(method %in% acceptable_methods)) {
      stop(paste0("Error: Provided method of `", method, "` was not valid."))
    }
    if (weight != "NULL") {
      if (!(weight %in% names(df))) {
        if (weight %in% names(evt_name_change)) {
          weight <- evt_name_change[[weight]]
        } else {
          stop(paste0("Error: Weighting column, `", weight, "` was not in the data."))
        }
      }
      if (method %in% c("rmean")) {
        method <- "weighted_mean" # update to the weighted mean if a weighting column is given
      } else if (method %in% c("mean", "sum", "xmean", "xsum")) {
        method <- paste0("weighted_", method)
      } else if (!(method %in% c("weighted_sum", "weighted_mean", "mean", "rmean"))) {
        warning(paste0("Warning: A weighting column was given, but the `", method, "` method does not support weighting."))
      }
    } else {
      if ((!time_table) && grepl("weighted_", method)) {
        stop("Error: No default weighting for event-count tables.")
      }
      weight <- "PYR"
    }
    if (method %in% c("count")) { # measuring instances, assumes the column is binary and returns the sum
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method %in% c("sum", "rsum")) { # sum of event across each category combination
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]]), .groups = "drop")
    } else if (method %in% c("weighted_sum")) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]] * .data[[weight]]), .groups = "drop")
    } else if (method %in% c("xsum")) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]] * .data[["L_AT_RISK"]]), .groups = "drop")
    } else if (method %in% c("xmean")) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := weighted.mean(.data[[evt]], .data[["L_AT_RISK"]]), .groups = "drop")
      df_temp[[col_name]] <- case_when(is.nan(df_temp[[col_name]]) ~ 0, .default = df_temp[[col_name]])
    } else if (method %in% c("weighted_xsum")) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := sum(.data[[evt]] * .data[["L_AT_RISK"]] * .data[[weight]]), .groups = "drop")
    } else if (method %in% c("weighted_xmean")) {
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := weighted.mean(.data[[evt]], .data[["L_AT_RISK"]] * .data[[weight]]), .groups = "drop")
      df_temp[[col_name]] <- case_when(is.nan(df_temp[[col_name]]) ~ 0, .default = df_temp[[col_name]])
    } else if (method %in% c("mean", "rmean")) { # mean value of event across each category combination
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := mean(.data[[evt]]), .groups = "drop")
    } else if (method == "weighted_mean") { # mean value weighted by person-years
      df_temp <- df %>%
        group_by(across(all_of(categ_cols))) %>%
        summarize("{col_name}" := weighted.mean(.data[[evt]], .data[[weight]]), .groups = "drop")
      df_temp[[col_name]] <- case_when(is.nan(df_temp[[col_name]]) ~ 0, .default = df_temp[[col_name]])
    } else {
      stop(paste0("Error: method, `", method, "` was missed"))
    }
    df_group[[col_name]] <- df_temp[[col_name]]
  }
  if (time_table) {
    df_group <- df_group[df_group$PYR > 0, ]
  }
  df_group
}

#' splits a table by non-time categories
#'
#' \code{Category_Process} generates table with categories
#'
#' @inheritParams R_template
#'
#' @return returns table
#' @family Table Generation Functions
#' @noRd
Category_Process <- function(df, table_names, categ, categ_cols, categ_bounds) {
  `%>%` <- dplyr::`%>%`
  #
  names(categ) <- lapply(names(categ), function(x) tryCatch(match.arg(x, choices = names(table)), error = function(error_message) x)) # match against columns in the table
  for (cat in names(categ)) { # for each non-time category
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
      categ_bounds[[cat_col]] <- cat_str # update the list of category boundaries
    } else { # boundary as string
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
    if (nrow(df) < 1) {
      stop("Error: The table passed had no data at risk after applying categorization.")
    }
  }
  list(df = df, categ_cols = categ_cols, categ_bounds = categ_bounds)
}
