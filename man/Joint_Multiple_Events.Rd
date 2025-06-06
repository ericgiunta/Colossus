% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{Joint_Multiple_Events}
\alias{Joint_Multiple_Events}
\title{Automates creating data for a joint competing risks analysis}
\usage{
Joint_Multiple_Events(
  df,
  events,
  name_list,
  term_n_list = list(),
  tform_list = list(),
  keep_constant_list = list(),
  a_n_list = list()
)
}
\arguments{
\item{df}{a data.table containing the columns of interest}

\item{events}{vector of event column names}

\item{name_list}{list of vectors for columns for event specific or shared model elements, required}

\item{term_n_list}{list of vectors for term numbers for event specific or shared model elements, defaults to term 0}

\item{tform_list}{list of vectors for subterm types for event specific or shared model elements, defaults to loglinear}

\item{keep_constant_list}{list of vectors for constant elements for event specific or shared model elements, defaults to free (0)}

\item{a_n_list}{list of vectors for parameter values for event specific or shared model elements, defaults to term 0}
}
\value{
returns the updated dataframe and model inputs
}
\description{
\code{Joint_Multiple_Events} generates input for a regression with multiple non-independent events and models
}
\examples{
library(data.table)
a <- c(0, 0, 0, 1, 1, 1)
b <- c(1, 1, 1, 2, 2, 2)
c <- c(0, 1, 2, 2, 1, 0)
d <- c(1, 1, 0, 0, 1, 1)
e <- c(0, 1, 1, 1, 0, 0)
df <- data.table("t0" = a, "t1" = b, "e0" = c, "e1" = d, "fac" = e)
time1 <- "t0"
time2 <- "t1"
df$pyr <- df$t1 - df$t0
pyr <- "pyr"
events <- c("e0", "e1")
names_e0 <- c("fac")
names_e1 <- c("fac")
names_shared <- c("t0", "t0")
term_n_e0 <- c(0)
term_n_e1 <- c(0)
term_n_shared <- c(0, 0)
tform_e0 <- c("loglin")
tform_e1 <- c("loglin")
tform_shared <- c("quad_slope", "loglin_top")
keep_constant_e0 <- c(0)
keep_constant_e1 <- c(0)
keep_constant_shared <- c(0, 0)
a_n_e0 <- c(-0.1)
a_n_e1 <- c(0.1)
a_n_shared <- c(0.001, -0.02)
name_list <- list("shared" = names_shared, "e0" = names_e0, "e1" = names_e1)
term_n_list <- list("shared" = term_n_shared, "e0" = term_n_e0, "e1" = term_n_e1)
tform_list <- list("shared" = tform_shared, "e0" = tform_e0, "e1" = tform_e1)
keep_constant_list <- list(
  "shared" = keep_constant_shared,
  "e0" = keep_constant_e0, "e1" = keep_constant_e1
)
a_n_list <- list("shared" = a_n_shared, "e0" = a_n_e0, "e1" = a_n_e1)
val <- Joint_Multiple_Events(
  df, events, name_list, term_n_list,
  tform_list, keep_constant_list, a_n_list
)

}
\seealso{
Other Data Cleaning Functions: 
\code{\link{Check_Dupe_Columns}()},
\code{\link{Check_Trunc}()},
\code{\link{Check_Verbose}()},
\code{\link{Convert_Model_Eq}()},
\code{\link{Correct_Formula_Order}()},
\code{\link{Date_Shift}()},
\code{\link{Def_Control}()},
\code{\link{Def_Control_Guess}()},
\code{\link{Def_model_control}()},
\code{\link{Def_modelform_fix}()},
\code{\link{Event_Count_Gen}()},
\code{\link{Event_Time_Gen}()},
\code{\link{Replace_Missing}()},
\code{\link{Time_Since}()},
\code{\link{factorize}()},
\code{\link{factorize_par}()},
\code{\link{gen_time_dep}()},
\code{\link{interact_them}()}
}
\concept{Data Cleaning Functions}
