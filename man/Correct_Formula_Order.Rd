% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{Correct_Formula_Order}
\alias{Correct_Formula_Order}
\title{Corrects the order of terms/formula/etc}
\usage{
Correct_Formula_Order(
  term_n,
  tform,
  keep_constant,
  a_n,
  names,
  der_iden = 0,
  Cons_Mat = matrix(c(0)),
  Cons_Vec = c(0),
  verbose = FALSE,
  model_control = list(para_number = 0)
)
}
\arguments{
\item{term_n}{term numbers}

\item{tform}{list of string function identifiers, used for linear/step}

\item{keep_constant}{binary values to denote which parameters to change}

\item{a_n}{list of initial parameter values, used to determine number of parameters}

\item{names}{columns for elements of the model, used to identify data columns}

\item{der_iden}{number for the subterm to test derivative at, only used for testing runs with a single varying parameter, should be smaller than total number of parameters}

\item{Cons_Mat}{Matrix containing coefficients for system of linear constraints, formatted as matrix}

\item{Cons_Vec}{Vector containing constants for system of linear constraints, formatted as vector}

\item{verbose}{boolean to control if additional information is printed to the console, also accepts 0/1 integer}

\item{model_control}{controls which alternative model options are used, see Def_model_control() for options and vignette("Alt_Run_opt") for further details}
}
\value{
returns the corrected lists
}
\description{
\code{Correct_Formula_Order} checks the order of formulas given and corrects any ordering issues, orders alphabetically, by term number, etc.
}
\examples{
library(data.table)
## basic example code reproduced from the starting-description vignette

term_n <- c(0,1,1,0,0)
tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
keep_constant <- c(0,0,0,1,0)
a_n <- c(1,2,3,4,5)
names <- c("a","a","a","a","a")
val <- Correct_Formula_Order(term_n, tform, keep_constant,
                             a_n, names, Cons_Mat=matrix(c(0)),
                             Cons_Vec=c(0))
term_n <- val$term_n
tform <- val$tform
keep_constant <- val$keep_constant
a_n <- val$a_n
der_iden <- val$der_iden
names <- val$names

}
\seealso{
Other Data Cleaning Functions: 
\code{\link{Check_Dupe_Columns}()},
\code{\link{Check_Trunc}()},
\code{\link{Date_Shift}()},
\code{\link{Def_Control}()},
\code{\link{Def_Control_Guess}()},
\code{\link{Def_model_control}()},
\code{\link{Def_modelform_fix}()},
\code{\link{Joint_Multiple_Events}()},
\code{\link{Replace_Missing}()},
\code{\link{Time_Since}()},
\code{\link{factorize}()},
\code{\link{factorize_par}()},
\code{\link{gen_time_dep}()},
\code{\link{interact_them}()}
}
\concept{Data Cleaning Functions}
