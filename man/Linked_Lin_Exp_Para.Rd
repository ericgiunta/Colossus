% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{Linked_Lin_Exp_Para}
\alias{Linked_Lin_Exp_Para}
\title{Calculates The Additional Parameter For a linear-exponential formula with known maximum}
\usage{
Linked_Lin_Exp_Para(y, a0, a1_goal, verbose = 0)
}
\arguments{
\item{y}{point formula switch}

\item{a0}{linear slope}

\item{a1_goal}{exponential maximum desired}

\item{verbose}{integer valued 0-4 controlling what information is printed to the terminal. Each level includes the lower levels. 0: silent, 1: errors printed, 2: warnings printed, 3: notes printed, 4: debug information printed. Errors are situations that stop the regression, warnings are situations that assume default values that the user might not have intended, notes provide information on regression progress, and debug prints out C++ progress and intermediate results. The default level is 2 and True/False is converted to 3/0.}
}
\value{
returns parameter used by Colossus
}
\description{
\code{Linked_Lin_Exp_Para} Calculates what the additional parameter would be for a desired maximum
}
\examples{
library(data.table)
y <- 7.6
a0 <- 1.2
a1_goal <- 15
full_paras <- Linked_Lin_Exp_Para(y, a0, a1_goal)

}
