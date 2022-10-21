#' Performs basic Cox Proportional Hazards regression
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return prints the final results and return null
#' @export
#'
RunCoxRegression <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    dfend <- df[get(event)==1, ]
    tu <- unlist(unique(dfend[,time2, with = FALSE]))
    if (length(tu)==0){
        stop()
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    return (e)
}
