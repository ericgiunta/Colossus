#' Performs basic poisson regression
#' \code{RunPoissonRegression} uses user provided data, person-year/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
RunPoissonRegression <- function(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            ;
        } else {
            df$CONST <- 1
        }
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)
    if (sum(df[,event, with = FALSE])==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(pyr,event)
    #
    control <- Def_Control(control)
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(a_n)-length(names)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
    return (e)
}

#' Performs basic poisson regression with strata effect
#' \code{RunPoissonRegression_STRATA} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' @param Strat_Cols column to stratify by
#'
#' @return returns a list of the final results
#' @export
#'
RunPoissonRegression_STRATA <- function(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols){
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            ;
        } else {
            df$CONST <- 1
        }
    }
    #
    val <- factorize(df, Strat_Cols)
    df <- val$df
    names <- c(names, val$cols)
    Term_max <- fir
    Term_n <- c(Term_n, rep(Term_max,length(val$cols)))
    tform <- c(tform, rep("loglin",length(val$cols)))
    keep_constant <- c(keep_constant, rep(0,length(val$cols)))
    a_n <- c(a_n, runif(length(val$cols),min=-0.01,0.01))
    #
    all_names <- unique(names)
    dfc <- match(names,all_names)
    if (sum(df[,event, with = FALSE])==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(pyr,event)
    #
    control <- Def_Control(control)
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(a_n)-length(names)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
    return (e)
}
