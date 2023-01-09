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
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
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
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
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

#' Performs basic poisson regression, with multiple guesses
#' \code{RunPoissonRegression_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
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
#' @param guesses_control list of parameters to control how the guessing works
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data
RunPoissonRegression_Guesses <- function(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols, guesses_control){
    if ("verbose" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$verbose <- FALSE
    }
    if (guesses_control$stata==FALSE){
        if ("CONST" %in% names){
            if ("CONST" %in% names(df)){
                ;
            } else {
                df$CONST <- 1
            }
        }
        all_names <- unique(names)
        #
        df <- Replace_Missing(df,all_names,0.0,control$verbose)
        #
        dfc <- match(names,all_names)
        if (sum(df[,event, with = FALSE])==0){
            if (guesses_control$verbose){
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
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$Iterations
        #
        df_res <- data.table()
        e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
#        print(e$beta_0)
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$Deviation]
        for (it in 1:guesses_control$guesses){
            for (i in 1:length(tform)){
                if (grepl("log",tform[i],fixed=FALSE)){
                    if (guesses_control$loglin_method == "uniform"){
                        a_n[i] <- runif(1,min=guesses_control$loglin_min,max=guesses_control$loglin_max)
                    } else {
                        print("bad")
                        stop()
                    }
                } else {
                    if (guesses_control$lin_method == "uniform"){
                        a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max)
                    } else {
                        print("bad")
                        stop()
                    }
                }
            }
            df_res0 <- data.table()
            e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$Deviation]
            if (is.na(e$Deviation)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            }
            
        }
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        a_n_ind <- which.min(df_res[,get(paste(length(e$beta_0)+1))])
#        print(a_n_ind)
#        stop()
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
        #
        control$maxiter <- iteration0
        control <- Def_Control(control)
        e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
        ;
        return (e)
    } else {
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
        #
        df <- Replace_Missing(df,all_names,0.0,control$verbose)
        #
        dfc <- match(names,all_names)
        if (sum(df[,event, with = FALSE])==0){
            if (guesses_control$verbose){
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
        #
        control <- Def_Control(control)
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$Iterations
        keep_str <- control$keep_strata
        control$keep_strata <- TRUE
        #
        df_res <- data.table()
        e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
#        print(e$beta_0)
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$Deviation]
        for (it in 1:guesses_control$guesses){
            for (i in 1:length(tform)){
                if (grepl("log",tform[i],fixed=FALSE)){
                    if (guesses_control$loglin_method == "uniform"){
                        a_n[i] <- runif(1,min=guesses_control$loglin_min,max=guesses_control$loglin_max)
                    } else {
                        print("bad")
                        stop()
                    }
                } else {
                    if (guesses_control$lin_method == "uniform"){
                        a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max)
                    } else {
                        print("bad")
                        stop()
                    }
                }
            }
            df_res0 <- data.table()
            e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$Deviation]
            if (is.na(e$Deviation)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            }
            
        }
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        #
        a_n_ind <- which.min(df_res[,get(paste(length(e$beta_0)+1))])
#        print(a_n_ind)
#        stop()
        control$keep_strata <- keep_str
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
        e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
        return (e)
    }
    return (NULL)
}
