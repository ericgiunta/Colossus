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
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data

RunCoxRegression <- function(df, time1="age_start", time2="age_exit", event="cases", names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), keep_constant=rep(0,length(names)), a_n=rep(0.01,length(names)), modelform="M", fir=0, der_iden=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)){
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
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
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    a_n <- a_n0
    ;
    return (e)
}

#' Performs basic Cox Proportional Hazards calculation with no derivative
#' \code{RunCoxRegression_Single} uses user provided data, time/event columns, vectors specifying the model, and options and returns the log-likelihood
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data

RunCoxRegression_Single <- function(df, time1="age_start", time2="age_exit", event="cases", names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), a_n=rep(0.01,length(names)), modelform="M", fir=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1),keep_constant=rep(0,length(names))){
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
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
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition_single(Term_n,tform,a_n,dfc,x_all, fir, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    a_n <- a_n0
    ;
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with a basic model
#' \code{RunCoxRegression_Basic} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data

RunCoxRegression_Basic <- function(df, time1, time2, event, names, keep_constant, a_n, der_iden, control){
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition_basic(a_n,dfc,x_all, der_iden, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant)
    a_n <- a_n0
    ;
    return (e)
}

#' Performs basic Cox Proportional Hazards regression, Allows for multiple starting guesses
#' \code{RunCoxRegression_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions. Has additional options to starting with several initial guesses
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
#' @param guesses_control list of parameters to control how the guessing works
#' @param Strat_Col column to stratify by if needed
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data
RunCoxRegression_Guesses <- function(df, time1="age_start", time2="age_exit", event="cases", names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), keep_constant=rep(0,length(names)), a_n=rep(0.01,length(names)), modelform="M", fir=0, der_iden=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1), guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1,"loglin_max"=1,"lin_method"="uniform","loglin_method"="uniform",stata=FALSE),Strat_Col='cell'){
    if ("verbose" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$verbose <- FALSE
    }
    if ("guess_constant" %in% names(guesses_control)){
        if (length(guesses_control$guess_constant)<length(a_n)){
            guesses_control$guess_constant <- c(guesses_control$guess_constant, rep(0,length(a_n)-length(guesses_control$guess_constant)))
        }
    } else {
        guesses_control$guess_constant <- rep(0,length(a_n))
    }
    a_n_default <- rep(0,length(a_n))
    for (i in 1:length(a_n)){
        a_n_default[i] = a_n[i]
    }
    if (guesses_control$stata==FALSE){
        setkeyv(df, c(time2, event))
        dfend <- df[get(event)==1, ]
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
        if (length(tu)==0){
            if (guesses_control$verbose){
                print("no events")
            }
            stop()
        }
        if (guesses_control$verbose){
            print(paste(length(tu)," risk groups",sep=""))
        }
        all_names <- unique(names)
        #
        df <- Replace_Missing(df,all_names,0.0,control$verbose)
        #
        dfc <- match(names,all_names)

        term_tot <- max(Term_n)+1
        x_all=as.matrix(df[,all_names, with = FALSE])
        ce <- c(time1,time2,event)
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        if (length(a_n)<length(names)){
            print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
            a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
        } else if (length(a_n)>length(names)){
            print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(Term_n)<length(names)){
            print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        } else if (length(Term_n)>length(names)){
            print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(tform)<length(names)){
            print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
            stop()
        } else if (length(tform)>length(names)){
            print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(keep_constant)<length(names)){
            keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
        } else if (length(keep_constant)>length(names)){
            keep_constant <- keep_constant[1:length(names)]
        }
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$Iterations
        control <- Def_Control(control)
        #
        df_res <- data.table()
        e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
#        print(e$beta_0)
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$LogLik]
        for (it in 1:guesses_control$guesses){
            for (i in 1:length(tform)){
                if (guesses_control$guess_constant[i]==0){
                    if (grepl("log",tform[i],fixed=FALSE)){
                        if (guesses_control$loglin_method == "uniform"){
                            a_n[i] <- runif(1,min=guesses_control$loglin_min,max=guesses_control$loglin_max) + a_n_default[i]
                        } else {
                            print("bad")
                            stop()
                        }
                    } else {
                        if (guesses_control$lin_method == "uniform"){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print("bad")
                            stop()
                        }
                    }
                } else {
                    a_n[i] <- a_n_default[i]
                }
            }
#            print(a_n)
            df_res0 <- data.table()
            e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$LogLik]
            if (is.na(e$LogLik)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            } 
            
        }
        setnames(df_res,c(names,"Deviation"))
        setkeyv(df_res, get("Deviation"))
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
            
        }
#        stop()
        a_n_ind <- which.min(df_res[,get("Deviation")])
#        print(a_n_ind)
#        stop()
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
        #
        control$maxiter <- iteration0
        a_n0 <- copy(a_n)
        control <- Def_Control(control)
        e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
        a_n <- a_n0
        ;
        return (e)
    } else {
        setkeyv(df, c(time2, event, Strat_Col))
        dfend <- df[get(event)==1, ]
        #
        ce <- c(time1,time2,event,Strat_Col)
        all_names <- unique(names)
        #
        df <- Replace_Missing(df,all_names,0.0,control$verbose)
        #
        dfc <- match(names,all_names)

        term_tot <- max(Term_n)+1
        x_all=as.matrix(df[,all_names, with = FALSE])
        #
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
        if (length(tu)==0){
            if (guesses_control$verbose){
                print("no events")
            }
            stop()
        }
        if (guesses_control$verbose){
            print(paste(length(tu)," risk groups",sep=""))
        }
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        #
        for (i in 1:length(uniq)){
            df0 <- dfend[get(Strat_Col)==uniq[i],]
            tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
            if (length(tu0)==0){
                if (control$verbose){
                    print(paste("no events for strata group:",uniq[i],sep=" "))
                }
                stop()
            }
        }
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        control <- Def_Control(control)
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$Iterations
        #
        df_res <- data.table()
        e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
#        print(e$beta_0)
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$LogLik]
        for (it in 1:guesses_control$guesses){
            for (i in 1:length(tform)){
                if (guesses_control$guess_constant[i]==0){
                    if (grepl("log",tform[i],fixed=FALSE)){
                        if (guesses_control$loglin_method == "uniform"){
                            a_n[i] <- runif(1,min=guesses_control$loglin_min,max=guesses_control$loglin_max) + a_n_default[i]
                        } else {
                            print("bad")
                            stop()
                        }
                    } else {
                        if (guesses_control$lin_method == "uniform"){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print("bad")
                            stop()
                        }
                    }
                } else {
                    a_n[i] <- a_n_default[i]
                }
            }
            df_res0 <- data.table()
            e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$LogLik]
            if (is.na(e$LogLik)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            } 
            
        }
        setnames(df_res,c(names,"Deviation"))
        setkeyv(df_res, get("Deviation"))
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        #
        a_n_ind <- which.min(df_res[,get("Deviation")])
#        print(a_n_ind)
#        stop()
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
        e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
        return (e)
    }
    return (NULL)
}

#' Performs basic Cox Proportional Hazards regression with strata effect
#' \code{RunCoxRegression_STRATA} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
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
#' @param Strat_Col column to stratify by
#'
#' @return returns a list of the final results
#' @export
#'
RunCoxRegression_STRATA <- function(df, time1="age_start", time2="age_exit", event="cases",  names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), keep_constant=rep(0,length(names)), a_n=c(0.01), modelform="M", fir=0, der_iden=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1), Strat_Col="cell"){
    setkeyv(df, c(time2, event, Strat_Col))
    dfend <- df[get(event)==1, ]
    #
    ce <- c(time1,time2,event,Strat_Col)
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    #
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
    #
    for (i in 1:length(uniq)){
        df0 <- dfend[get(Strat_Col)==uniq[i],]
        tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
        if (length(tu0)==0){
            if (control$verbose){
                print(paste("no events for strata group:",uniq[i],sep=" "))
            }
            stop()
        }
    }
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
    return (e)
}

#' Calculates proportional hazard for a reference vector
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to calculate risk for a reference
#'
#' @param df data used to calculate PH
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
#' @param control list of parameters for verbosity and tie method
#'
#' @return returns a list of the final results
#' @export
#'
Cox_Relative_Risk <- function(df, time1="age_start", time2="age_exit", event="cases",  names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), keep_constant=rep(0,length(names)), a_n=rep(0.01,length(names)), modelform="M", fir=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)){
    setkeyv(df, c(time2, event))
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    control <- Def_Control(control)
    e <- cox_ph_risk_sub(Term_n, tform, a_n, dfc, x_all,  fir, modelform, control, term_tot)
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with the null model
#' \code{RunCoxRegression} uses user provided data and time/event columns to calculate the log-likelihood with constant hazard ratio
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event column used for event status
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
RunCoxNull <- function(df, time1="age_start", time2="age_exit", event="cases",control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)){
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    ce <- c(time1,time2,event)
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    control <- Def_Control(control)
    e <- cox_ph_null( control, as.matrix(df[,ce, with = FALSE]), tu)
    return (e)

}


#' Performs Cox Proportional Hazard model plots
#' \code{RunCoxPlots} uses user provided data, time/event columns, vectors specifying the model, and options to choose plots and saves plots
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
#' @param control list of parameters controlling the convergence
#' @param plot_options list of parameters controlling the plot options
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#'
RunCoxPlots <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options){
    if (plot_options$verbose){
        print("Starting Plot Function")
    }
    setkeyv(df, c(time2, event))
    base  <- NULL
    der_iden <- 0
    Plot_Type <- plot_options$type
    if (plot_options$verbose){
        print("Getting Plot Info")
    }
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (plot_options$verbose){
            print("no events")
        }
        stop()
    }
    if (plot_options$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    if ("type" %in% names(plot_options)){
        ;
    } else {
        if (plot_options$verbose){
            print("Plot type not given")
        }
        stop()
    }
    if ("age_unit" %in% names(plot_options)){
        ;
    } else {
        plot_options$age_unit <- "unitless"
    }
    if ("strat_haz" %in% names(plot_options)){
        if (plot_options$strat_haz){
            if ("Strat_Col" %in% names(plot_options)){
                if (plot_options$Strat_Col %in% names(df)){
                    ;
                } else {
                    if (plot_options$verbose){
                        print("Stratification Column not in the dataframe")
                    }
                    stop()
                }
            } else {
                if (plot_options$verbose){
                    print("Stratification Column not given")
                }
                stop()
            }
        }
    } else {
        plot_options$strat_haz <- FALSE
    }
    if ("smooth_haz" %in% names(plot_options)){
        if (plot_options$smooth_haz){
            if ("bw" %in% names(plot_options)){
                if (plot_options$bw > 0){
                    ;
                } else {
                    if (plot_options$verbose){
                        print("bandwidth needs to be above zero")
                    }
                    stop()
                }
            } else {
                if (plot_options$verbose){
                    print("bandwidth not given")
                }
                stop()
            }
        }
    } else {
        plot_options$smooth_haz <- FALSE
    }
    if ("Martingale" %in% names(plot_options)){
        if (plot_options$Martingale){
            if ("cov_cols" %in% names(plot_options)){
                for (cov_i in 1:length(plot_options$cov_cols)){
                    dose_col <- unlist(plot_options$cov_cols,use.names=FALSE)[cov_i]
                    if (dose_col%in% names(df)){
                        ;
                    } else {
                        if (plot_options$verbose){
                            print("Covariate column "+dose_col+" is not in the dataframe")
                        }
                        stop()
                    }
                }
            } else {
                if (plot_options$verbose){
                    print("dose column not given")
                }
                stop()
            }
        }
    } else {
        plot_options$Martingale <- FALSE
    }
    if ("KM" %in% names(plot_options)){
        if (plot_options$KM){
            if ("studyID" %in%  names(plot_options)){
                if (plot_options$studyID%in% names(df)){
                    ;
                } else {
                    if (plot_options$verbose){
                        print("ID column is not in the dataframe")
                    }
                    stop()
                }
            } else {
                if (plot_options$verbose){
                    print("ID column not given")
                }
                stop()
            }
        }
    }
    if (Plot_Type[1]=="SURV"){
        if ("time_lims" %in% names(plot_options)){
            ;
        } else {
            plot_options$time_lims <- c(min(tu),max(tu))
        }
    }
    for (iden_col in c("verbose","Martingale","surv_curv","strat_haz","smooth_haz","KM")){
        if (iden_col %in% names(plot_options)){
            ;
        } else {
            plot_options[iden_col] <- FALSE
        }
    }
    control <- Def_Control(control)
    verbose <- copy(plot_options$verbose)
    verbosec <- copy(control$verbose)
    maxiterc <- copy(control$maxiter)
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all=as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    time1 <- ce[1]
    time2 <- ce[2]
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    control$maxiter <- -1
#    control$verbose <- verbose
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
#    control$verbose <- verbosec
    control$maxiter <- maxiterc
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
    if (Plot_Type[1]=="SURV"){
        if (verbose){
            print("starting ph_plot")
        }
        #
        e <- cox_ph_plot(Term_n, tform, a_n,er, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        #
        t <- c()
        h <- c()
        ch <- c()
        surv <- c()
        dt <- 1
        if (verbose){
            print("writing survival data")
        }
        dft=data.table("time"=tu,"base"=e$baseline,"basehaz"=e$standard_error)
        for (i in tu){
            t <- c(t,i)
            temp <- sum(dft[time<i, base])
            ch <- c(ch, temp)
            if (length(h)==0){
                h <- c(temp)
            } else {
                h <- c(h, ch[length(ch)]-ch[length(ch)-1])
            }
            surv <- c(surv, exp(-1*temp))
        }
        #
        age_unit <- plot_options$age_unit
        if (plot_options$Martingale==TRUE){
            #
            CoxMartingale(verbose, df, time1, time2, event, e, t, ch, plot_options$cov_cols, Plot_Type[2], age_unit,plot_options$studyID)
            #
        }
        if (plot_options$surv_curv==TRUE){
            CoxSurvival(t,h,ch,surv,Plot_Type[2],verbose,plot_options$time_lims, age_unit)
            if (plot_options$strat_haz==TRUE){
                CoxStratifiedSurvival(verbose, df, event, time1, time2, all_names,Term_n, tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type, plot_options$Strat_Col,plot_options$time_lims,age_unit)
            }
            if (plot_options$smooth_haz==TRUE){
                CoxSmoothHazard(dft,Plot_Type[2],verbose,plot_options$bw,plot_options$time_lims,age_unit)
            }
        }
        if (plot_options$KM==TRUE){
            #
            CoxKaplanMeier(verbose, verbosec, plot_options$studyID,all_names,df,event,time1,time2,tu,Term_n, tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type,age_unit)
        }
    } else if (Plot_Type[1]=="RISK"){
        CoxRisk(verbose, df, event, time1, time2, names,Term_n, tform, a_n, fir, der_iden, modelform, control,keep_constant, Plot_Type, b, er)
        #      (verbose, df, event, time1, time2, names,    Term_n, tform, a_n, fir, der_iden, modelform, control,keep_constant, Plot_Type, b, er)
    } else if (Plot_Type[1]=="SCHOENFELD"){
        age_unit <- plot_options$age_unit
        PlotCox_Schoenfeld_Residual(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,age_unit,Plot_Type[2])
    }
    ;
    return ("Passed")
}

#' Performs basic cox regression, with multiple guesses, starts with a single term
#' \code{RunCoxRegression_Tier_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
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
#' @param guesses_control list of parameters to control how the guessing works
#' @param Strat_Col column to stratify by if needed
#'
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data
RunCoxRegression_Tier_Guesses <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col){
    if ("verbose" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$verbose <- FALSE
    }
    if ("guess_constant" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$guess_constant <- rep(0,length(a_n))
    }
    if ("guesses_start" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$guesses_start <- guesses_control$guesses
    }
    t_initial <- guesses_control$term_initial
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    name_initial <- c()
    term_n_initial <- c()
    tform_initial <- c()
    constant_initial <- c()
    a_n_initial <- c()
    guess_constant <- c()
    #
    for (i in 1:length(a_n)){
        if (Term_n[i] %in% t_initial){
            name_initial <- c(name_initial, names[i])
            term_n_initial <- c(term_n_initial, Term_n[i])
            tform_initial <- c(tform_initial, tform[i])
            constant_initial <- c(constant_initial, keep_constant[i])
            a_n_initial <- c(a_n_initial, a_n[i])
            guess_constant <- c(guess_constant, 0)
        }
    }
    guesses_control$guess_constant <- guess_constant
    guess_second <- guesses_control$guesses
    guesses_control$guesses <- guesses_control$guesses_start
    e <- RunCoxRegression_Guesses(df, time1, time2, event, name_initial, term_n_initial, tform_initial, constant_initial, a_n_initial, modelform, fir, der_iden, control, guesses_control,Strat_Col)
    #
    a_n_initial <- e$beta_0
    guess_constant <- c()
    j <- 0
    for (i in 1:length(a_n)){
        if (Term_n[i] %in% t_initial){
            a_n[i] <- a_n_initial[j]
            j = j+1
            guess_constant <- c(guess_constant, 1)
        } else {
            guess_constant <- c(guess_constant, 0)
        }
    }
    guesses_control$guess_constant <- guess_constant
    guesses_control$guesses <- guess_second
    e <- RunCoxRegression_Guesses(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col)
    #
    return(e)
}

#' Calculates Schoenfeld residuals for a Cox Proportional Hazards regression
#' \code{RunCox_Schoenfeld_Residual} uses user provided data, time/event columns, vectors specifying the model, and options to calculate the residuals
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
#' @return returns a list of the final results
#' @export
#'
#' @importFrom rlang .data

RunCox_Schoenfeld_Residual <- function(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    setkeyv(df, c(time2, event))
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
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
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(Term_n)<length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(Term_n)>length(names)){
        print(paste("Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(tform)<length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    } else if (length(tform)>length(names)){
        print(paste("Term types used: ",length(tform),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[1:length(names)]
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_schoenfeld_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    ;
    return (e)
}
