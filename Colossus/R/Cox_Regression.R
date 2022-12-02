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
        a_n <- c(a_n, rep(0.01,length(a_n)-length(names)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    #
    a_n0 <- copy(a_n)
    control <- Def_Control(control)
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
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
    if (guesses_control$stata==FALSE){
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
        e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
#        print(e$beta_0)
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$LogLik]
        for (it in 1:guesses_control$guesses){
            for (i in 1:length(tform)){
                if ("log" %in% tform[i]){
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
            e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$LogLik]
            df_res <- rbindlist(list(df_res, df_res0)) 
            
        }
        print(df_res)
#        stop()
        a_n_ind <- which.max(df_res[,get(paste(length(e$beta_0)+1))])
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
                if ("log" %in% tform[i]){
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
            e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
            }
            df_res0[,paste(length(e$beta_0)+1):=e$LogLik]
            df_res <- rbindlist(list(df_res, df_res0)) 
            
        }
        print(df_res)
        #
        a_n_ind <- which.max(df_res[,get(paste(length(e$beta_0)+1))])
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
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    control <- Def_Control(control)
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
RunCoxPlots <- function(df, time1="age_start", time2="age_exit", event="cases", names=c("dose"), Term_n=rep(0,length(names)), tform=rep("loglin",length(names)), keep_constant=rep(0,length(names)), a_n=rep(0.01,length(names)), modelform="M", fir=0, control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1), plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=TRUE)){
    setkeyv(df, c(time2, event))
    base  <- NULL
    der_iden <- 0
    Plot_Type <- plot_options$type
    dfend <- df[get(event)==1, ]
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
    if ("type" %in% names(plot_options)){
        ;
    } else {
        print("Plot type not given")
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
                    print("Stratification Column not in the dataframe")
                    stop()
                }
            } else {
                print("Stratification Column not given")
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
                    print("bandwidth needs to be above zero")
                    stop()
                }
            } else {
                print("bandwidth not given")
                stop()
            }
        }
    } else {
        plot_options$smooth_haz <- FALSE
    }
    if ("Martingale" %in% names(plot_options)){
        if (plot_options$Martingale){
            if ("dose_col" %in% names(plot_options)){
                if (plot_options$dose_col%in% names(df)){
                    ;
                } else {
                    print("Dose column is not in the dataframe")
                    stop()
                }
            } else {
                print("dose column not given")
                stop()
            }
        }
    } else {
        plot_options$Martingale <- FALSE
    }
    if (Plot_Type[1]=="SURV"){
        if ("studyID" %in%  names(plot_options)){
            if (plot_options$studyID%in% names(df)){
                ;
            } else {
                print("ID column is not in the dataframe")
                stop()
            }
        } else {
            print("ID column not given")
            stop()
        }
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
    verbose <- copy(plot_options$verbose)
    verbosec <- copy(control$verbose)
    dfend <- df[get(event)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
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
    time1 <- ce[1]
    time2 <- ce[2]
    #
#    control$verbose <- verbose
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
#    control$verbose <- verbosec
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
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
    if (Plot_Type[1]=="SURV"){
        if (plot_options$Martingale==TRUE){
            #
            CoxMartingale(verbose, df, time1, time2, event, e, t, ch, plot_options$dose_col, Plot_Type[2], age_unit)
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
        CoxRisk(verbose, df, event, time1, time2, all_names,Term_n, tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type)
    }
    ;
    return ("Passed")
}
