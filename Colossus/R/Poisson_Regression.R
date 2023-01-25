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

#' Performs basic poisson calculation
#' \code{RunPoissonRegression_Single} uses user provided data, person-year/event columns, vectors specifying the model, and returns the score
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
#' @param event column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the intial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#'
#' @return prints the final results and return null
#' @export
#'
RunPoissonRegression_Single <- function(df, pyr, event, names, Term_n, tform, a_n, modelform, fir, control){
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
    e <- poisson_transition_single(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir, modelform, control,term_tot)
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
    if ("guess_constant" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$guess_constant <- rep(0,length(a_n))
    }
    a_n_default <- rep(0,length(a_n))
    for (i in 1:length(a_n)){
        a_n_default[i] = a_n[i]
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
        a_n0 <- rep(0,length(a_n))
        for (i in 1:length(a_n)){
            a_n0[i] = a_n[i]
        }
        df_res <- data.table()
        e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
#        print(e$beta_0)
        radius <- 0
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
            radius = radius + (e$beta_0[i]-a_n0[i])^2
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$Deviation]
        df_res[,paste(length(e$beta_0)+2):=sqrt(radius)]
        df_res[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
        df_res[,paste(length(e$beta_0)+4):=e$Converged]
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
            a_n0 <- rep(0,length(a_n))
            for (i in 1:length(a_n)){
                a_n0[i] = a_n[i]
            }
            df_res0 <- data.table()
            e <- poisson_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot)
            radius <- 0
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
                radius = radius + (e$beta_0[i]-a_n0[i])^2
            }
            df_res0[,paste(length(e$beta_0)+1):=e$Deviation]
            df_res0[,paste(length(e$beta_0)+2):=sqrt(radius)]
            df_res0[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
            df_res0[,paste(length(e$beta_0)+4):=e$Converged]
            if (is.na(e$Deviation)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            }
            
        }
        setnames(df_res,c(names,"Deviation","Radius","Iterations","Converged"))
        setkeyv(df_res, "Deviation")
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        a_n_ind <- which.min(df_res[,get("Deviation")])
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
        if (length(guesses_control$guess_constant)<length(a_n)){
            guesses_control$guess_constant <- c(guesses_control$guess_constant, rep(0,length(a_n)-length(guesses_control$guess_constant)))
        }
        a_n_default <- rep(0,length(a_n))
        for (i in 1:length(a_n)){
            a_n_default[i] = a_n[i]
        }
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
        a_n0 <- rep(0,length(a_n))
        for (i in 1:length(a_n)){
            a_n0[i] = a_n[i]
        }
        e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
#        print(e$beta_0)
        radius <- 0
        for (i in 1:length(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
            radius = radius + (e$beta_0[i]-a_n0[i])^2
#            print(df_res)
        }
        df_res[,paste(length(e$beta_0)+1):=e$Deviation]
        df_res[,paste(length(e$beta_0)+2):=sqrt(radius)]
        df_res[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
        df_res[,paste(length(e$beta_0)+4):=e$Converged]
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
            a_n0 <- rep(0,length(a_n))
            for (i in 1:length(a_n)){
                a_n0[i] = a_n[i]
            }
            df_res0 <- data.table()
            e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
            radius <- 0
            for (i in 1:length(e$beta_0)){
                df_res0[,paste(i):=e$beta_0[i]]
                radius = radius + (e$beta_0[i]-a_n0[i])^2
            }
            df_res0[,paste(length(e$beta_0)+1):=e$Deviation]
            df_res0[,paste(length(e$beta_0)+2):=sqrt(radius)]
            df_res0[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
            df_res0[,paste(length(e$beta_0)+4):=e$Converged]
            if (is.na(e$Deviation)){
                ;
            } else {
                df_res <- rbindlist(list(df_res, df_res0)) 
            }
            
        }
        setnames(df_res,c(names,"Deviation","Radius","Iterations","Converged"))
        setkeyv(df_res, "Deviation")
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        #
        a_n_ind <- which.min(df_res[,get("Deviation")])
#        print(a_n_ind)
#        stop()
        control$keep_strata <- keep_str
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
        e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
        return (e)
    }
    return (NULL)
}

#' Performs basic poisson regression, with multiple guesses, starts with a single term
#' \code{RunPoissonRegression_Tier_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
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
RunPoissonRegression_Tier_Guesses <- function(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols, guesses_control){
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
    e <- RunPoissonRegression_Guesses(df, pyr, event, name_initial, term_n_initial, tform_initial, constant_initial, a_n_initial, modelform, fir, der_iden, control, Strat_Cols, guesses_control)
    if (guesses_control$verbose){
        print("INITIAL TERM COMPLETE")
        print(e)
    }
    #
    a_n_initial <- unlist(e$beta_0,use.names=FALSE)
    guess_constant <- c()
    j <- 1
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
    e <- RunPoissonRegression_Guesses(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols, guesses_control)
    #
    return(e)
}
    
#' Performs basic poisson regression with stratification, starts with regression on the stratification values
#' \code{RunPoissonRegression_Strata_First} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
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
RunPoissonRegression_Strata_First <- function(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols, guesses_control){
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            ;
        } else {
            df$CONST <- 1
        }
    }
    if ("verbose" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$verbose <- FALSE
    }
    if ("guesses_start" %in% names(guesses_control)){
        ;
    } else {
        guesses_control$guesses_start <- guesses_control$guesses
    }
    #
    val <- factorize(df, Strat_Cols,guesses_control$verbose)
    df <- val$df
    names_strata <- c( val$cols)
    Term_max <- 0
    Term_n_strata <- rep(Term_max,length(val$cols))
    tform_strata <- rep("loglin",length(val$cols))
    keep_constant_strata <- rep(0,length(val$cols))
    a_n_strata <- runif(length(val$cols),min=-0.01,0.01)
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
    #
    guess_constant <- rep(0,length(a_n_strata))
    #
    guesses_control$guess_constant <- guess_constant
    guesses_control$stata=FALSE
    #
    guess_second <- guesses_control$guesses
    guesses_control$guesses <- guesses_control$guesses_start
    e <- RunPoissonRegression_Guesses(df, pyr, event, names_strata, Term_n_strata, tform_strata, keep_constant_strata, a_n_strata, modelform, fir, der_iden, control, Strat_Cols, guesses_control)
    if (guesses_control$verbose){
        print("INITIAL TERM COMPLETE")
        print(e)
    }
    guesses_control$guesses <- guess_second
    #
    names <- c(names,names_strata)
    Term_n <- c(Term_n,Term_n_strata)
    tform <- c(tform,tform_strata)
    keep_constant <- c(keep_constant,keep_constant_strata)
    a_n <- c(a_n,a_n_strata)
    guesses_control$guess_constant <- c(guesses_control$guess_constant,rep(0,length(names)-length(names_strata)))
    #
    a_n_default <- rep(0,length(a_n))
    for (i in 1:length(a_n)){
        a_n_default[i] = a_n[i]
    }
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
    a_n0 <- rep(0,length(a_n))
    for (i in 1:length(a_n)){
        a_n0[i] = a_n[i]
    }
    e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
#        print(e$beta_0)
    radius <- 0
    for (i in 1:length(e$beta_0)){
        df_res[,paste(i):=e$beta_0[i]]
        radius = radius + (e$beta_0[i]-a_n0[i])^2
#            print(df_res)
    }
    df_res[,paste(length(e$beta_0)+1):=e$Deviation]
    df_res[,paste(length(e$beta_0)+2):=sqrt(radius)]
    df_res[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
    df_res[,paste(length(e$beta_0)+4):=e$Converged]
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
        a_n0 <- rep(0,length(a_n))
        for (i in 1:length(a_n)){
            a_n0[i] = a_n[i]
        }
        df_res0 <- data.table()
        e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
        radius <- 0
        for (i in 1:length(e$beta_0)){
            df_res0[,paste(i):=e$beta_0[i]]
            radius = radius + (e$beta_0[i]-a_n0[i])^2
        }
        df_res0[,paste(length(e$beta_0)+1):=e$Deviation]
        df_res0[,paste(length(e$beta_0)+2):=sqrt(radius)]
        df_res0[,paste(length(e$beta_0)+3):=e$Control_List$Iteration]
        df_res0[,paste(length(e$beta_0)+4):=e$Converged]
        if (is.na(e$Deviation)){
            ;
        } else {
            df_res <- rbindlist(list(df_res, df_res0)) 
        }
        
    }
    setnames(df_res,c(names,"Deviation","Radius","Iterations","Converged"))
    setkeyv(df_res, "Deviation")
    if (guesses_control$verbose){
        print(df_res)
        fwrite(df_res,"last_guess_full.csv")
    }
    #
    a_n_ind <- which.min(df_res[,get("Deviation")])
#        print(a_n_ind)
#        stop()
    control$keep_strata <- keep_str
    a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[1:length(a_n)]
    e <- poisson_strata_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,keep_constant,term_tot,rep(1,length(val$cols)))
    return (e)
}
    
    
