#' Performs basic Cox Proportional Hazards regression
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(0.1, 0.1, 0.1, 0.1)
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform,
#'                      keep_constant, a_n, modelform, fir, der_iden, control)
#' @importFrom rlang .data

RunCoxRegression <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
            ", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    a_n0 <- copy(a_n)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    a_n <- a_n0
    return (e)
}

#' Performs basic Cox Proportional Hazards calculation with no derivative
#' \code{RunCoxRegression_Single} uses user provided data, time/event columns, vectors specifying the model, and options and returns the log-likelihood
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' 
#' keep_constant <- c(0,0,0,0)
#' 
#' control=list("Ncores"=2,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression_Single(df, time1, time2, event, names, Term_n, tform,
#'                              keep_constant, a_n, modelform, fir, control)
#'
#' @importFrom rlang .data

RunCoxRegression_Single <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, control){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
            ", Remaining filled with 0.01",sep=""))
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
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition_single(Term_n,tform,a_n,dfc,x_all, fir, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    a_n <- a_n0
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with a multiplicative log-linear model
#' \code{RunCoxRegression_Basic} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,
#'    'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression_Basic(df, time1, time2, event, names, keep_constant, a_n, der_iden, control)
#'
#' @importFrom rlang .data

RunCoxRegression_Basic <- function(df, time1, time2, event0, names, keep_constant, a_n, der_iden, control){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
            ", Remaining filled with 0.01",sep=""))
        a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
    } else if (length(a_n)>length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),sep=""))
        stop()
    }
    if (length(keep_constant)<length(names)){
        keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
    } else if (length(keep_constant)>length(names)){
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    a_n0 <- copy(a_n)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition_basic(a_n,dfc,x_all, der_iden, control,
         as.matrix(df[,ce, with = FALSE]),tu,keep_constant)
    a_n <- a_n0
    return (e)
}

#' Performs basic Cox Proportional Hazards regression, Allows for multiple starting guesses
#' \code{RunCoxRegression_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions. Has additional options to starting with several initial guesses
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param guesses_control list of parameters to control how the guessing works
#' @param Strat_Col column to stratify by if needed
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1),
#'                       "e"=c(0,   0,   1,   0,   0,   0,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001,"lin_max"=1,
#'     "loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform","loglin_method"="uniform",strata=TRUE)
#' Strat_Col='e'
#' 
#' e <- RunCoxRegression_Guesses(df, time1, time2, event, names, Term_n,
#'                               tform, keep_constant, a_n, modelform, fir,
#'                               der_iden, control,guesses_control,Strat_Col)
#' @importFrom rlang .data
RunCoxRegression_Guesses <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col){
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    a_n_default <- rep(0,length(a_n))
    for (i in seq_along(a_n)){
        a_n_default[i] <- a_n[i]
    }
    if (guesses_control$strata==FALSE){
        setkeyv(df, c(time2, event0))
        dfend <- df[get(event0)==1, ]
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
        if (length(tu)==0){
            print("no events")
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
        x_all <- as.matrix(df[,all_names, with = FALSE])
        ce <- c(time1,time2,event0)
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        if (length(a_n)<length(names)){
            print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
                ", Remaining filled with 0.01",sep=""))
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
            keep_constant <- keep_constant[seq_len(length(names))]
        }
        rmin <- guesses_control$rmin
        rmax <- guesses_control$rmax
        if (length(rmin)!=length(rmax)){
            if (control$verbose){
                print("rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
            }
        }
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$maxiter
        control <- Def_Control(control)
        #
        df_res <- data.table()
        e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
        for (i in seq_along(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
        }
        df_res[,paste(length(e$beta_0)+1):=-2*e$LogLik]
        for (it in 1:guesses_control$guesses){
            if (length(rmin)!=length(rmax)){
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        if (grepl("_int",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$intercept_min,max=guesses_control$intercept_max) + a_n_default[i]
                        } else if (grepl("lin_exp_exp_slope",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$exp_slope_min,max=guesses_control$exp_slope_max) + a_n_default[i]
                        } else if (grepl("_slope",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else if (grepl("loglin",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$exp_min,max=guesses_control$exp_max) + a_n_default[i]
                        } else if (grepl("lin",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print(paste("tform not implemented ", tform[i],sep=""))
                            stop()
                        }
                    } else {
                        a_n[i] <- a_n_default[i]
                    }
                }
            } else {
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        a_n[i] <- runif(1,min=guesses_control$rmin[i],max=guesses_control$rmax[i]) + a_n_default[i]
                    } else {
                        a_n[i] <- a_n_default[i]
                    }
                }
            }
            e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
                as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
            if (is.na(e$LogLik)){
                #fine
            } else {
                df_res0 <- data.table()
                for (i in seq_along(e$beta_0)){
                    df_res0[,paste(i):=e$beta_0[i]]
                }
                df_res0[,paste(length(e$beta_0)+1):=-2*e$LogLik]
                df_res <- rbindlist(list(df_res, df_res0)) 
            } 
            
        }
        setnames(df_res,c(names,"Deviation"))
        setkeyv(df_res, "Deviation")
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
            
        }
        a_n_ind <- which.min(df_res[,get("Deviation")])
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[seq_along(a_n)]
        #
        control$maxiter <- iteration0
        a_n0 <- copy(a_n)
        control <- Def_Control(control)
        e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
        a_n <- a_n0
        #fine
        return (e)
    } else {
        dfend <- df[get(event0)==1, ]
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        #
        for (i in seq_along(uniq)){
            df0 <- dfend[get(Strat_Col)==uniq[i],]
            tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
            if (length(tu0)==0){
                if (control$verbose){
                    print(paste("no events for strata group:",uniq[i],sep=" "))
                }
                df <- df[get(Strat_Col)!=uniq[i],]
            }
        }
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        if (control$verbose){
            print(paste(length(uniq)," strata used",sep=" "))
        }
        setkeyv(df, c(time2, event0, Strat_Col))
        dfend <- df[get(event0)==1, ]
        #
        ce <- c(time1,time2,event0,Strat_Col)
        all_names <- unique(names)
        #
        df <- Replace_Missing(df,all_names,0.0,control$verbose)
        #
        dfc <- match(names,all_names)

        term_tot <- max(Term_n)+1
        x_all <- as.matrix(df[,all_names, with = FALSE])
        #
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
        if (length(tu)==0){
            print("no events")
            stop()
        }
        if (guesses_control$verbose){
            print(paste(length(tu)," risk groups",sep=""))
        }
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        control <- Def_Control(control)
        #
        if (length(a_n)<length(names)){
            print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
                ", Remaining filled with 0.01",sep=""))
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
            keep_constant <- keep_constant[seq_len(length(names))]
        }
        rmin <- guesses_control$rmin
        rmax <- guesses_control$rmax
        if (length(rmin)!=length(rmax)){
            if (control$verbose){
                print("rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
            }
        }
        #
        iteration0 <- control$maxiter
        control$maxiter <- guesses_control$maxiter
        #
        df_res <- data.table()
        e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
        for (i in seq_along(e$beta_0)){
            df_res[,paste(i):=e$beta_0[i]]
        }
        df_res[,paste(length(e$beta_0)+1):=-2*e$LogLik]
        for (it in 1:guesses_control$guesses){
            if (length(rmin)!=length(rmax)){
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        if (grepl("_int",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$intercept_min,max=guesses_control$intercept_max) + a_n_default[i]
                        } else if (grepl("lin_exp_exp_slope",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$exp_slope_min,max=guesses_control$exp_slope_max) + a_n_default[i]
                        } else if (grepl("_slope",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else if (grepl("loglin",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$exp_min,max=guesses_control$exp_max) + a_n_default[i]
                        } else if (grepl("lin",tform[i],fixed=FALSE)){
                            a_n[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print(paste("tform not implemented ", tform[i],sep=""))
                            stop()
                        }
                    } else {
                        a_n[i] <- a_n_default[i]
                    }
                }
            } else {
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        a_n[i] <- runif(1,min=guesses_control$rmin[i],max=guesses_control$rmax[i]) + a_n_default[i]
                    } else {
                        a_n[i] <- a_n_default[i]
                    }
                }
            }
            e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
                as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
            if (is.na(e$LogLik)){
                #fine
            } else {
                df_res0 <- data.table()
                for (i in seq_along(e$beta_0)){
                    df_res0[,paste(i):=e$beta_0[i]]
                }
                df_res0[,paste(length(e$beta_0)+1):=-2*e$LogLik]
                df_res <- rbindlist(list(df_res, df_res0)) 
            } 
            
        }
        setnames(df_res,c(names,"Deviation"))
        setkeyv(df_res, "Deviation")
        if (guesses_control$verbose){
            print(df_res)
            fwrite(df_res,"last_guess.csv")
        }
        #
        a_n_ind <- which.min(df_res[,get("Deviation")])
        a_n <- unlist(df_res[a_n_ind],use.names=FALSE)[seq_along(a_n)]
        e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
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
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param Strat_Col column to stratify by
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1),
#'                       "e"=c(0,   0,   0,   0,   1,   0,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' Strat_Col='e'
#' 
#' e <- RunCoxRegression_STRATA(df, time1, time2, event, names, Term_n, tform, keep_constant,
#'                              a_n, modelform, fir, der_iden, control,Strat_Col)
#'
RunCoxRegression_STRATA <- function(df, time1, time2, event0,  names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Col){
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
    #
    for (i in seq_along(uniq)){
        df0 <- dfend[get(Strat_Col)==uniq[i],]
        tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
        if (length(tu0)==0){
            if (control$verbose){
                print(paste("no events for strata group:",uniq[i],sep=" "))
            }
            df <- df[get(Strat_Col)!=uniq[i],]
        }
    }
    uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
    if (control$verbose){
        print(paste(length(uniq)," strata used",sep=" "))
    }
    setkeyv(df, c(time2, event0, Strat_Col))
    #
    ce <- c(time1,time2,event0,Strat_Col)
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
#    print(uniq)
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",
            length(names),", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_STRATA(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
    return (e)
}

#' Calculates proportional hazard for a reference vector
#' \code{RunCoxRegression} uses user provided data, time/event columns, vectors specifying the model, and options to calculate risk for a reference
#'
#' @param df data used to calculate PH
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters for verbosity and tie method
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' fir <- 0
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' a_n <- c(1.1, 0.1, 0.2, 0.5) #used to test at a specific point
#' 
#' keep_constant <- c(0,0,0,0)
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n, tform,
#'      keep_constant, a_n, modelform, fir, control)
#'
Cox_Relative_Risk <- function(df, time1, time2, event0,  names, Term_n, tform, keep_constant, a_n, modelform, fir, control){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    if (control$verbose){
        print("Starting")
        print(keep_constant)
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",
            length(names),", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    e <- cox_ph_risk_sub(Term_n, tform, a_n, dfc, x_all,  fir, modelform, control, term_tot)
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with the null model
#' \code{RunCoxRegression} uses user provided data and time/event columns to calculate the log-likelihood with constant hazard ratio
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' 
#' control=list("Ncores"=2,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxNull(df, time1, time2, event, control)
#'
RunCoxNull <- function(df, time1, time2, event0,control){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    ce <- c(time1,time2,event0)
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    e <- cox_ph_null( control, as.matrix(df[,ce, with = FALSE]), tu)
    return (e)

}


#' Performs Cox Proportional Hazard model plots
#' \code{RunCoxPlots} uses user provided data, time/event columns, vectors specifying the model, and options to choose and saves plots
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#' @param plot_options list of parameters controlling the plot options
#'
#' @return saves the plots in the current directory and returns a string that it passed
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(-0.1, 0.5, 1.1, -0.3)
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' #setting maxiter below 0 forces the function to calculate the score and return
#' plot_options=list("type"=c("SURV","run_0"), "studyID"="UserID",'verbose'=FALSE)
#' 
#' RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant,
#'             a_n, modelform, fir, control, plot_options)
#'
RunCoxPlots <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options){
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if (plot_options$verbose){
        print("Starting Plot Function")
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    setkeyv(df, c(time2, event0))
    base  <- NULL
    der_iden <- 0
    Plot_Type <- plot_options$type
    if (plot_options$verbose){
        print("Getting Plot Info")
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (plot_options$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    if ("type" %in% names(plot_options)){
        #fine
    } else {
        if (plot_options$verbose){
            print("Plot type not given")
        }
        stop()
    }
    if ("age_unit" %in% names(plot_options)){
        #fine
    } else {
        plot_options$age_unit <- "unitless"
    }
    if ("strat_haz" %in% names(plot_options)){
        if (plot_options$strat_haz){
            if ("Strat_Col" %in% names(plot_options)){
                if (plot_options$Strat_Col %in% names(df)){
                    #fine
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
    if ("Martingale" %in% names(plot_options)){
        if (plot_options$Martingale){
            if ("cov_cols" %in% names(plot_options)){
                for (cov_i in seq_along(plot_options$cov_cols)){
                    dose_col <- unlist(plot_options$cov_cols,use.names=FALSE)[cov_i]
                    if (dose_col%in% names(df)){
                        #fine
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
                    #fine
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
            #fine
        } else {
            plot_options$time_lims <- c(min(tu),max(tu))
        }
    }
    for (iden_col in c("verbose","Martingale","surv_curv","strat_haz","KM")){
        if (iden_col %in% names(plot_options)){
            #fine
        } else {
            plot_options[iden_col] <- FALSE
        }
    }
    control <- Def_Control(control)
    verbose <- copy(plot_options$verbose)
    verbosec <- copy(control$verbose)
    maxiterc <- copy(control$maxiter)
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    all_names <- unique(names)
    #
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    time1 <- ce[1]
    time2 <- ce[2]
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",
            length(names),", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    #
    control$maxiter <- -1
    e <- cox_ph_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    control$maxiter <- maxiterc
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
    #
    if (Plot_Type[1]=="SURV"){
        if (verbose){
            print("starting ph_plot")
        }
        #
        e <- cox_ph_plot(Term_n, tform, a_n,er, dfc, x_all, fir, der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, Plot_Type , 0)
        #
        t <- c()
        h <- c()
        ch <- c()
        surv <- c()
        if (verbose){
            print("writing survival data")
        }
        dft <- data.table("time"=tu,"base"=e$baseline,"basehaz"=e$standard_error)
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
            CoxMartingale(verbose, df, time1, time2, event0, e, t, ch, plot_options$cov_cols,
                          Plot_Type[2], age_unit,plot_options$studyID)
            #
        }
        if (plot_options$surv_curv==TRUE){
            CoxSurvival(t,h,ch,surv,Plot_Type[2],verbose,plot_options$time_lims, age_unit)
            if (plot_options$strat_haz==TRUE){
                CoxStratifiedSurvival(verbose, df, event0, time1, time2, all_names,Term_n, tform, a_n, er,
                     fir, der_iden, modelform, control,keep_constant, Plot_Type, plot_options$Strat_Col,plot_options$time_lims,age_unit)
            }
        }
        if (plot_options$KM==TRUE){
            #
            CoxKaplanMeier(verbose, verbosec, plot_options$studyID,all_names,df,event0,time1,time2,tu,Term_n,
                 tform, a_n, er, fir, der_iden, modelform, control,keep_constant, Plot_Type,age_unit)
        }
    } else if (Plot_Type[1]=="RISK"){
        CoxRisk(verbose, df, event0, time1, time2, names,Term_n, tform, a_n, fir,
                der_iden, modelform, control,keep_constant, Plot_Type, b, er)
    } else if (Plot_Type[1]=="SCHOENFELD"){
        age_unit <- plot_options$age_unit
        PlotCox_Schoenfeld_Residual(df, time1, time2, event0, names, Term_n, tform, keep_constant,
                                    a_n, modelform, fir, der_iden, control,age_unit,Plot_Type[2])
    }
    return ("Passed")
}

#' Performs basic cox regression, with multiple guesses, starts with solving for a single term
#' \code{RunCoxRegression_Tier_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param guesses_control list of parameters to control how the guessing works
#' @param Strat_Col column to stratify by if needed
#'
#' @return returns a list of the final results
#' @export
#'
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1),
#'                       "e"=c(0,   0,   0,   0,   1,   0,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001,
#'    "lin_max"=1,"loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform",
#'    "loglin_method"="uniform",strata=TRUE,term_initial = c(0,1))
#' Strat_Col='e'
#' 
#' e <- RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, Term_n, tform, keep_constant,
#'                                    a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)
#'
#' @importFrom rlang .data
RunCoxRegression_Tier_Guesses <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col){
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    rmin <- guesses_control$rmin
    rmax <- guesses_control$rmax
    if (length(rmin)!=length(rmax)){
        if (control$verbose){
            print("rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
        }
    }
    #
    name_initial <- c()
    term_n_initial <- c()
    tform_initial <- c()
    constant_initial <- c()
    a_n_initial <- c()
    guess_constant <- c()
    #
    for (i in seq_along(a_n)){
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
    e <- RunCoxRegression_Guesses(df, time1, time2, event0, name_initial, term_n_initial, tform_initial,
         constant_initial, a_n_initial, modelform, fir, der_iden, control, guesses_control,Strat_Col)
    #
    a_n_initial <- e$beta_0
    guess_constant <- c()
    j <- 1
    for (i in seq_along(a_n)){
        if (Term_n[i] %in% t_initial){
            a_n[i] <- a_n_initial[j]
            j <- j+1
            guess_constant <- c(guess_constant, 1)
        } else {
            guess_constant <- c(guess_constant, 0)
        }
    }
    guesses_control$guess_constant <- guess_constant
    guesses_control$guesses <- guess_second
    e <- RunCoxRegression_Guesses(df, time1, time2, event0, names, Term_n, tform, keep_constant,
         a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col)
    #
    return(e)
}

#' Calculates Schoenfeld residuals for a Cox Proportional Hazards regression
#' \code{RunCox_Schoenfeld_Residual} uses user provided data, time/event columns, vectors specifying the model, and options to calculate the residuals
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#'
#' @return returns a list of the final results
#' @export
#'
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(0.1, 0.1, 0.1, 0.1)
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCox_Schoenfeld_Residual(df, time1, time2, event, names, Term_n, tform,
#'                                 keep_constant, a_n, modelform, fir, der_iden, control)
#' @importFrom rlang .data

RunCox_Schoenfeld_Residual <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    setkeyv(df, c(time2, event0))
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
            ", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    control <- Def_Control(control)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_schoenfeld_transition(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
    return (e)
}


#' Performs basic Cox Proportional Hazards regression with competing risks
#' \code{RunCoxRegression_CR} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence, starting positions, and competing risks
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param cens_weight list of weights for censoring rate
#'
#' @return returns a list of the final results
#' @export
#'
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   2,   1,   2,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(0.1, 0.1, 0.1, 0.1)
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' #weights the probability that a row would continue to extend without censoring,
#' #    for risk group calculation 
#' cens_weight <- c(0.83, 0.37, 0.26, 0.34, 0.55, 0.23, 0.27)
#' #censoring weight is generated by the survival library finegray function, or by hand.
#' #The ratio of weight at event end point to weight at row endpoint is used.
#' e <- RunCoxRegression_CR(df, time1, time2, event, names, Term_n, tform,
#'      keep_constant, a_n, modelform, fir, der_iden, control, cens_weight)
#'
#' @importFrom rlang .data

RunCoxRegression_CR <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,cens_weight){
    setkeyv(df, c(time2, event0))
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        print("no events")
        stop()
    }
    if (control$verbose){
        print(paste(length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(time1,time2,event0)
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    if (length(a_n)<length(names)){
        print(paste("Parameters used: ",length(a_n),", Covariates used: ",length(names),
            ", Remaining filled with 0.01",sep=""))
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
        keep_constant <- keep_constant[seq_len(length(names))]
    }
    #
    a_n0 <- copy(a_n)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    e <- cox_ph_transition_CR(Term_n,tform,a_n,dfc,x_all, fir,der_iden, modelform, control,
        as.matrix(df[,ce, with = FALSE]),tu,cens_weight,keep_constant,term_tot)
    a_n <- a_n0
    return (e)
}

#' Performs basic Cox Proportional Hazards regression, Allows for multiple starting guesses on c++ side
#' \code{RunCoxRegression_Guesses_CPP} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions. Has additional options to starting with several initial guesses
#'
#' @param df data used for regression
#' @param time1 column used for time period starts
#' @param time2 column used for time period end
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param der_iden number for the subterm to test derivative at, only used for testing runs with a single varying parameter
#' @param control list of parameters controlling the convergence
#' @param guesses_control list of parameters to control how the guessing works
#' @param Strat_Col column to stratify by if needed
#'
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
#'                       "a"=c(0,   1,   1,   0,   1,   0,   1),
#'                       "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
#'                       "c"=c(10,  11,  10,  11,  12,  9,   11),
#'                       "d"=c(0,   0,   0,   1,   1,   1,   1),
#'                       "e"=c(0,   0,   1,   0,   0,   0,   1))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' a_n <- c(1.1, -0.1, 0.2, 0.5) #used to test at a specific point
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' guesses_control=list("maxiter"=10,"guesses"=10,"lin_min"=0.001,"lin_max"=1,
#'     "loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform","loglin_method"="uniform",strata=FALSE)
#' Strat_Col='e'
#' 
#' e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n,
#'                               tform, keep_constant, a_n, modelform, fir,
#'                               der_iden, control,guesses_control,Strat_Col)
#' @importFrom rlang .data
RunCoxRegression_Guesses_CPP <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col){
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n[[1]])
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    a_n_default <- rep(0,length(a_n[[1]]))
    for (i in seq_along(a_n[[1]])){
        a_n_default[i] <- a_n[[1]][i]
    }
    if (guesses_control$strata==FALSE){
        setkeyv(df, c(time2, event0))
        dfend <- df[get(event0)==1, ]
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
        x_all <- as.matrix(df[,all_names, with = FALSE])
        ce <- c(time1,time2,event0)
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        a_ns <- c(NaN)
        maxiters <- c(NaN)
        #
        for (i in 1:length(a_n)){
            a_n0 <- a_n[[i]]
            if (length(a_n0)<length(names)){
                print(paste("Parameters used: ",length(a_n0),", Covariates used: ",length(names),
                    ", Remaining filled with 0.01",sep=""))
                a_n0 <- c(a_n0, rep(0.01,length(names)-length(a_n0)))
            } else if (length(a_n0)>length(names)){
                print(paste("Parameters used: ",length(a_n0),", Covariates used: ",length(names),sep=""))
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
                keep_constant <- keep_constant[seq_len(length(names))]
            }
            #
            rmin <- guesses_control$rmin
            rmax <- guesses_control$rmax
            if (length(rmin)!=length(rmax)){
                if (control$verbose){
                    print("rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
                }
            }
            #
            keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir, modelform, control,keep_constant,term_tot)
            if (keep){
                if (is.nan(maxiters[1])){
                    a_ns <- c(a_n0)
                    maxiters <- c(guesses_control$maxiter)
                } else {
                    a_ns <- c(a_ns, a_n0)
                    maxiters <- c(maxiters, guesses_control$maxiter)
                }
            }
        }
        while (length(maxiters) - length(a_n) < guesses_control$guesses){
            if (guesses_control$verbose){
                print(paste(length(maxiters)," valid guesses",sep=""))
            }
            if (length(rmin)!=length(rmax)){
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        if (grepl("_int",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$intercept_min,max=guesses_control$intercept_max) + a_n_default[i]
                        } else if (grepl("lin_exp_exp_slope",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$exp_slope_min,max=guesses_control$exp_slope_max) + a_n_default[i]
                        } else if (grepl("_slope",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else if (grepl("loglin",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$exp_min,max=guesses_control$exp_max) + a_n_default[i]
                        } else if (grepl("lin",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print(paste("tform not implemented ", tform[i],sep=""))
                            stop()
                        }
                    } else {
                        a_n0[i] <- a_n_default[i]
                    }
                }
            } else {
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        a_n0[i] <- runif(1,min=guesses_control$rmin[i],max=guesses_control$rmax[i]) + a_n_default[i]
                    } else {
                        a_n0[i] <- a_n_default[i]
                    }
                }
            }
            keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir, modelform, control,keep_constant,term_tot)
            if (keep){
                if (is.nan(maxiters[1])){
                    a_ns <- c(a_n0)
                    maxiters <- c(guesses_control$maxiter)
                } else {
                    a_ns <- c(a_ns,a_n0)
                    maxiters <- c(maxiters,guesses_control$maxiter)
                }
            }
        }
        #
        control$maxiters <- c(maxiters,control$maxiter)
        control$guesses <- length(maxiters)-1
        #
        e <- cox_ph_transition_guess(Term_n,tform,matrix(a_ns,nrow=length(control$maxiters)-1,byrow=T),dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot)
        #fine
        return (e)
    } else {
        #
        dfend <- df[get(event0)==1, ]
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        #
        for (i in seq_along(uniq)){
            df0 <- dfend[get(Strat_Col)==uniq[i],]
            tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
            if (length(tu0)==0){
                if (control$verbose){
                    print(paste("no events for strata group:",uniq[i],sep=" "))
                }
                df <- df[get(Strat_Col)!=uniq[i],]
            }
        }
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        if (control$verbose){
            print(paste(length(uniq)," strata used",sep=" "))
        }
        #
        setkeyv(df, c(time2, event0, Strat_Col))
        dfend <- df[get(event0)==1, ]
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
        x_all <- as.matrix(df[,all_names, with = FALSE])
        ce <- c(time1,time2,event0,Strat_Col)
        #
        #
        t_check <- Check_Trunc(df,ce)
        df <- t_check$df
        ce <- t_check$ce
        #
        a_ns <- c(NaN)
        maxiters <- c(NaN)
        #
        for (i in 1:length(a_n)){
            a_n0 <- a_n[[i]]
            if (length(a_n0)<length(names)){
                print(paste("Parameters used: ",length(a_n0),", Covariates used: ",length(names),
                    ", Remaining filled with 0.01",sep=""))
                a_n0 <- c(a_n0, rep(0.01,length(names)-length(a_n0)))
            } else if (length(a_n0)>length(names)){
                print(paste("Parameters used: ",length(a_n0),", Covariates used: ",length(names),sep=""))
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
                keep_constant <- keep_constant[seq_len(length(names))]
            }
            #
            rmin <- guesses_control$rmin
            rmax <- guesses_control$rmax
            if (length(rmin)!=length(rmax)){
                if (control$verbose){
                    print("rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
                }
            }
            #
            keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir, modelform, control,keep_constant,term_tot)
            if (keep){
                if (is.nan(maxiters[1])){
                    a_ns <- c(a_n0)
                    maxiters <- c(guesses_control$maxiter)
                } else {
                    a_ns <- c(a_ns, a_n0)
                    maxiters <- c(maxiters, guesses_control$maxiter)
                }
            }
        }
        while (length(maxiters) - length(a_n) < guesses_control$guesses){
            if (guesses_control$verbose){
                print(paste(length(maxiters)," valid guesses",sep=""))
            }
            if (length(rmin)!=length(rmax)){
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        if (grepl("_int",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$intercept_min,max=guesses_control$intercept_max) + a_n_default[i]
                        } else if (grepl("lin_exp_exp_slope",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$exp_slope_min,max=guesses_control$exp_slope_max) + a_n_default[i]
                        } else if (grepl("_slope",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else if (grepl("loglin",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$exp_min,max=guesses_control$exp_max) + a_n_default[i]
                        } else if (grepl("lin",tform[i],fixed=FALSE)){
                            a_n0[i] <- runif(1,min=guesses_control$lin_min,max=guesses_control$lin_max) + a_n_default[i]
                        } else {
                            print(paste("tform not implemented ", tform[i],sep=""))
                            stop()
                        }
                    } else {
                        a_n0[i] <- a_n_default[i]
                    }
                }
            } else {
                for (i in seq_along(tform)){
                    if (guesses_control$guess_constant[i]==0){
                        a_n0[i] <- runif(1,min=guesses_control$rmin[i],max=guesses_control$rmax[i]) + a_n_default[i]
                    } else {
                        a_n0[i] <- a_n_default[i]
                    }
                }
            }
            keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir, modelform, control,keep_constant,term_tot)
            if (keep){
                if (is.nan(maxiters[1])){
                    a_ns <- c(a_n0)
                    maxiters <- c(guesses_control$maxiter)
                } else {
                    a_ns <- c(a_ns,a_n0)
                    maxiters <- c(maxiters,guesses_control$maxiter)
                }
            }
        }
        #
        control$maxiters <- c(maxiters,control$maxiter)
        control$guesses <- length(maxiters)-1
        #
        e <- cox_ph_transition_guess_strata(Term_n,tform,matrix(a_ns,nrow=length(control$maxiters)-1,byrow=T),dfc,x_all, fir,der_iden, modelform, control,
            as.matrix(df[,ce, with = FALSE]),tu,keep_constant,term_tot, uniq)
        #fine
        return (e)
    }
    return (NULL)
}
