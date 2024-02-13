#' Performs Cox Proportional Hazards regression using the omnibus function
#'
#' \code{RunCoxRegression_Omnibus} uses user provided data, time/event columns,
#'       vectors specifying the model, and options to control the convergence
#'       and starting positions. Has additional options for starting with several
#'       initial guesses, using stratification, multiplicative loglinear 1-term,
#'       competing risks, and calculation without derivatives
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' a_n <- list(c(1.1, -0.1, 0.2, 0.5),c(1.6, -0.12, 0.3, 0.4))
#' #used to test at a specific point
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' 
#' keep_constant <- c(0,0,0,0)
#' der_iden <- 0
#' 
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiters' = c(5,5,5),
#'    'halfmax' = 5,'epsilon' = 1e-3, 'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE, 'dose_abs_max'=100.0,'verbose'=FALSE,
#'    'ties'='breslow','double_step'=1, "guesses"=2)
#' 
#' e <- RunCoxRegression_Omnibus(df, time1, time2, event,
#'                               names, Term_n, tform, keep_constant,
#'                               a_n, modelform, fir, der_iden, control,
#'                               model_control=list("single"=FALSE,
#'                               "basic"=FALSE, "CR"=FALSE, 'null'=FALSE))
#' @importFrom rlang .data
RunCoxRegression_Omnibus <- function(df, time1="start", time2="end", event0="event", names=c("CONST"), Term_n=c(0), tform="loglin", keep_constant=c(0), a_n=c(0), modelform="M", fir=0, der_iden=0, control=list(),Strat_Col="null", cens_weight=c(1), model_control=list(),Cons_Mat=as.matrix(c(0)),Cons_Vec=c(0)){
    control <- Def_Control(control)
    val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n,
                                 names, der_iden, Cons_Mat, Cons_Vec,control$verbose)
    Term_n <- val$Term_n
    tform <- val$tform
    keep_constant <- val$keep_constant
    a_n <- val$a_n
    der_iden <- val$der_iden
    names <- val$names
    Cons_Mat <- as.matrix(val$Cons_Mat)
    Cons_Vec <- val$Cons_Vec
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    if (control$verbose){
        if (any(val$Permutation != seq_along(tform))){
            message("Warning: model covariate order changed")
        }
    }
    model_control <- Def_model_control(model_control)
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    if (model_control$strata==FALSE){
        data.table::setkeyv(df, c(time2, event0))
        uniq <- c(0)
        ce <- c(time1,time2,event0)
    } else {
        #
        dfend <- df[get(event0)==1, ]
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]),
                            use.names=FALSE))
        #
        for (i in seq_along(uniq)){
            df0 <- dfend[get(Strat_Col)==uniq[i],]
            tu0 <- unlist(unique(df0[,time2,with=FALSE]), use.names=FALSE)
            if (length(tu0)==0){
                if (control$verbose){
                    message(paste("Warning: no events for strata group:",
                                 uniq[i],sep=" "))
                }
                df <- df[get(Strat_Col)!=uniq[i],]
            }
        }
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]),
                            use.names=FALSE))
        if (control$verbose){
            message(paste("Note:",length(uniq)," strata used",sep=" "))
        }
        #
        data.table::setkeyv(df, c(time2, event0, Strat_Col))
        ce <- c(time1,time2,event0,Strat_Col)
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            message("Error: no events")
        }
        stop()
    }
    if (control$verbose){
        message(paste("Note: ",length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    a_ns <- c()
    for (i in a_n){
        a_ns <- c(a_ns, i)
    }
    #
    if ("maxiters" %in% names(control)){
    	if (length(control$maxiters) == length(a_n)+1){
    		#all good, it matches
    	} else {
    		if (control$verbose){
                message(paste("Note: Initial starts:",length(a_n),
                      ", Number of iterations provided:",length(control$maxiters),". Colossus requires one more iteration counts than number of guesses (for best guess)",sep=" "))
            }
            if (length(control$maxiters) < length(a_n)+1){
		        additional <- length(a_n)+1 - length(control$maxiters)
		        control$maxiters <- c(control$maxiters, rep(1, additional))
	        } else {
	        	additional <- length(a_n)+1
	        	control$maxiters <- control$maxiters[1:additional]
	        }
    	}
	    if ("guesses" %in% names(control)){
	        #both are in
	        if (control$guesses+1 == length(control$maxiters)){
	            #all good, it matches
	        } else {
	            if (control$verbose){
                    message(paste("Error: guesses:",control["guesses"],
                          ", iterations per guess:",control["maxiters"],sep=" "))
                }
                stop()
	        }
	    } else {
	        control$guesses = length(control$maxiters)-1
	    }
	} else {
	    if ("guesses" %in% names(control)){
	    	if (control$guesses == length(a_n)){
	    		#both match, all good
    		} else {
    			control$guesses = length(a_n)
    		}
            control$maxiters = rep(1,control$guesses+1)
        } else {
            control$guesses = length(a_n)
            control$maxiters = c(rep(1,length(a_n)),control$maxiter)
        }
    }
    if (model_control$null){
        a_ns <- matrix(a_ns)
    } else {
        a_ns <- matrix(a_ns,nrow=length(control$maxiters)-1,byrow=TRUE)
    }
    #
    e <- cox_ph_Omnibus_transition(Term_n,tform,a_ns,dfc,x_all, fir,der_iden,
         modelform, control, as.matrix(df[,ce, with = FALSE]),tu,
         keep_constant,term_tot, uniq, cens_weight, model_control,
         Cons_Mat, Cons_Vec)
	if (is.nan(e$LogLik)){
		if (control$verbose){message("Invalid risk")}
		stop()
	}
    e$Parameter_Lists$names <- names
    return (e)
}

#' Performs basic Cox Proportional Hazards regression without special options
#'
#' \code{RunCoxRegression} uses user provided data, time/event columns,
#' vectors specifying the model, and options to control the convergence
#' and starting position
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform,
#'                      keep_constant, a_n, modelform, fir, der_iden, control)
#' @importFrom rlang .data

RunCoxRegression <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    #
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n,
                                  tform, keep_constant, a_n, modelform,
                                  fir, der_iden, control, model_control=list())
    #
    return (e)
}

#' Approximates how many events were due to baseline vs excess risk
#'
#' \code{RunCoxEventAssignment} uses user provided data, time/event columns,
#' vectors specifying the model, and options to control the convergence
#' and starting position, calculates approximated background and excess events
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform,
#'                      keep_constant, a_n, modelform, fir, der_iden, control)
#' @importFrom rlang .data

RunCoxEventAssignment <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    #
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1

    control <- Def_Control(control)
    val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n,
                                 names, der_iden, as.matrix(c(0)), c(0),control$verbose)
    Term_n <- val$Term_n
    tform <- val$tform
    keep_constant <- val$keep_constant
    a_n <- val$a_n
    der_iden <- val$der_iden
    names <- val$names
    model_control <- Def_model_control(list())
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    data.table::setkeyv(df, c(time2, event0))
    ce <- c(time1,time2,event0)
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            message("Error: no events")
        }
        stop()
    }
    if (control$verbose){
        message(paste("Note: ",length(tu)," risk groups",sep=""))
    }
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    #
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    e <- Assigned_Event_transition(matrix(c(0)),Term_n, tform, a_n, dfc, x_all, fir, der_iden, modelform, control, as.matrix(df[,ce, with = FALSE]), tu, keep_constant, term_tot, model_control)
    #
    return (e)
}

#' Performs basic Cox Proportional Hazards calculation with no derivative
#'
#' \code{RunCoxRegression_Single} uses user provided data, time/event columns, vectors specifying the model, and options and returns the log-likelihood
#'
#' @inheritParams R_template
#' @family Cox Wrapper Functions
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression_Single(df, time1, time2, event, names, Term_n, tform,
#'                              keep_constant, a_n, modelform, fir, control)
#'
#' @importFrom rlang .data

RunCoxRegression_Single <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, control){
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n,
         tform, keep_constant, a_n, modelform, fir, 0, control,
         model_control=list('single'=TRUE))
    #
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with a multiplicative log-linear model
#'
#' \code{RunCoxRegression_Basic} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @inheritParams R_template
#' @family Cox Wrapper Functions
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,
#'    'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE,
#'    'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxRegression_Basic(df, time1, time2, event, names, keep_constant,
#'                             a_n, der_iden, control)
#'
#' @importFrom rlang .data

RunCoxRegression_Basic <- function(df, time1, time2, event0, names, keep_constant, a_n, der_iden, control){
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names,
         rep(0,length(names)), rep('loglin',length(names)), keep_constant, a_n,
         "M", 0, der_iden, control, model_control=list("basic"=TRUE))
    #
    return (e)
}


#' Performs basic Cox Proportional Hazards regression with strata effect
#'
#' \code{RunCoxRegression_STRATA} uses user provided data,
#' time/event columns, vectors specifying the model, and options to control
#' the convergence and starting positions
#'
#' @inheritParams R_template
#' @family Cox Wrapper Functions
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' Strat_Col <- 'e'
#' 
#' e <- RunCoxRegression_STRATA(df, time1, time2, event, names, Term_n,
#'                              tform, keep_constant, a_n, modelform,
#'                              fir, der_iden, control,Strat_Col)
#'
RunCoxRegression_STRATA <- function(df, time1, time2, event0,  names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Col){
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n,
                                  tform, keep_constant, a_n, modelform, fir,
                                  der_iden, control,Strat_Col=Strat_Col,
                                  model_control=list("strata"=TRUE))
    return (e)
}


#' Calculates hazard ratios for a reference vector
#'
#' \code{RunCoxRegression} uses user provided data,  vectors specifying the model,
#' and options to calculate relative risk for every row in the provided data
#'
#' @inheritParams R_template
#' @family Plotting Wrapper Functions
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n, tform,
#'      keep_constant, a_n, modelform, fir, control)
#'
Cox_Relative_Risk <- function(df, time1, time2, event0,  names, Term_n, tform, keep_constant, a_n, modelform, fir, control, model_control=list()){
    control <- Def_Control(control)
    model_control <- Def_model_control(model_control)
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    #
    val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names)
    Term_n <- val$Term_n
    tform <- val$tform
    keep_constant <- val$keep_constant
    a_n <- val$a_n
    der_iden <- val$der_iden
    names <- val$names
    #
    all_names <- unique(names)
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    model_control$Risk_Subset <- TRUE
    e <- Plot_Omnibus_transition(Term_n, tform, a_n, dfc, x_all, fir,
                                 0, modelform, control, matrix(c(0)),
                                 c(1), keep_constant, term_tot, c(0),
                                 c(0), model_control)
    return (e)
}

#' Performs basic Cox Proportional Hazards regression with the null model
#'
#' \code{RunCoxRegression} uses user provided data and time/event columns
#' to calculate the log-likelihood with constant hazard ratio
#'
#' @inheritParams R_template
#' @family Cox Wrapper Functions
#' @return returns a list of the final results
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' event <- "Cancer_Status"
#' 
#' control <- list("Ncores"=2,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' 
#' e <- RunCoxNull(df, time1, time2, event, control)
#'
RunCoxNull <- function(df, time1, time2, event0,control){
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0,  control=control,
                                  model_control=list("null"=TRUE))
    #
    return (e)

}


#' Performs Cox Proportional Hazard model plots
#'
#' \code{RunCoxPlots} uses user provided data, time/event columns,
#' vectors specifying the model, and options to choose and save plots
#'
#' @inheritParams R_template
#'
#' @return saves the plots in the current directory and returns a string
#' @export
#' @family Plotting Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' #setting maxiter below 0 forces the function to calculate the score
#' # and return
#' plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")), "studyID"="UserID",
#'                   'verbose'=FALSE)
#' 
#' RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant,
#'             a_n, modelform, fir, control, plot_options)
#'
RunCoxPlots <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options, model_control=list()){
    control <- Def_Control(control)
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
        stop()
    }
    if (plot_options$verbose){
        message("Note: Starting Plot Function")
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    data.table::setkeyv(df, c(time2, event0))
    base  <- NULL
    der_iden <- 0
    Plot_Type <- plot_options$type
    if (plot_options$verbose){
        message("Note: Getting Plot Info")
    }
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        message("Error: no events")
        stop()
    }
    if (plot_options$verbose){
        message(paste("Note: ",length(tu)," risk groups",sep=""))
    }
    if ("type" %in% names(plot_options)){
        #fine
    } else {
        if (plot_options$verbose){
            message("Error: Plot type not given")
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
                        message("Error: Stratification Column not in dataframe")
                    }
                    stop()
                }
            } else {
                if (plot_options$verbose){
                    message("Error: Stratification Column not given")
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
                    dose_col <- unlist(plot_options$cov_cols,
                                       use.names=FALSE)[cov_i]
                    if (dose_col%in% names(df)){
                        #fine
                    } else {
                        if (plot_options$verbose){
                            message("Error: Covariate column "+
                                   dose_col+" is not in the dataframe")
                        }
                        stop()
                    }
                }
            } else {
                if (plot_options$verbose){
                    message("Error: dose column not given")
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
                        message("Error: ID column is not in the dataframe")
                    }
                    stop()
                }
            } else {
                if (plot_options$verbose){
                    message("Error: ID column not given")
                }
                stop()
            }
        }
    }
    model_control <- Def_model_control(model_control)
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    if (tolower(Plot_Type[1])=="surv"){
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
    verbose <- data.table::copy(plot_options$verbose)
    verbosec <- data.table::copy(control$verbose)
    maxiterc <- data.table::copy(control$maxiter)
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (length(tu)==0){
        if (control$verbose){
            message("Error: no events")
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
    #
    #
    control$maxiters <- c(-1, -1)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n,
                                  tform, keep_constant, a_n, modelform, fir,
                                  der_iden, control, model_control)
    control$maxiter <- maxiterc
    b <- e$beta_0
    er <- e$Standard_Deviation
    #
    #
    if (tolower(Plot_Type[1])=="surv"){
        if (verbose){
            message("Note: starting ph_plot")
        }
        #
        if (plot_options$strat_haz==FALSE){
        	if (verbose){
		        message("Note: nonStratified survival curve calculation")
		    }
		    model_control$Surv <- TRUE
		    e <- Plot_Omnibus_transition(Term_n, tform, a_n, dfc, x_all, fir,
		                                 der_iden, modelform, control,
		                                 as.matrix(df[,ce, with = FALSE]),tu,
		                                 keep_constant, term_tot, c(0), c(0),
		                                 model_control)
		    #
		    t <- c()
		    h <- c()
		    ch <- c()
		    surv <- c()
		    if (verbose){
		        message("Note: writing survival data")
		    }
		    dft <- data.table::data.table("time"=tu,"base"=e$baseline,
		                      "basehaz"=e$standard_error)
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
		        CoxMartingale(verbose, df, time1, time2, event0, e, t, ch,
		                      plot_options$cov_cols,
		                      Plot_Type[2], age_unit,plot_options$studyID)
		        #
		    }
		    if (plot_options$surv_curv==TRUE){
		        CoxSurvival(t,h,ch,surv,Plot_Type[2],verbose,
		                    plot_options$time_lims, age_unit)
            }
        } else {
        	age_unit <- plot_options$age_unit
        	if (verbose){
		        message("Note: Stratified survival curve calculation")
		    }
            if (plot_options$surv_curv==TRUE){
                model_control$strata <- TRUE
                CoxStratifiedSurvival(verbose, df, event0, time1, time2,
                     all_names,Term_n, tform, a_n, er, fir, der_iden,
                     modelform, control, keep_constant, Plot_Type,
                     plot_options$Strat_Col, plot_options$time_lims,age_unit)
            }
        }
        if (plot_options$KM==TRUE){
            #
            CoxKaplanMeier(verbose, verbosec, plot_options$studyID,
                           all_names,df,event0,time1, time2,tu,Term_n,
                           tform, a_n, er, fir, der_iden, modelform,
                           control,keep_constant, Plot_Type,age_unit)
        }
    } else if (tolower(Plot_Type[1])=="risk"){
        CoxRisk(verbose, df, event0, time1, time2, names,Term_n, tform,
                a_n, fir, der_iden, modelform, control,keep_constant,
                Plot_Type, b, er)
    } else if (tolower(Plot_Type[1])=="schoenfeld"){
        age_unit <- plot_options$age_unit
        PlotCox_Schoenfeld_Residual(df, time1, time2, event0, names, Term_n,
                                    tform, keep_constant, a_n, modelform, fir,
                                    der_iden, control,age_unit,Plot_Type[2])
    }
    return ("Passed")
}

#' Performs basic cox regression, with multiple guesses, starts with
#' solving for a single term
#'
#' \code{RunCoxRegression_Tier_Guesses} uses user provided data, time/event
#' columns, vectors specifying the model, and options to control the
#' convergence and starting positions, with additional guesses
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,
#'    'epsilon' = 1e-3,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3,
#'    'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'    'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001,
#'    "lin_max"=1,"loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform",
#'    "loglin_method"="uniform",strata=TRUE,term_initial = c(0,1))
#' Strat_Col <- 'e'
#' 
#' e <- RunCoxRegression_Tier_Guesses(df, time1, time2, event, names,
#'                                    Term_n, tform, keep_constant,
#'                                    a_n, modelform, fir, der_iden,
#'                                    control,guesses_control,
#'                                    Strat_Col)
#'
#' @importFrom rlang .data
RunCoxRegression_Tier_Guesses <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col,model_control=list(),cens_weight=c(0)){
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n)
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
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
    rmin <- guesses_control$rmin
    rmax <- guesses_control$rmax
    if (length(rmin)!=length(rmax)){
        if (control$verbose){
            message("Warning: rmin/rmax not equal size, lin/loglin min/max used")
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
    e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event0, name_initial,
                                      term_n_initial, tform_initial,
                                      constant_initial, a_n_initial,
                                      modelform, fir, der_iden, control,
                                      guesses_control,Strat_Col)
    #
    if (guesses_control$verbose){
        message("Note: INITIAL TERM COMPLETE")
        message(e)
    }
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
    e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event0, names,
         Term_n, tform,keep_constant, a_n, modelform, fir, der_iden, control,
         guesses_control,Strat_Col)
    #
    return(e)
}


#' Performs basic Cox Proportional Hazards regression with competing risks
#'
#' \code{RunCoxRegression_CR} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence, starting positions, and censoring adjustment
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
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
    control <- Def_Control(control)
    control$maxiters <- c(1, control$maxiter)
    control$guesses <- 1
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n, tform, keep_constant,
                                  a_n, modelform, fir, der_iden, control,cens_weight=cens_weight,
                                  model_control=list("CR"=TRUE))
    #
    return (e)
}

#' Performs basic Cox Proportional Hazards regression, Allows for multiple starting guesses on c++ side
#'
#' \code{RunCoxRegression_Guesses_CPP} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions. Has additional options to starting with several initial guesses
#'
#' @inheritParams R_template
#'
#' @return returns a list of the final results
#' @export
#' @family Cox Wrapper Functions
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
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
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,
#'    'dbeta_max' = 0.5,'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' guesses_control=list("maxiter"=10,"guesses"=10,"lin_min"=0.001,"lin_max"=1,
#'     "loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform","loglin_method"="uniform",strata=FALSE)
#' Strat_Col <- 'e'
#' 
#' e <- RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n,
#'                               tform, keep_constant, a_n, modelform, fir,
#'                               der_iden, control,guesses_control,Strat_Col)
#' @importFrom rlang .data
RunCoxRegression_Guesses_CPP <- function(df, time1, time2, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col,model_control=list(),cens_weight=c(0)){
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    control <- Def_Control(control)
    if ("strata" %in% names(guesses_control)){
        if ("strata" %in% names(model_control)){
            if (guesses_control$strata != model_control$strata){
                if (guesses_control$verbose){
                    message("Error: guesses_control and model_control have different strata options")
                }
                stop()
            }
        } else {
            model_control$strata <- guesses_control$strata
        }
    } else if ("strata" %in% names(model_control)){
        guesses_control$strata <- model_control$strata
    }
    guesses_control <- Def_Control_Guess(guesses_control, a_n[[1]])
    model_control <- Def_model_control(model_control)
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
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
        data.table::setkeyv(df, c(time2, event0))
        dfend <- df[get(event0)==1, ]
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
        if (length(tu)==0){
            if (guesses_control$verbose){
                message("Error: no events")
            }
            stop()
        }
        if (guesses_control$verbose){
            message(paste("Note: ",length(tu)," risk groups",sep=""))
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
                    message(paste("Warning: no events for strata group:",uniq[i],sep=" "))
                }
                df <- df[get(Strat_Col)!=uniq[i],]
            }
        }
        uniq <- sort(unlist(unique(df[,Strat_Col, with = FALSE]), use.names=FALSE))
        if (control$verbose){
            message(paste("Note:",length(uniq)," strata used",sep=" "))
        }
        #
        data.table::setkeyv(df, c(time2, event0, Strat_Col))
        dfend <- df[get(event0)==1, ]
        tu <- sort(unlist(unique(dfend[,time2, with = FALSE]),use.names=FALSE))
        if (length(tu)==0){
            if (guesses_control$verbose){
                message("Error: no events")
            }
            stop()
        }
        if (guesses_control$verbose){
            message(paste("Note: ",length(tu)," risk groups",sep=""))
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
    }
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    #
    dat_val <- Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant,
                                  a_n, x_all, a_n_default, modelform, fir, control,
                                  guesses_control, model_control)
    a_ns <- dat_val$a_ns
    maxiters <- dat_val$maxiters
    #
    control$maxiters <- c(maxiters,control$maxiter)
    control$guesses <- length(maxiters)
    #
    #
    #fine
    a_n_mat <- matrix(a_ns,nrow=length(control$maxiters)-1,byrow=TRUE)
    a_n <- lapply(seq_len(nrow(a_n_mat)), function(i) a_n_mat[i,])
    e <- RunCoxRegression_Omnibus(df, time1, time2, event0, names, Term_n, tform, keep_constant,
                                  a_n, modelform, fir, der_iden, control, Strat_Col=Strat_Col,
                                  model_control=model_control,cens_weight=cens_weight)
    return (e)
}
