#' Performs basic Poisson regression using the omnibus function
#' \code{RunPoissonRegression_Omnibus} uses user provided data, time/event columns,
#'  vectors specifying the model, and options to control the convergence and starting positions.
#'  Has additional options to starting with several initial guesses
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' @param Strat_Col column to stratify by if needed
#' @param model_control controls which alternative model options are used
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
#' pyr <- "Ending_Age"
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
#' e <- RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant,
#'                                   a_n, modelform, fir, der_iden, control,Strat_Col)
#' @importFrom rlang .data
RunPoissonRegression_Omnibus <- function(df, pyr="pyr", event0="event", names=c("CONST"), Term_n=c(0), tform="loglin", keep_constant=c(0), a_n=c(0), modelform="M", fir=0, der_iden=0, control=list(),Strat_Col="null",model_control=list()){
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    control <- Def_Control(control)
    model_control <- Def_model_control(model_control)
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
    if (sum(df[,event0, with = FALSE])==0){
        if (control$verbose){
            print("no events")
        }
        stop()
    }
    if ("CONST" %in% names){
        if ("CONST" %in% names(df)){
            #fine
        } else {
            df$CONST <- 1
        }
    }
    if (model_control$strata==TRUE){
        val <- factorize(df, Strat_Col)
        df0 <- val$df
        #
    } else {
        df0 <- data.table("a"=c(0,0))
        val <- list(cols=c("a"))
    }
    #
    setkeyv(df, c(pyr, event0))
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(pyr,event0)
    #
    a_ns <- c()
    for (i in a_n){
        a_ns <- c(a_ns, i)
    }
    #
    e <- pois_Omnibus_transition(as.matrix(df[,ce, with = FALSE]),Term_n,tform,
                                 matrix(a_ns,nrow=length(control$maxiters)-1,byrow=TRUE),
                                 dfc,x_all, fir,der_iden, modelform, control,keep_constant,
                                 term_tot,as.matrix(df0[,val$cols, with=FALSE]),model_control)
    #fine
    return (e)
}





#' Performs basic poisson regression
#' \code{RunPoissonRegression} uses user provided data, person-year/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' df$pyr <- df$Ending_Age - df$Starting_Age
#' pyr <- 'pyr'
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
#'    'dose_abs_max'=100.0,'verbose'=FALSE, 'double_step'=1)
#' 
#' e <- RunPoissonRegression(df, pyr, event, names, Term_n, tform, keep_constant,
#'                           a_n, modelform, fir, der_iden, control)
#' @export
#'
RunPoissonRegression <- function(df, pyr, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control){
    control <- Def_Control(control)
    control$maxiters <- c(1,control$maxiter)
    control$guesses <- 1
    e <- RunPoissonRegression_Omnibus(df, pyr, event0, names, Term_n, tform, keep_constant,
                                      a_n, modelform, fir, der_iden, control)
    #
    return (e)
}

#' Performs poisson regression with no derivative calculations
#' \code{RunPoissonRegression_Single} uses user provided data, person-year/event columns, vectors specifying the model, and returns the results
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#'
#' @return returns a list of the final results
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
#' df$pyr <- df$Ending_Age - df$Starting_Age
#' pyr <- 'pyr'
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(0.1, 0.1, 0.1, 0.1)
#' 
#' keep_constant <- c(0,0,0,0)
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,'dbeta_max' = 0.5,
#'             'deriv_epsilon' = 1e-3, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'             'verbose'=FALSE, 'double_step'=1)
#' 
#' e <- RunPoissonRegression_Single(df, pyr, event, names, Term_n, tform, a_n, modelform, fir, control)
#' @export
#'
RunPoissonRegression_Single <- function(df, pyr, event0, names, Term_n, tform, a_n, modelform, fir, control,keep_constant=rep(0,length(names))){
    control <- Def_Control(control)
    control$maxiters <- c(1,control$maxiter)
    control$guesses <- 1
    e <- RunPoissonRegression_Omnibus(df, pyr, event0, names, Term_n, tform, keep_constant,
                                      a_n, modelform, fir, 0, control, model_control=list("single"=TRUE))
    return (e)
}

#' Performs poisson regression with strata effect
#' \code{RunPoissonRegression_STRATA} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' @param Strat_Cols column to stratify by
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
#' df$pyr <- df$Ending_Age - df$Starting_Age
#' pyr <- 'pyr'
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
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,'dbeta_max' = 0.5,
#'              'deriv_epsilon' = 1e-3,'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'              'verbose'=FALSE, 'double_step'=1)
#' Strat_Col <- c("e")
#' e <- RunPoissonRegression_STRATA(df, pyr, event, names, Term_n, tform, keep_constant,
#'      a_n, modelform, fir, der_iden, control, Strat_Col)
#'
RunPoissonRegression_STRATA <- function(df, pyr, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, Strat_Cols){
    control <- Def_Control(control)
    control$maxiters <- c(1,control$maxiter)
    control$guesses <- 1
    e <- RunPoissonRegression_Omnibus(df, pyr, event0, names, Term_n, tform, keep_constant, a_n,
                                      modelform, fir, der_iden, control,Strat_Cols, model_control=list("strata"=TRUE))
    return (e)
}

#' Performs poisson regression with strata effect and no derivative calculation
#' \code{RunPoissonRegression_STRATA_Single} uses user provided data, time/event columns, and vectors specifying the model
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
#' @param event0 column used for event status
#' @param names columns names for elements of the model, used to identify data columns
#' @param Term_n term numbers for each element of the model
#' @param tform subterm type for each element of the model
#' @param keep_constant vector of 0/1 to identify parameters to force to be constant
#' @param a_n starting parameters for regression
#' @param modelform string specifying the model type
#' @param fir term number for the initial term, used for models of the form T0*f(Ti) in which the order matters
#' @param control list of parameters controlling the convergence
#' @param Strat_Cols column to stratify by
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
#' df$pyr <- df$Ending_Age - df$Starting_Age
#' pyr <- 'pyr'
#' event <- "Cancer_Status"
#' names <- c('a','b','c','d')
#' Term_n <- c(0,1,1,2)
#' tform <- c("loglin","lin","lin","plin")
#' modelform <- "M"
#' fir <- 0
#' a_n <- c(0.1, 0.1, 0.1, 0.1)
#' 
#' keep_constant <- c(0,0,0,0)
#' 
#' control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 5,'halfmax' = 5,'epsilon' = 1e-3,'dbeta_max' = 0.5,
#'              'deriv_epsilon' = 1e-3,'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'              'verbose'=FALSE, 'double_step'=1)
#' Strat_Col <- c("e")
#' e <- RunPoissonRegression_STRATA_Single(df, pyr, event, names, Term_n, tform, keep_constant,
#'      a_n, modelform, fir, control, Strat_Col)
#'
RunPoissonRegression_STRATA_Single <- function(df, pyr, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, Strat_Cols){
    control <- Def_Control(control)
    control$maxiters <- c(1,control$maxiter)
    control$guesses <- 1
    e <- RunPoissonRegression_Omnibus(df, pyr, event0, names, Term_n, tform, keep_constant, a_n,
                                      modelform, fir, 0, control,Strat_Cols, model_control=list("strata"=TRUE,'single'=TRUE))
    return (e)
}

#' Performs basic poisson regression, with multiple guesses, starts with a single term
#' \code{RunPoissonRegression_Tier_Guesses} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions, with additional guesses
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' @param Strat_Cols column to stratify by
#' @param guesses_control list of parameters to control how the guessing works
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
#' df$pyr <- df$Ending_Age - df$Starting_Age
#' pyr <- 'pyr'
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
#'    'dose_abs_max'=100.0,'verbose'=FALSE,'double_step'=1)
#' guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001,"lin_max"=1,
#'   "loglin_min"=-1,"loglin_max"=1,"lin_method"="uniform",
#'   "loglin_method"="uniform",strata=TRUE,term_initial = c(0,1))
#' Strat_Cols=c('e')
#' 
#' e <- RunPoissonRegression_Tier_Guesses(df, pyr, event, names,
#'      Term_n, tform, keep_constant, a_n, modelform,
#'      fir, der_iden, control, guesses_control, Strat_Cols)
#'
#' @importFrom rlang .data
RunPoissonRegression_Tier_Guesses <- function(df, pyr, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, Strat_Cols){
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n)
    t_initial <- guesses_control$term_initial
    if (min(keep_constant)>0){
        print("Atleast one parameter must be free")
        stop()
    }
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
    name_initial <- c()
    term_n_initial <- c()
    tform_initial <- c()
    constant_initial <- c()
    a_n_initial <- c()
    guess_constant <- c()
    #
    for (i in seq_len(length(a_n))){
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
    e <- RunPoissonRegression_Guesses_CPP(df, pyr, event0, name_initial, term_n_initial, tform_initial,
         constant_initial, a_n_initial, modelform, fir, der_iden, control, guesses_control, Strat_Cols)
    if (guesses_control$verbose){
        print("INITIAL TERM COMPLETE")
        print(e)
    }
    #
    a_n_initial <- unlist(e$beta_0,use.names=FALSE)
    guess_constant <- c()
    j <- 1
    for (i in seq_len(length(a_n))){
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
    e <- RunPoissonRegression_Guesses_CPP(df, pyr, event0, names, Term_n, tform,
         keep_constant, a_n, modelform, fir, der_iden,
         control, guesses_control, Strat_Cols)
    #
    return(e)
}

#' Performs basic Poisson regression, Allows for multiple starting guesses on c++ side
#' \code{RunPoissonRegression_Guesses_CPP} uses user provided data, time/event columns, vectors specifying the model, and options to control the convergence and starting positions. Has additional options to starting with several initial guesses
#'
#' @param df data used for regression
#' @param pyr column used for person-years per row
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
#' @param model_control controls which alternative model options are used
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
#' pyr <- "Ending_Age"
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
#' e <- RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n,
#'                               tform, keep_constant, a_n, modelform, fir,
#'                               der_iden, control,guesses_control,Strat_Col)
#' @importFrom rlang .data
RunPoissonRegression_Guesses_CPP <- function(df, pyr, event0, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control,Strat_Col=c("null"),model_control=list()){
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    control <- Def_Control(control)
    if ("strata" %in% names(guesses_control)){
        if ("strata" %in% names(model_control)){
            if (guesses_control$strata != model_control$strata){
                if (guesses_control$verbose){
                    print("guesses_control and model_control have different strata options")
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
    #
    setkeyv(df, c(pyr, event0))
    all_names <- unique(names)
    #
    df <- Replace_Missing(df,all_names,0.0,control$verbose)
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    ce <- c(pyr,event0)
    #
    #
    a_ns <- c(NaN)
    maxiters <- c(NaN)
    #
    for (i in seq_len(length(a_n))){
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
    #
    a_n_mat <- matrix(a_ns,nrow=length(control$maxiters)-1,byrow=TRUE)
    a_n <- lapply(seq_len(nrow(a_n_mat)), function(i) a_n_mat[i,])
    e <- RunPoissonRegression_Omnibus(df, pyr, event0, names, Term_n, tform, keep_constant, a_n,
                                      modelform, fir, der_iden, control,model_control=model_control, Strat_Col=Strat_Col)
    #fine
    return (e)
}
    
