#' Performs checks to gather a list of guesses and iterations
#' \code{Gather_Guesses_CPP} called from within R, uses a list of options and the model definition to generate a list of parameters and iterations that do not produce errors
#'
#' @inheritParams R_template
#' @param a_n_default center of parameter distribution guessing scope
#' @param dfc vector matching subterm number to matrix column
#' @param x_all covariate matrix
#'
#' @return returns a list of the final results
#' @importFrom rlang .data
Gather_Guesses_CPP <- function(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control, model_control=list()){
    if (typeof(a_n)!="list"){
        a_n <- list(a_n)
    }
    term_tot <- max(Term_n)+1
    a_ns <- c(NaN)
    maxiters <- c(NaN)
    model_control <- Def_model_control(model_control)
    val <- Def_modelform_fix(control,model_control,modelform,Term_n)
    modelform <- val$modelform
    model_control <- val$model_control
    #
    for (i in seq_len(length(a_n))){
        a_n0 <- a_n[[i]]
        if (length(a_n0)<length(names)){
            message(paste("Warning: Parameters used: ",length(a_n0),", Covariates used: ",length(names),
                ", Remaining filled with 0.01",sep=""))
            a_n0 <- c(a_n0, rep(0.01,length(names)-length(a_n0)))
        } else if (length(a_n0)>length(names)){
            message(paste("Error: Parameters used: ",length(a_n0),", Covariates used: ",length(names),
                        sep=""))
            stop()
        }
        if (length(Term_n)<length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        } else if (length(Term_n)>length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(tform)<length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
            stop()
        } else if (length(tform)>length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
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
                message("Warning: rmin and rmax lists not equal size, defaulting to lin and loglin min/max values")
            }
        }
        #
        keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir, modelform,
                                      control,model_control,keep_constant,term_tot)
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
            message(paste("Note: ",length(maxiters)," valid guesses",sep=""))
        }
        if (length(rmin)!=length(rmax)){
            for (i in seq_along(tform)){
                if (guesses_control$guess_constant[i]==0){
                    if (grepl("_int",tform[i],fixed=FALSE)){
                        a_n0[i] <- runif(1,min=guesses_control$intercept_min,
                                           max=guesses_control$intercept_max) + a_n_default[i]
                    } else if (grepl("lin_exp_exp_slope",tform[i],fixed=FALSE)){
                        a_n0[i] <- runif(1,min=guesses_control$exp_slope_min,
                                           max=guesses_control$exp_slope_max) + a_n_default[i]
                    } else if (grepl("_slope",tform[i],fixed=FALSE)){
                        a_n0[i] <- runif(1,min=guesses_control$lin_min,
                                           max=guesses_control$lin_max) + a_n_default[i]
                    } else if (grepl("loglin",tform[i],fixed=FALSE)){
                        a_n0[i] <- runif(1,min=guesses_control$exp_min,
                                           max=guesses_control$exp_max) + a_n_default[i]
                    } else if (grepl("lin",tform[i],fixed=FALSE)){
                        a_n0[i] <- runif(1,min=guesses_control$lin_min,
                                           max=guesses_control$lin_max) + a_n_default[i]
                    } else {
                        message(paste("Error: tform not implemented ", tform[i],sep=""))
                        stop()
                    }
                } else {
                    a_n0[i] <- a_n_default[i]
                }
            }
        } else {
            for (i in seq_along(tform)){
                if (guesses_control$guess_constant[i]==0){
                    a_n0[i] <- runif(1,min=guesses_control$rmin[i],
                                       max=guesses_control$rmax[i]) + a_n_default[i]
                } else {
                    a_n0[i] <- a_n_default[i]
                }
            }
        }
        keep <- risk_check_transition(Term_n,tform,a_n0,dfc,x_all, fir,
                                      modelform, control,model_control,keep_constant,term_tot)
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
    list("a_ns"=a_ns, "maxiters"=maxiters)
}

#' Corrects the order of terms/formula/etc
#' \code{Correct_Formula_Order} checks the order of formulas given and corrects any ordering issues
#'
#' @inheritParams R_template
#'
#' @return returns the corrected lists
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' Term_n <- c(0,1,1,0,0)
#' tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
#' keep_constant <- c(0,0,0,1,0)
#' a_n <- c(1,2,3,4,5)
#' names <- c("a","a","a","a","a")
#' val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names)
#' Term_n <- val$Term_n
#' tform <- val$tform
#' keep_constant <- val$keep_constant
#' a_n <- val$a_n
#' der_iden <- val$der_iden
#' names <- val$names
#'
Correct_Formula_Order <- function(Term_n, tform, keep_constant, a_n, names,der_iden=0,verbose=FALSE){
    #
    if (min(keep_constant)>0){
        message("Error: Atleast one parameter must be free")
        stop()
    }
    if (((typeof(a_n)=="list")&&(length(a_n)==1))||(typeof(a_n)!="list")){
        #
        if (typeof(a_n)=="list"){
            a_n <- a_n[[1]]
        }
        if (length(a_n)<length(names)){
            message(paste("Warning: Parameters used: ",length(a_n),", Covariates used: ",length(names),
                ", Remaining filled with 0.01",sep=""))
            a_n <- c(a_n, rep(0.01,length(names)-length(a_n)))
        } else if (length(a_n)>length(names)){
            message(paste("Error: Parameters used: ",length(a_n),", Covariates used: ",length(names),
                        sep=""))
            stop()
        }
        if (length(Term_n)<length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        } else if (length(Term_n)>length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(tform)<length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
            stop()
        } else if (length(tform)>length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
            stop()
        }
        if (length(keep_constant)<length(names)){
            keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
        } else if (length(keep_constant)>length(names)){
            keep_constant <- keep_constant[seq_len(length(names))]
        }
        #
        df <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant,
                         "a_n"=a_n, "names"=names, "iden_const"=rep(0,length(names)),"current_order"=1:(length(tform)))
        df$iden_const[[der_iden+1]] <- 1
        tform_order <- c("loglin", "lin", "plin", "loglin_slope", "loglin_top",
                         "lin_slope", "lin_int", "quad_slope", "step_slope",
                         "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope",
                         "lin_exp_int", "lin_exp_exp_slope")
        tform_iden <- match(tform,tform_order)
        df$tform_order <- tform_iden
        if (verbose){
            message(df)
        }
        keycol <-c("Term_n","names","tform_order")
        setorderv(df, keycol)
        a_n <- df$a_n
    } else {
        a_0 <- a_n[[1]]
        for (a_i in a_n){
            if (length(a_i)!=length(a_0)){
                message(paste("Error: Parameters used in first option: ",length(a_0),
                            ", Parameters used in different option: ",length(a_i),
                            ", please fix parameter length",sep=""))
                stop()
            }
        }
        #
        if (length(a_0)<length(names)){
            message(paste("Warning: Parameters used: ",length(a_0),", Covariates used: ",length(names),
                ", Remaining filled with 0.01",sep=""))
            for (i in seq_len(length(a_n))){
                a_n[[i]] <- c(a_n[[i]], rep(0.01,length(names)-length(a_n[[i]])))
            }
        } else if (length(a_0)>length(names)){
            message(paste("Error: Parameters used: ",length(a_0),", Covariates used: ",length(names),
                        sep=""))
            stop()
        }
        if (length(Term_n)<length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        } else if (length(Term_n)>length(names)){
            message(paste("Error: Terms used: ",length(Term_n),", Covariates used: ",length(names),sep=""))
            stop()
        }
        if (length(tform)<length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
            stop()
        } else if (length(tform)>length(names)){
            message(paste("Error: Term types used: ",length(tform),", Covariates used: ",length(names),
                        sep=""))
            stop()
        }
        if (length(keep_constant)<length(names)){
            keep_constant <- c(keep_constant, rep(0.01,length(names)-length(keep_constant)))
        } else if (length(keep_constant)>length(names)){
            keep_constant <- keep_constant[seq_len(length(names))]
        }
        #
        df <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant,
                         "names"=names,"iden_const"=rep(0,length(names)),"current_order"=1:(length(tform)))
        df$iden_const[[der_iden+1]] <- 1
        for (i in seq_len(length(a_n))){
            df[[paste("a_",i,sep="")]] <- a_n[[i]]
        }
        tform_order <- c("loglin", "lin", "plin", "loglin_slope", "loglin_top",
                         "lin_slope", "lin_int", "quad_slope", "step_slope",
                         "step_int", "lin_quad_slope", "lin_quad_int", "lin_exp_slope",
                         "lin_exp_int", "lin_exp_exp_slope")
        tform_iden <- match(tform,tform_order)
        df$tform_order <- tform_iden
        if (verbose){
            message(df)
        }
        keycol <-c("Term_n","names","tform_order")
        setorderv(df, keycol)
        for (i in seq_len(length(a_n))){
            a_n[[i]] <- df[[paste("a_",i,sep="")]]
        }
    }
    a <- df$tform
    for (i in seq_len(length(a))){
        if (i<length(a)){
            if ((a[i]=="loglin_slope")){
                if (a[i+1]!="loglin_top"){
                    message("Error: loglin_top missing")
                    stop()
                }
            } else if  ((a[i]=="lin_slope")){
                if (a[i+1]!="lin_int"){
                    message("Error: lin_int missing")
                    stop()
                }
            } else if ((a[i]=="step_slope")){
                if (a[i+1]!="step_int"){
                    message("Error: step_int missing")
                    stop()
                }
            } else if ((a[i]=="lin_quad_slope")){
                if (a[i+1]!="lin_quad_int"){
                    message("Error: lin_quad_int missing")
                    stop()
                }
            } else if ((a[i]=="lin_exp_slope")){
                if (a[i+1]!="lin_exp_int"){
                    message("Error: lin_exp_int missing")
                    stop()
                }
            } else if ((a[i]=="lin_exp_int")){
                if (a[i+1]!="lin_exp_exp_slope"){
                    message("Error: lin_exp_exp_slope missing")
                    stop()
                }
            }
        } 
        if (i>1){
            if  ((a[i]=="lin_int")){
                if (a[i-1]!="lin_slope"){
                    message("Error: lin_slope missing")
                    stop()
                }
            } else if ((a[i]=="step_int")){
                if (a[i-1]!="step_slope"){
                    message("Error: step_slope missing")
                    stop()
                }
            } else if ((a[i]=="lin_quad_int")){
                if (a[i-1]!="lin_quad_slope"){
                    message("Error: lin_quad_slope missing")
                    stop()
                }
            } else if ((a[i]=="lin_exp_int")){
                if (a[i-1]!="lin_exp_slope"){
                    message("Error: lin_exp_slope missing")
                    stop()
                }
            } else if ((a[i]=="lin_exp_exp_slope")){
                if (a[i-1]!="lin_exp_int"){
                    message("Error: lin_exp_int missing")
                    stop()
                }
            }
        }
    }
    a_temp <- df$iden_const
    der_iden <- which(a_temp==1)
    list("Term_n"=df$Term_n, "tform"=df$tform, "keep_constant"=df$keep_constant,
         "a_n"=a_n, "der_iden"=der_iden, "names"=df$names,"Permutation"=df$current_order)
}

#' Automatically assigns missing values in listed columns
#' \code{Replace_Missing} checks each column and fills in NA values
#'
#' @inheritParams R_template
#'
#' @return returns a filled datatable
#' @export
#' @examples
#' library(data.table)
#' ## basic example code reproduced from the starting-description vignette
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  NA,  47,  36,  NA,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0))
#' df <- Replace_Missing(df, c("Starting_Age","Ending_Age"), 70)
Replace_Missing <- function(df,name_list,MSV,verbose=FALSE){
    if (is.na(MSV)){
        if (verbose){
            message("Error: The missing-value replacement is also NA")
        }
        stop()
    }
    for (j in name_list){
        #
        if (j %in% names(df)){
            #fine
        } else {
            if (verbose){
                message("Error: ",paste(j," missing from column names",sep=""))
            }
            stop()
        }
        #
        if (sum(is.na(df[[j]]))){
            set(df,which(is.na(df[[j]])),j,MSV)
            if (verbose){
                message(paste("Note: Column ",j," had replaced values",sep=""))
            }
        }
    }
    return (df)
}


#' Automatically assigns missing control values
#' \code{Def_Control} checks and assigns default values
#'
#' @inheritParams R_template
#'
#' @return returns a filled list
#' @export
#' @examples
#' library(data.table)
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5, 'ties'='breslow','double_step'=1)
#' control <- Def_Control(control)
#'
Def_Control <- function(control){
    control_def <- list('verbose'=FALSE,'lr' = 0.75,'maxiter' = 20,'maxiters' = c(1,20),
        'guesses'=1,'halfmax' = 5,'epsilon' = 1e-9, 'dbeta_max' = 0.5,
        'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
        'ties'='breslow','double_step'=1,"keep_strata"=FALSE,
        "Ncores"=as.numeric(detectCores()), "cens_thres"=0)
	if ((identical(Sys.getenv("TESTTHAT"), "true"))||(identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))){
		control_def$Ncores <- 2
	}
    for (nm in names(control_def)){
        if (nm %in% names(control)){
            if (nm=="Ncores"){
                if (control$Ncores>control_def$Ncores){
                    if (control$verbose){
                        message(paste("Error: Cores Requested:",control["Ncores"],
                              ", Cores Avaliable:",control_def["Ncores"],sep=" "))
                    }
                    stop()
                }
            }
        } else {
            control[nm] <- control_def[nm]
        }
    }
    return (control)
}

#' Automatically assigns geometric-mixture values
#' \code{Def_model_control} checks and assigns default values
#'
#' @inheritParams R_template
#'
#' @return returns a filled list
#' @export
#' @examples
#' library(data.table)
#' control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 5, 'ties'='breslow','double_step'=1)
#' control <- Def_Control(control)
#' model_control <- list("single"=TRUE)
#' model_control <- Def_model_control(model_control)
#' Term_n <- c(0,1,1)
#' modelform <- 'a'
#' val <- Def_modelform_fix(control,model_control,modelform,Term_n)
#' model_control <- val$model_control
#' modelform <- val$modelform
#'
Def_modelform_fix <- function(control,model_control,modelform,Term_n){
    term_tot <- max(Term_n)+1
    modelform <- toupper(modelform)
    acceptable <- c('A','PA','PAE','M','ME','GMIX','GMIX-R','GMIX-E')
    if (modelform %in% acceptable){
        if (modelform=='ME'){
            modelform <- 'M'
        } else if (modelform=="GMIX-R"){
            model_control$gmix_term <- rep(0,term_tot)
            modelform <- "GMIX"
        } else if (modelform=="GMIX-E"){
            model_control$gmix_term <- rep(1,term_tot)
            modelform <- "GMIX"
        }
    } else {
        if (control$verbose){
            message(paste("Error: Model formula ",modelform," not in acceptable list",sep=""))
        }
        stop()
    }
    if (modelform == "GMIX"){
        gmix_term <- model_control$gmix_term
        if (length(gmix_term) != term_tot){
            if (control$verbose){
                message(paste("Error: Terms used:",term_tot,", Terms with gmix types avaliable:",length(gmix_term),sep=" "))
            }
            stop()
        }
    }
    return(list('modelform'=modelform,'model_control'=model_control))
}

#' Automatically assigns missing model control values
#' \code{Def_model_control} checks and assigns default values
#'
#' @inheritParams R_template
#'
#' @return returns a filled list
#' @export
#' @examples
#' library(data.table)
#' control <- list("single"=TRUE)
#' control <- Def_model_control(control)
#'
Def_model_control <- function(control){
    control_def_names <- c('single','basic','null','CR','strata','Surv','Schoenfeld','Risk','Risk_Subset')
    for (nm in control_def_names){
        if (nm %in% names(control)){
            #fine
        } else {
            control[nm] <- FALSE
        }
    }
    if ("Unique_Values" %in% names(control)){
        #fine
    } else {
        control["Unique_Values"] <- 2
    }
    if ("gmix_theta" %in% names(control)){
        #fine
    } else {
        control["gmix_theta"] <- 0.5
    }
    if ("gmix_term" %in% names(control)){
        #fine
    } else {
        control["gmix_term"] <- c(0)
    }
    return (control)
}

#' Automatically assigns missing guessing control values
#' \code{Def_Control_Guess} checks and assigns default values
#'
#' @inheritParams R_template
#'
#' @return returns a filled list
#' @export
#' @examples
#' library(data.table)
#' guesses_control <- list("maxiter"=10,"guesses"=10,
#'     "loglin_min"=-1,"loglin_max"=1,"loglin_method"="uniform")
#' a_n <- c(0.1,2,1.3)
#' guesses_control <- Def_Control_Guess(guesses_control, a_n)
#'
Def_Control_Guess <- function(guesses_control, a_n){
    if ("verbose" %in% names(guesses_control)){ #determines extra printing
        #fine
    } else {
        guesses_control$verbose <- FALSE
    }
    if ("maxiter" %in% names(guesses_control)){ #determines the iterations for each guess
        #fine
    } else {
        guesses_control$maxiter <- 5
    }
    if ("guesses" %in% names(guesses_control)){ #determines the number of guesses
        #fine
    } else {
        guesses_control$guesses <- 10
    }
    # -------------------------------------------------------------------- #
    if ("exp_min" %in% names(guesses_control)){ #minimum exponential parameter change
        #fine
    } else {
        guesses_control$exp_min <- -1
    }
    if ("exp_max" %in% names(guesses_control)){ #maximum exponential parameter change
        #fine
    } else {
        guesses_control$exp_max <- 1
    }
    # -------------------------------------------------------------------- #
    if ("intercept_min" %in% names(guesses_control)){ #minimum intercept parameter change
        #fine
    } else {
        guesses_control$intercept_min <- -1
    }
    if ("intercept_max" %in% names(guesses_control)){ #maximum intercept parameter change
        #fine
    } else {
        guesses_control$intercept_max <- 1
    }
    # -------------------------------------------------------------------- #
    if ("lin_min" %in% names(guesses_control)){ #minimum linear parameter change
        #fine
    } else {
        guesses_control$lin_min <- 0.1
    }
    if ("lin_max" %in% names(guesses_control)){ #maximum linear parameter change
        #fine
    } else {
        guesses_control$lin_max <- 1
    }
    # -------------------------------------------------------------------- #
    if ("exp_slope_min" %in% names(guesses_control)){ #minimum exp_slope parameter change
        #fine
    } else {
        guesses_control$exp_slope_min <- 0.001
    }
    if ("exp_slope_max" %in% names(guesses_control)){ #maximum exp_slope parameter change
        #fine
    } else {
        guesses_control$exp_slope_max <- 2
    }
    # -------------------------------------------------------------------- #
    if ("term_initial" %in% names(guesses_control)){ #list of term numbers for tiered guessing
        #fine
    } else {
        guesses_control$term_initial <- c(0)
    }
    #
    if ("guess_constant" %in% names(guesses_control)){ 
        #binary values for if a parameter is distributed (0) or not (1)
        if (length(guesses_control$guess_constant)<length(a_n)){
            guesses_control$guess_constant <- c(guesses_control$guess_constant, 
                rep(0,length(a_n)-length(guesses_control$guess_constant)))
        }
    } else {
        guesses_control$guess_constant <- rep(0,length(a_n))
    }
    if ("guesses_start" %in% names(guesses_control)){
        #number of guesses for first part of tiered guessing
        #fine
    } else {
        guesses_control$guesses_start <- guesses_control$guesses
    }
    if ("strata" %in% names(guesses_control)){ #if stratification is used
        #fine
    } else {
        guesses_control$strata <- FALSE
    }
    if (("rmin" %in% names(guesses_control))&&("rmax" %in% names(guesses_control))){
        #if rmin/rmax used is used
        #fine
    } else {
        # if they are unequal length, the default is used
        guesses_control$rmin <- c(-1)
        guesses_control$rmax <- c(-1,-1)
    }
    return (guesses_control)
}


#' Calculates Full Parameter list for Special Dose Formula
#' \code{Linked_Dose_Formula} Calculates all parameters for linear-quadratic and linear-exponential linked formulas
#'
#' @inheritParams R_template
#'
#' @return returns list of full parameters
#' @export
#' @examples
#' library(data.table)
#' tforms <- list("cov_0"="quad", "cov_1"="exp")
#' paras <- list("cov_0"=c(1,3.45), "cov_1"=c(1.2, 4.5, 0.1))
#' full_paras <- Linked_Dose_Formula(tforms, paras)
#'
Linked_Dose_Formula <- function(tforms,paras,verbose=FALSE){
    full_paras <- list()
    for (nm in names(tforms)){
        if (tforms[nm]=="quad"){
            plist <- unlist(paras[nm],use.names=FALSE)
            a0 <- plist[1]
            y <- plist[2]
            if (a0 < 0 ){
                if (verbose){
                    message("Error: a0 arguement was negative")
                }
                stop()
            }
            if (is.numeric(a0)){
                #fine
            } else {
                if (verbose){
                    message("Error: a0 arguement was not a number")
                }
                stop()
            }
            if (is.numeric(y)){
                #fine
            } else {
                if (verbose){
                    message("Error: threshold arguement was not a number")
                }
                stop()
            }
            a1 <- a0/2/y
            b1 <- a0*y/2
            full_paras[[nm]] <- c(y,a0,a1,b1)
        } else if (tforms[nm]=="exp"){
            plist <- unlist(paras[nm],use.names=FALSE)
            a0 <- plist[1]
            y <- plist[2]
            b1 <-plist[3]
            if (a0 < 0 ){
                if (verbose){
                    message("Error: a0 arguement was negative")
                }
                stop()
            }
            if (is.numeric(a0)){
                #fine
            } else {
                if (verbose){
                    message("Error: a0 arguement was not a number")
                }
                stop()
            }
            if (is.numeric(y)){
                #fine
            } else {
                if (verbose){
                    message("Error: threshold arguement was not a number")
                }
                stop()
            }
            if (is.numeric(b1)){
                #fine
            } else {
                if (verbose){
                    message("Error: exponential arguement was not a number")
                }
                stop()
            }
            c1 <- log(a0)-log(b1)+b1*y
            a1 <- a0*y+exp(c1-b1*y)
            full_paras[[nm]] <- c(y,a0,a1,b1,c1)
        }
    }
    return (full_paras)
}

#' Calculates The Additional Parameter For a linear-exponential formula with known maximum
#' \code{Linked_Lin_Exp_Para} Calculates what the additional parameter would be for a desired maximum
#'
#' @inheritParams R_template
#'
#' @return returns parameter used by Colossus
#' @export
#' @examples
#' library(data.table)
#' y <- 7.6
#' a0 <- 1.2
#' a1_goal <- 15
#' full_paras <- Linked_Lin_Exp_Para(y,a0,a1_goal)
#'
Linked_Lin_Exp_Para <- function(y,a0,a1_goal,verbose=FALSE){
    b1 <- 10
    lr <- 1.0
    if (a0 < 0 ){
        if (verbose){
            message("Error: a0 arguement was negative")
        }
        stop()
    }
    if (a1_goal > y*a0){
        #fine
    } else {
        if (verbose){
            message("Error: goal is too low")
        }
        stop()
    }
    iter_i <- 0
    while (iter_i<100){
        iter_i <- iter_i + 1
        c1 <- log(a0/b1)+b1*y
        a1 <- a0*y+exp(c1-b1*y)
        a_dif <- a1-a1_goal
        if (abs(a_dif)<1e-3){
            break   
        }
        da1 <- -1/b1*exp(c1-b1*y)
        db1 <- (a1_goal-a1)/da1
        if (-1*db1 > b1){
            db1 <- -0.9*b1
        }
        b1 <- b1 + lr*db1
    }
    return (b1)
}

#' Splits a parameter into factors
#' \code{factorize} uses user provided list of columns to define new parameter for each unique value and update the data.table.
#' Not for interaction terms
#'
#' @inheritParams R_template
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#' @examples
#' library(data.table)
#' a <- c(0,1,2,3,4,5,6)
#' b <- c(1,2,3,4,5,6,7)
#' c <- c(0,1,2,1,0,1,0)
#' df <- data.table("a"=a,"b"=b,"c"=c)
#' col_list <- c("c")
#' val <- factorize(df,col_list)
#' df <- val$df
#' new_col <- val$cols
#'
factorize <-function(df,col_list,verbose=FALSE){
    cols <- c()
    col0 <- names(df)
    tnum <- c()
    for (i in seq_len(length(col_list))){
        col <- col_list[i]
        x <- sort(unlist(as.list(unique(df[,col, with = FALSE])),use.names=FALSE))
        for (j in x){
            newcol <- c(paste(col,j,sep="_"))
            df[, newcol] <- 1*(df[,col, with = FALSE]==j)
            cols <- c(cols, newcol)
            tnum <- c(tnum, i)
        }
    }
    #
    cols <- setdiff(names(df), col0)
    #
    cols <- Check_Dupe_Columns(df,cols,rep(0,length(cols)),verbose,TRUE)
    if (verbose){
        message(paste("Note: Number of factors:",length(cols),sep=""))
    }
    list('df'=df, 'cols'=cols)
}


#' Splits a parameter into factors in parallel
#' \code{factorize_par} uses user provided list of columns to define new parameter for each unique value and update the data.table.
#' Not for interaction terms
#'
#' @inheritParams R_template
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#' @examples
#' library(data.table)
#' a <- c(0,1,2,3,4,5,6)
#' b <- c(1,2,3,4,5,6,7)
#' c <- c(0,1,2,1,0,1,0)
#' df <- data.table("a"=a,"b"=b,"c"=c)
#' col_list <- c("c")
#' val <- factorize_par(df,col_list,FALSE,2)
#' df <- val$df
#' new_col <- val$cols
#'
factorize_par <-function(df,col_list,verbose=FALSE, nthreads=as.numeric(detectCores())){
    cols <- c()
    vals <- c()
    names <- c()
    if ((identical(Sys.getenv("TESTTHAT"), "true"))||(identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))){
		nthreads <- 2
	}
    for (i in seq_len(length(col_list))){
        col <- col_list[i]
        x <- sort(unlist(as.list(unique(df[,col, with = FALSE])),use.names=FALSE))
        for (j in x){
            newcol <- c(paste(col,j,sep="_"))
            ##
            names <- c(names,newcol)
            vals <- c(vals, j)
            cols <- c(cols, i-1)
            ##
        }
    }
    #
    df0 <- Gen_Fac_Par(as.matrix(df[, col_list,with=FALSE]), vals, cols, nthreads)
    df0 <- as.data.table(df0)
    names(df0) <- names
    #
    col_keep <- Check_Dupe_Columns(df0,names,rep(0,length(cols)),verbose,TRUE)
    #
    #
    list('df'=cbind(df,df0), 'cols'=col_keep)
}


#' Defines Interactions
#' \code{interact_them} uses user provided interactions define interaction terms and update the data.table. assumes interaction is "+" or "*" and applies basic anti-aliasing to avoid duplicates
#'
#' @inheritParams R_template
#'
#' @return returns a list with two named fields. df for the updated dataframe, and cols for the new column names
#' @export
#' @examples
#' library(data.table)
#' a <- c(0,1,2,3,4,5,6)
#' b <- c(1,2,3,4,5,6,7)
#' c <- c(0,1,2,1,0,1,0)
#' df <- data.table("a"=a,"b"=b,"c"=c)
#' interactions <- c("a?+?b","a?*?c")
#' new_names <- c("ab","ac")
#' vals <- interact_them(df, interactions, new_names)
#' df <- vals$df
#' new_col <- vals$cols
#' 
interact_them <- function(df,interactions,new_names,verbose=FALSE){
    cols <- c()
    for (i in seq_len(length(interactions))){
        interac <- interactions[i]
        formula <- unlist(strsplit(interac,"[?]"),use.names=FALSE)
        if (length(formula)!=3){
            if (verbose){
                message(paste("Error: Iteration:",interac,"has incorrect length of",
                      length(formula),"but should be 3."))
            }
            stop()
        }
        newcol <- paste(formula[1],formula[2],formula[3],sep="")
        if (new_names[i]!=''){
            newcol <- new_names[i]
        }
        col1 <- formula[1]
        col2 <- formula[3]
        if (paste(formula[1],"?",formula[2],"?",formula[3],sep="") %in% interactions[i+seq_len(length(interactions))]){
            "its duped"
        } else if (paste(formula[3],"?",formula[2],"?",formula[1],sep="") %in% interactions[i+seq_len(length(interactions))]){
            "the reverse is duped"
        } else {
            if (formula[2]=="+"){
                df[, newcol] <- df[,col1, with = FALSE] + df[,col2, with = FALSE]
                cols <- c(cols, newcol)
            } else if (formula[2]=="*"){
                df[, newcol] <- df[,col1, with = FALSE] * df[,col2, with = FALSE]
                cols <- c(cols, newcol)
            } else {
                if (verbose){
                    message(paste("Error: Incorrect operation of",formula[2]))
                }
                stop()
            }
        }
    }
    list('df'=df, 'cols'=cols)
}

#' Defines the likelihood ratio test
#' \code{Likelihood_Ratio_Test} uses two models and calculates the ratio
#'
#' @inheritParams R_template
#'
#' @return returns the score statistic
#' @export
#' @examples
#' library(data.table)
#' #In an actual example, one would run two seperate RunCoxRegression regressions,
#' #    assigning the results to e0 and e1
#' e0 <- list("name"="First Model","LogLik"=-120)
#' e1 <- list("name"="New Model","LogLik"=-100)
#' score <- Likelihood_Ratio_Test(e1, e0)
#'
Likelihood_Ratio_Test <- function(alternative_model, null_model){
    if (("LogLik" %in% names(alternative_model))&&("LogLik" %in% names(null_model))){
        return (2*(unlist(alternative_model["LogLik"],use.names=FALSE) - unlist(null_model["LogLik"],use.names=FALSE)))
    }
    stop()
    return (NULL)
}



#' checks for duplicated column names
#' \code{Check_Dupe_Columns} checks for duplicated columns, columns with the same values, and columns with 1 value. Currently not updated for multi-terms
#'
#' @inheritParams R_template
#'
#' @return returns the usable columns
#' @export
#' @examples
#' library(data.table)
#' a <- c(0,1,2,3,4,5,6)
#' b <- c(1,2,3,4,5,6,7)
#' c <- c(0,1,2,1,0,1,0)
#' df <- data.table("a"=a,"b"=b,"c"=c)
#' cols <- c("a","b","c")
#' term_n <- c(0,0,1)
#' unique_cols <- Check_Dupe_Columns(df, cols, term_n)
#'
Check_Dupe_Columns <- function(df,cols,term_n,verbose=FALSE, factor_check=FALSE){
    ##
    if (length(cols)>1){
        features_pair <- combn(cols, 2, simplify = FALSE) # list all column pairs
        terms_pair <- combn(term_n, 2, simplify = FALSE) # list all term pairs
        toRemove <- c() # init a vector to store duplicates
        for(pair_n in seq_len(length(features_pair))) { 
            # put the pairs for testing into temp objects
            pair <- unlist(features_pair[pair_n])
            term <- unlist(terms_pair[pair_n])
            f1 <- pair[1]
            f2 <- pair[2]
            if (!(f1 %in% names(df))){
                if (verbose){
                    message(paste("Error: ",f1," not in data.table",sep=""))
                }
                stop()
            }
            if (!(f2 %in% names(df))){
                if (verbose){
                    message(paste("Error: ",f2," not in data.table",sep=""))
                }
                stop()
            }
            t1 <- term[1]
            t2 <- term[2]
            #
            checked_factor <- TRUE
            if (factor_check){ #
                if ((is.numeric(df[[f1]]))&(is.numeric(df[[f2]]))){
                    if (sum(df[[f1]]*df[[f2]])==0){
                        checked_factor <- FALSE
                    }
                } else if (is.numeric(df[[f1]]) != is.numeric(df[[f2]])){
                    checked_factor <- FALSE
                }
            }
            if ((t1==t2)&(checked_factor)){
                df[,get(f1)]
                df[,get(f2)]
                if (!(f1 %in% toRemove) & !(f2 %in% toRemove)) {
                    if (all(df[[f1]] == df[[f2]])) { # test for duplicates
                        if (verbose){
                            message(paste("Warning: ",f1, " and ", f2, " are equal",sep=""))
                        }
                        toRemove <- c(toRemove, f2) # build the list of duplicates
                    }
                    if (min(df[[f2]])==max(df[[f2]])){
                        if (min(df[[f2]])==0){
                            toRemove <- c(toRemove, f2) # remove zero values
                        }
                    }
                }
            }
        }
        newcol <- setdiff(cols, toRemove)
        if (length(newcol)==1){
            if (min(df[,newcol, with = FALSE])==max(df[,newcol, with = FALSE])){
                return(c())
            } else {
                return(newcol)
            }
        }
        return(newcol)
    } else if (length(cols)==1){
        f1 <- cols[1]
        if (!(f1 %in% names(df))){
            if (verbose){
                message(paste("Error: ",f1," not in data.table",sep=""))
            }
            stop()
        }
        if (min(df[,cols, with = FALSE])==max(df[,cols, with = FALSE])){
            if (min(df[,cols, with = FALSE])==0){
                return(c())
            } else {
                return(cols)
            }
        } else {
            return(cols)
        }
    } else {
        return(c())
    }
    return(c())
}

#' Applies time duration truncation limits
#' \code{Check_Trunc} creates columns to use for truncation
#'
#' @inheritParams R_template
#'
#' @return returns the updated data and time period columns
#' @export
#' @examples
#' library(data.table)
#' 
#' df <- data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
#'            "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
#'              "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
#'           "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0))
#' # For the interval case
#' time1 <- "Starting_Age"
#' time2 <- "Ending_Age"
#' ce <- c("%trunc%","Ending_Age")
#' val <- Check_Trunc(df, ce)
#' df <- val$df
#' ce <- val$ce
#'
Check_Trunc <- function(df,ce,verbose=FALSE){
    if (ce[1]=="%trunc%"){
        if (ce[2]=="%trunc%"){
            if (verbose){
                message("Error: Both endpoints are truncated, not acceptable")
            }
            stop()
        }
        tname <- ce[2]
        tmin <- min(df[,get(tname)])-1
        if (!("right_trunc" %in% names(df))){
            df[,':='(right_trunc=tmin)]
        }
        ce[1] <- "right_trunc"
    } else if (ce[2]=="%trunc%") {
        tname <- ce[1]
        tmax <- max(df[,get(tname)])+1
        if (!("left_trunc" %in% names(df))){
            df[,':='(left_trunc=tmax)]
        }
        ce[2] <- "left_trunc"
    }
    return (list('df'=df,'ce'=ce))
}

#' Applies time dependence to parameters
#' \code{gen_time_dep} generates a new dataframe with time dependent covariates by applying a grid in time
#'
#' @inheritParams R_template
#'
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#' #Adapted from the tests
#' a <- c(20,20,5,10,15)
#' b <- c(1,2,1,1,2)
#' c <- c(0,0,1,1,1)
#' df <- data.table("a"=a,"b"=b,"c"=c)
#' time1 <- "%trunc%"
#' time2 <- "a"
#' event <- "c"
#' control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,
#'            'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
#'            'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#' grt_f <- function(df,time_col){
#'     return ((df[,"b"] * df[,get(time_col)])[[1]])
#' }
#' func_form <- c("lin")
#' df_new <- gen_time_dep(df,time1,time2,event,TRUE,0.01,c("grt"),c(),
#'        c(grt_f),paste("test","_new.csv",sep=""), func_form,2)
#' file.remove('test_new.csv')
#'
gen_time_dep <- function(df, time1, time2, event0, iscox, dt, new_names, dep_cols, func_form,fname, tform,nthreads=as.numeric(detectCores())){
    dfn <- names(df)
    ce <- c(time1,time2,event0)
    t_check <- Check_Trunc(df,ce)
    df <- t_check$df
    ce <- t_check$ce
    time1 <- ce[1]
    time2 <- ce[2]
    dfn_same <- dfn[!(dfn %in% dep_cols)]
    dfn_dep <- c()
    for (i in seq_len(length(new_names))){
        name0 <- paste(new_names[i],0,sep="_")
        name1 <- paste(new_names[i],1,sep="_")
        func <- func_form[i]
        df[, name0] <- lapply(func, function(f) f(df, time1))
        df[, name1] <- lapply(func, function(f) f(df, time2))
        dfn_dep <- c(dfn_dep, name0, name1)
    }
    #
    #
    dfn_time <- c(time1, time2)
    dfn_event <- c(event0)
    dfn_same <- dfn_same[!(dfn_same %in% dfn_time)]
    dfn_same <- dfn_same[!(dfn_same %in% dfn_event)]
    ##
    dfend <- df[get(event0)==1, ]
    tu <- sort(unlist(unique(dfend[,time2, with = FALSE]), use.names=FALSE))
    if (iscox){
        #
        df <- df[get(time2)>=min(tu),]
        df <- df[get(time1)<=max(tu),]
        #
    }
    ##
    x_time <- as.matrix(df[,dfn_time, with = FALSE])
    x_dep <- as.matrix(df[,dfn_dep, with = FALSE])
    x_same <- as.matrix(df[,dfn_same, with = FALSE])
    x_event <- as.matrix(df[,dfn_event, with = FALSE])
    #
    if (grepl(".csv", fname, fixed = TRUE)){
        #fine
    } else {
        fname <- paste(fname,".csv",sep="_")
    }
    if ((identical(Sys.getenv("TESTTHAT"), "true"))||(identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true"))){
		nthreads <- 2
	}
    Write_Time_Dep(x_time, x_dep, x_same, x_event, dt, fname,tform,tu,iscox, nthreads)
    df_new <- fread(fname,data.table=TRUE,header=FALSE,nThread=nthreads,
                    col.names=c(time1,time2,new_names,dfn_same,event0))
    setkeyv(df_new, c(time2, event0))
    return (df_new)
}

#' Automates creating a date difference column
#' \code{Date_Shift} generates a new dataframe with a column containing time difference in a given unit
#'
#' @inheritParams R_template
#' @param dcol0 list of starting month, day, and year
#' @param dcol1 list of ending month, day, and year
#'
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#' m0 <- c(1,1,2,2)
#' m1 <- c(2,2,3,3)
#' d0 <- c(1,2,3,4)
#' d1 <- c(6,7,8,9)
#' y0 <- c(1990,1991,1997,1998)
#' y1 <- c(2001,2003,2005,2006)
#' df <- data.table("m0"=m0,"m1"=m1,"d0"=d0,"d1"=d1,"y0"=y0,"y1"=y1)
#' df <- Date_Shift(df,c("m0","d0","y0"),c("m1","d1","y1"),"date_since")
#'
Date_Shift <- function(df, dcol0, dcol1, col_name, units="days"){
    def_cols <- names(df)
    #
    df$dt0 <- paste(df[[match(dcol0[1],names(df))]],df[[match(dcol0[2],names(df))]],
              df[[match(dcol0[3],names(df))]],sep="-")
    df$dt1 <- paste(df[[match(dcol1[1],names(df))]],df[[match(dcol1[2],names(df))]],
              df[[match(dcol1[3],names(df))]],sep="-")
    #
    # TO NOT ENCOUNTER DAYLIGHT SAVINGS ISSUES, THE UTC TIMEZONE IS USED
    # IF NOT USED THEN RESULTS MAY HAVE TIMES OFF BY 1/24 decimals
    df[, col_name] <- difftime(strptime(df$dt1, format = "%m-%d-%Y",tz = 'UTC'), 
                   strptime(df$dt0,  format = "%m-%d-%Y"), units = units,tz = 'UTC')
    def_cols <- c(def_cols, col_name)
    return (df[,def_cols,with=FALSE])
}

#' Automates creating a date since a reference column
#' \code{Time_Since} generates a new dataframe with a column containing time since a reference in a given unit
#'
#' @inheritParams R_template
#' @param dcol0 list of ending month, day, and year
#'
#' @return returns the updated dataframe
#' @export
#' @examples
#' library(data.table)
#'m0 <- c(1,1,2,2)
#'m1 <- c(2,2,3,3)
#'d0 <- c(1,2,3,4)
#'d1 <- c(6,7,8,9)
#'y0 <- c(1990,1991,1997,1998)
#'y1 <- c(2001,2003,2005,2006)
#'df <- data.table("m0"=m0,"m1"=m1,"d0"=d0,"d1"=d1,"y0"=y0,"y1"=y1)
#'tref <- strptime( "3-22-1997", format = "%m-%d-%Y",tz = 'UTC')
#'df <- Time_Since(df,c("m1","d1","y1"),tref,"date_since")
#'
Time_Since <- function(df, dcol0, tref, col_name, units="days"){
    def_cols <- names(df)
    #
    df$dt0 <- paste(df[[match(dcol0[1],names(df))]],df[[match(dcol0[2],
                    names(df))]],df[[match(dcol0[3],names(df))]],sep="-")
    #
    #
    df[, col_name] <- lapply(df$dt0, function(x) (difftime(strptime(x,
                             format = "%m-%d-%Y",tz = 'UTC'), tref, units = units)))
    def_cols <- c(def_cols, col_name)
    return (df[,def_cols,with=FALSE])
}

