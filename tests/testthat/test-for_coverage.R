test_that("Gather Guesses list, incorrect keep_constant length and rmin/rmax not used", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d","d","d","d")
    Term_n <- c(0,0,0,0)
    tform <- c("loglin",'loglin','loglin','loglin')
    keep_constant <- c(0,0)
    a_n <- list(c(-0.1,6,0.1,0.1))
    a_n_default <- unlist(a_n[1])
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    #
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    guesses_control <- list()
    model_control <- list()
    
    all_names <- unique(names(df))
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n_default)
    guesses_control$verbose <- TRUE
    model_control <- Def_model_control(model_control)
    #
    expect_no_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
    keep_constant <- c(1,0,0,0,0,0,0)
    guesses_control$rmin <- c(-0.1,-1,-0.1,0)
    guesses_control$rmax <- c(0.1, 1, 0.1, 0.1)
    expect_no_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
})
test_that("Gather Guesses list, bad tform", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d","d","d","d")
    Term_n <- c(0,0,0,0)
    tform <- c("loglin",'loglin','loglin','bad_bad')
    keep_constant <- c(0,0)
    a_n <- list(c(-0.1,6,0.1,0.1))
    a_n_default <- unlist(a_n[1])
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    #
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    guesses_control <- list()
    model_control <- list()
    
    all_names <- unique(names(df))
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n_default)
    guesses_control$verbose <- TRUE
    model_control <- Def_model_control(model_control)
    #
    guesses_control$rmin <- c(-0.1,-1,-0.1,0)
    guesses_control$rmax <- c(0.1, 1, 0.1, 0.1,0,0,0)
    expect_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
})


test_that("tform order, der_iden out of bounds", {
    Term_n <- c(0,0,0,0,0)
    tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
    keep_constant <- c(0,0,0,0,0)
    a_n <- c(1,2,3,4,5)
    names <- c("a","a","a","a","a")
    expect_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,der_iden=-1))
    expect_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,der_iden=100))
})
test_that("tform order, matrix errors", {
    Term_n <- c(0,0,0,0,0)
    tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
    keep_constant <- c(0,0,0,0,0)
    a_n <- c(1,2,3,4,5)
    names <- c("a","a","a","a","a")
    Cons_Mat <- matrix(c(1,2,3,4,5),nrow=1,byrow=T)
    Cons_Vec <- c(1)

    expect_no_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
    Cons_Mat <- matrix(c(1,2,3),nrow=1,byrow=T)
    expect_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
    Cons_Mat <- matrix(c(1,2,3,4,5),nrow=1,byrow=T)
    Cons_Vec <- c(1,1,1)
    expect_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
})
test_that("tform order, keep_constant errors", {
    Term_n <- c(0,0,0,0,0)
    tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
    keep_constant <- c(0,0,0)
    a_n <- c(1,2,3,4,5)
    names <- c("a","a","a","a","a")
    Cons_Mat <- matrix(c(1,2,3,4,5),nrow=1,byrow=T)
    Cons_Vec <- c(1)
    expect_no_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
    keep_constant <- c(0,0,0,0,0,0,0,0)
    expect_no_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
    #
    keep_constant <- c(0,0,0)
    a_n <- list(c(1,2,3,4,5),c(4,2,3,4,5),c(1,2,7,4,5))
    names <- c("a","a","a","a","a")
    Cons_Mat <- matrix(c(1,2,3,4,5),nrow=1,byrow=T)
    Cons_Vec <- c(1)
    expect_no_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
    keep_constant <- c(0,0,0,0,0,0,0,0)
    expect_no_error(Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec))
})

test_that("Missing Value verbose error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_no_error(Replace_Missing(df,c("a","b","c","d"),0.0,T))
    expect_error(Replace_Missing(df,c("a","b","c","d"),0.0,-1))
})

test_that("Pois various_fixes", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    df$pyr <- df$t1-df$t0
	pyr <- "pyr"
    event <- "lung"
    df$rand <- floor(runif(nrow(df), min=0, max=5))
    names <- c("dose","rand","rand")
    Term_n <- c(2,1,0)
    tform <- c("loglin","loglin","loglin")
    keep_constant <- c(0,0,0)
    a_n <- c(0.01,0.1,0.1)
    modelform <- "PAE"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=0)
    Strat_Col <- "fac"
    
    verbose <- FALSE
    model_control=list('strata'=F, 'single'=F)
    if (verbose){print(model_control)}
    a_n <- c(0.01,0.1,0.1)
    modelform <- "PAE"
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    keep_constant <- c(1,1,1)
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    keep_constant <- c(0,0,0)
    ev0 <- df$lung
    df$lung <- rep(0,length(ev0))
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    names <- c("dose","rand","CONST")
    df$lung <- ev0
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    control$guesses <- 100
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
})
test_that("Pois_tier_guess various_fixes", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    df$pyr <- df$t1-df$t0
	pyr <- "pyr"
    event <- "lung"
    names <- c("t0","a")
    Term_n <- c(1,2)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,-15)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=TRUE,'term_initial' = c(0,1),'verbose'=TRUE)
    Strat_Col=c('b')

    expect_no_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, Strat_Col))
    keep_constant <- c(1,1)
    expect_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, Strat_Col))
})
test_that("Poisson_basic_guess_cpp various_fixes", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    pyr <- "t1"
    event <- "lung"
    names <- c("a","b")
    Term_n <- c(0,1)
    tform <- c("loglin","loglin")
    Strat_Col <- 'a'
    keep_constant <- c(0,0)
    a_n <- c(0.01,-15)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=1,"guesses"=1,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=T)
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    model_control <- list('strata'=T)
    expect_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    guesses_control=list("Iterations"=1,"guesses"=1,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'term_initial' = c(0,1),'verbose'=T)
    model_control <- list('strata'=F)
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    keep_constant <- c(1,1)
    expect_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    keep_constant <- c(0,0)
    names <- c("a","CONST")
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
})

test_that("Coxph Martingale no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,2,2,3,3,3)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=e)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=TRUE,"cov_cols"="d","surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="Not_In",'verbose'=TRUE)
    if (system.file(package='ggplot2')!=""){
        expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    }
})

test_that("Coxph loglin_M CENSOR Default various_fixes", {
    fname <- 'll_cens_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","fac")
    Term_n <- c(0,0)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0,0)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = -1,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    plot_options <- list("name"=paste(tempfile(),"run",sep=""),"verbose"=T,"studyID"="studyID","age_unit"="years")
    expect_no_error(GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    keep_constant <- c(1,1)
    expect_error(GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    keep_constant <- c(0,0)
    df$lung <- rep(0,nrow(df))
    expect_error(GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    #
})
test_that("Coxph risk plotting above discrete step number limit", {
    a <- rep(c(0,1,2,3,4,5,6),20)
    b <- rep(c(1,2,3,4,5,6,7),20)
    c <- rep(c(1,0,1,0,1,0,0),20)
    d <- runif(length(a))
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("RISK",paste(tempfile(),"run",sep="")),"studyID"="a",'verbose'=TRUE)
    if (system.file(package='ggplot2')!=""){
        expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    }
})

test_that("Various CoxRegressionOmnibus options", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    set.seed(3742)
    df$rand <- floor(runif(nrow(df), min=0, max=5))

    time1 <- "t0"
    time2 <- "t1"


    #
    event <- "lung"
    names <- c("rand","fac","dose")
    Term_n <- c(0,0,1)
    tform <- c("loglin","loglin","loglin")
    keep_constant <- c(0,0,0)
    a_n <- c(-0.1,0.1,0.2)
    modelform <- "M"
    fir <- 0
    der_iden <- 0

    cens_weight <- c(0)

    verbose <- FALSE

    devs <- c()

    modelform <- "M"
    model_control=list('strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'CR'=FALSE)
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0,0,0)
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    keep_constant <- c(1,1,1)
    expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    lung_temp <- df$lung
    df$lung <- rep(0,length(lung_temp))
    keep_constant <- c(0,0,0)
    expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    df$lung <- lung_temp
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1,1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    a_n <- list(c(0.6465390, 0.4260961, 0.1572781),c(0.6465390, 0.4260961, 0.1572781),c(0.6465390, 0.4260961, 0.1572781))
    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=50)
    expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    #
    control=list("Ncores"=2,'lr' = 0.75,'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=1)
    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    control=list("Ncores"=2,'lr' = 0.75,'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=10)
    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
    #
    names <- c("rand","fac","dose")
    Term_n <- c(0,0,1)
    tform <- c("lin","lin","lin")
    keep_constant <- c(0,0,0)
    a_n <- c(-0.1,-0.1,0.2)
    expect_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
})

test_that("Various RunPoissonRegression_Omnibus options", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    set.seed(3742)
    df$rand <- floor(runif(nrow(df), min=0, max=5))

    time1 <- "t0"
    time2 <- "t1"
	df$pyr <- df$t1-df$t0
	pyr <- "pyr"

    #
    event <- "lung"
    names <- c("rand","fac","dose")
    Term_n <- c(0,0,1)
    tform <- c("loglin","loglin","loglin")
    keep_constant <- c(0,0,0)
    a_n <- c(-0.1,0.1,0.2)
    modelform <- "M"
    fir <- 0
    der_iden <- 0

    cens_weight <- c(0)

    verbose <- FALSE

    devs <- c()

    modelform <- "M"
    model_control=list('strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'CR'=FALSE)
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0,0,0)
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    Strat_Col="fac"
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    keep_constant <- c(1,1,1)
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    lung_temp <- df$lung
    df$lung <- rep(0,length(lung_temp))
    keep_constant <- c(0,0,0)
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    df$lung <- lung_temp
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1,1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    a_n <- list(c(0.6465390, 0.4260961, 0.1572781),c(0.6465390, 0.4260961, 0.1572781),c(0.6465390, 0.4260961, 0.1572781))
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=50)
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    #
    control=list("Ncores"=2,'lr' = 0.75,'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=1)
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    control=list("Ncores"=2,'lr' = 0.75,'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1,"guesses"=10)
    expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
    #
    names <- c("rand","fac","dose")
    Term_n <- c(0,0,1)
    tform <- c("lin","lin","lin")
    keep_constant <- c(0,0,0)
    a_n <- c(-0.1,-0.1,0.2)
    expect_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
})

test_that("Cox Assigned Events, combinations", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
		       "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
		         "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
		      "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
		                  "a"=c(0,   1,   1,   0,   1,   0,   1),
		                  "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
		                  "c"=c(10,  11,  10,  11,  12,  9,   11),
		                  "d"=c(0,   0,   0,   1,   1,   1,   1))
	# For the interval case
	time1 <- "Starting_Age"
	time2 <- "Ending_Age"
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(0.1, 0.1, 0.1, 0.1)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,
	   'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
	   'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
	   'verbose'=T, 'ties'='breslow','double_step'=1)
	keep_constant <- c(1,1,1,1)
	expect_error(RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
	names <- c('a','b','CONST','d')
	keep_constant <- c(0,0,0,0)
	expect_no_error(RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
	
	df$Cancer_Status <- rep(0,nrow(df))
	expect_error(RunCoxEventAssignment(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})

test_that("Poisson Assigned Events, combinations", {
    df <- data.table::data.table("UserID"=c(112, 114, 213, 214, 115, 116, 117),
		       "Starting_Age"=c(18,  20,  18,  19,  21,  20,  18),
		         "Ending_Age"=c(30,  45,  57,  47,  36,  60,  55),
		      "Cancer_Status"=c(0,   0,   1,   0,   1,   0,   0),
		                  "a"=c(0,   1,   1,   0,   1,   0,   1),
		                  "b"=c(1,   1.1, 2.1, 2,   0.1, 1,   0.2),
		                  "c"=c(10,  11,  10,  11,  12,  9,   11),
		                  "d"=c(0,   0,   0,   1,   1,   1,   1))
	# For the interval case
	time1 <- "Starting_Age"
	time2 <- "Ending_Age"
	df$pyr <- df$Ending_Age - df$Starting_Age
	pyr <- 'pyr'
	event <- "Cancer_Status"
	names <- c('a','b','c','d')
	Term_n <- c(0,1,1,2)
	tform <- c("loglin","lin","lin","plin")
	modelform <- "M"
	fir <- 0
	a_n <- c(0.1, 0.1, 0.1, 0.1)

	keep_constant <- c(0,0,0,0)
	der_iden <- 0

	control <- list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,
	   'epsilon' = 1e-3,'dbeta_max' = 0.5, 'deriv_epsilon' = 1e-3,
	   'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,
	   'verbose'=T, 'ties'='breslow','double_step'=1)
	keep_constant <- c(1,1,1,1)
	expect_error(RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
	names <- c('a','b','CONST','d')
	keep_constant <- c(0,0,0,0)
	expect_no_error(RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
	
	df$Cancer_Status <- rep(0,nrow(df))
	expect_error(RunPoissonEventAssignment(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})

test_that("Coxph relative risk combinations", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    keep_constant <- c(1)
    expect_error(Cox_Relative_Risk(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control))
    keep_constant <- c(0)
    names <- c("d","CONST")
    Term_n <- c(0,0)
    tform <- c("loglin", "loglin")
    keep_constant <- c(0,0)
    a_n <- c(-0.1,0.1)
    expect_no_error(Cox_Relative_Risk(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control))
})

test_that("Coxph Martingale combinations", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,2,2,3,3,3)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=e)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=TRUE,"cov_cols"="d","surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="Not_In",'verbose'=TRUE)
    keep_constant <- c(1)
    if (system.file(package='ggplot2')!=""){
        expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
        names <- c("d","CONST")
        Term_n <- c(0,0)
        tform <- c("loglin", "loglin")
        keep_constant <- c(0,0)
        a_n <- c(-0.1,0.1)
        expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
        plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=TRUE,"cov_cols"="Not_In","surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="Not_In",'verbose'=TRUE)
        expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
        plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=TRUE,"cov_cols"="d","surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="Not_In",'verbose'=TRUE)
        df$c <- rep(0,nrow(df))
        expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
    }
})

test_that("Cox_tier_guess combinations", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("a","b")
    Term_n <- c(0,1)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,-15)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 1,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=1,"guesses"=1,"lin_min"=0.001, "lin_max"=1,"loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform", "loglin_method"="uniform",strata=FALSE,term_initial = c(0))
    Strat_Col='a'
	keep_constant <- c(1,1)
    expect_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    names <- c("a","CONST")
    keep_constant <- c(0,0)
    guesses_control=list("Iterations"=1,"guesses"=1,"lin_min"=0.001, "lin_max"=1,"loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform", "loglin_method"="uniform",strata=FALSE,term_initial = c(0),rmin=c(1,1,1,1),rmax=c(1,1))
    expect_no_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})

test_that("Cox_basic_guess_cpp combinations", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("a","b")
    Term_n <- c(0,1)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,-15)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    Strat_Col='a'
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=T)
    model_control=list('strata'=T)
    expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'term_initial' = c(0,1),'verbose'=T)
    names <- c("a","CONST")
    expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    keep_constant <- c(1,1)
    expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    keep_constant <- c(0,0)
    model_control=list('strata'=F)
    lung_temp <- df$lung
    df$lung <- rep(0,nrow(df))
    expect_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    model_control=list('strata'=T)
})

test_that("Gather Guesses list combinations", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c("d","d","d","d")
    Term_n <- c(0,0,0,0)
    tform <- c("loglin",'loglin','loglin','loglin')
    keep_constant <- c(0,0)
    a_n <- list(c(-0.1,6,0.1,0.1))
    a_n_default <- unlist(a_n[1])
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    #
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    guesses_control <- list()
    model_control <- list()
    
    all_names <- unique(names(df))
    #
    dfc <- match(names,all_names)

    term_tot <- max(Term_n)+1
    x_all <- as.matrix(df[,all_names, with = FALSE])
    #
    control <- Def_Control(control)
    guesses_control <- Def_Control_Guess(guesses_control, a_n_default)
    guesses_control$verbose <- TRUE
    model_control <- Def_model_control(model_control)
    #
    names <- c("d","d","d","d","d","d")
    expect_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
    names <- c("d","d","d","d")
    tform <- c("plin",'plin','lin','lin')
    expect_no_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
    tform <- c("Not",'Implemented','Currently','Error')
    expect_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
    keep_constant <- c(0,1)
    names <- c("d","d","d","d")
    Term_n <- c(0,0,0,0)
    tform <- c("loglin",'loglin','loglin','loglin')
    expect_no_error(Gather_Guesses_CPP(df, dfc, names, Term_n, tform, keep_constant, a_n, x_all, a_n_default, modelform, fir, control, guesses_control))
})

test_that("Default control guess combinations", {
    control_def=list("verbose"=T)
    a_n <- c(1,2,3)
    expect_no_error(Def_Control_Guess(control_def,a_n))
    control_def=list("verbose"="p")
    expect_error(Def_Control_Guess(control_def,a_n))
    control_def=list("verbose"=T,"guess_constant"=c(1))
    expect_no_error(Def_Control_Guess(control_def,a_n))
})

test_that("linked formula combinations", {
    tforms <- list("first"="quad")
    paras  <- list("first"=c(0.1,10))
    expect_error(Linked_Dose_Formula(tforms,paras,verbose='p'))
    paras  <- list("first"=c(0.1,"10"))
    expect_error(Linked_Dose_Formula(tforms,paras,T))
    #
    tforms <- list("first"="exp")
    paras  <- list("first"=c(0.1,"10",5))
    expect_error(Linked_Dose_Formula(tforms,paras,TRUE))
    paras  <- list("first"=c(0.1,10,"5"))
    expect_error(Linked_Dose_Formula(tforms,paras,TRUE))
})

test_that("Iteract formula operation error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?++?b","a?*?b")
    new_names <- c("","")
    expect_error(interact_them(df,interactions,new_names,FALSE))
})

test_that("gmix omnibus use", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    set.seed(3742)
    df$rand <- floor(runif(nrow(df), min=0, max=5))

    time1 <- "t0"
    time2 <- "t1"
    df$censor <- (df$lung==0)
    #
    event <- "lung"
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)

    verbose <- FALSE

    model_list <- c('GMIX-R','GMIX-E','GMIX')
    names <- c("dose","fac","dose","fac","rand")
    Term_n <- c(0,0,1,1,2)
    tform <- c("loglin","loglin","plin","plin","loglin")
    keep_constant <- c(0,0,0,0,0)
    a_n <- c(-0.1,0.1,0.2,0.3,-0.5)
    df_order <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant, "a_n"=a_n, "names"=names, "order"=1:5)

    count <- 0
    der_iden <- 0
    cens_weight <- c(0)
    for (model_i in 1:3){
        modelform <- model_list[model_i]
        if (modelform=='GMIX'){
            for (fir in c(0,1,2)){
                for (term_i in 0:3){
                    model_control=list('strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'CR'=FALSE)
                    if (fir==0){
                        model_control$gmix_term <- c(0,term_i%%2, floor(term_i/2))
                    } else if (fir==1){
                        model_control$gmix_term <- c(term_i%%2,0, floor(term_i/2))
                    }  else if (fir==2){
                        model_control$gmix_term <- c(term_i%%2, floor(term_i/2),1)
                    }
                    #
                    df_order$order <- sample(df_order$order)
                    setorderv(df_order, c("order"))
                    Term_n <- df_order$Term_n
                    tform <- df_order$tform
                    keep_constant <- df_order$keep_constant
                    a_n <- df_order$a_n
                    names <- df_order$names
                    #
                    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
                    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="rand", model_control=model_control, cens_weight=cens_weight))
                }
            }
        } else {
            for (fir in c(0,1,2)){
                model_control=list('strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'CR'=FALSE)
                #
                df_order$order <- sample(df_order$order)
                setorderv(df_order, c("order"))
                Term_n <- df_order$Term_n
                tform <- df_order$tform
                keep_constant <- df_order$keep_constant
                a_n <- df_order$a_n
                names <- df_order$names
                #
                control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
                expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="rand", model_control=model_control, cens_weight=cens_weight))
            }
        }
    }
})