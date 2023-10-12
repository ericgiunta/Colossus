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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=T, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=T)
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    model_control <- list('strata'=T)
    expect_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col,model_control))
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'term_initial' = c(0,1),'verbose'=T)
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
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
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








