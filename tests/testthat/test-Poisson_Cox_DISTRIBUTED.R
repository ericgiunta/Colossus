test_that("Pois_tier_guess", {
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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=TRUE,'term_initial' = c(0,1),'verbose'=FALSE)
    Strat_Col=c('b')

    expect_no_error(RunPoissonRegression_Tier_Guesses(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control, guesses_control, Strat_Col))
})

test_that("Cox_basic_guess_cpp", {
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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=FALSE)
    expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    a_n <- list(c(0.01,-15),c(0.1,-5))
    expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})
test_that("Cox_tier_guess", {
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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=10,"guesses"=10,"lin_min"=0.001, "lin_max"=1,"loglin_min"=-1,"loglin_max"=1, "lin_method"="uniform", "loglin_method"="uniform",strata=FALSE,term_initial = c(0))
    Strat_Col='a'

    expect_no_error(RunCoxRegression_Tier_Guesses(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})

test_that("Cox_basic_guess_cpp_strata", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("b")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0.01)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=TRUE,'term_initial' = c(0,1),'verbose'=FALSE)
    Strat_Col='a'
    expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    a_n <- list(c(0.01),c(0.1))
    expect_no_error(RunCoxRegression_Guesses_CPP(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})

test_that("Poisson_basic_guess_cpp", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    pyr <- "t1"
    event <- "lung"
    names <- c("a","b")
    Term_n <- c(0,1)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,-15)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=FALSE)
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    a_n <- list(c(0.01,-15),c(0.1,-5))
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})
test_that("Poisson_basic_guess_cpp_strata", {
    fname <- 'MULTI_COV.csv'
    colTypes=c("double","double","integer","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    pyr <- "t1"
    event <- "lung"
    names <- c("b")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0.01)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=TRUE,'term_initial' = c(0,1),'verbose'=FALSE)
    Strat_Col='a'
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
    a_n <- list(c(0.01),c(0.1))
    expect_no_error(RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col))
})
