test_that("Coxph loglin_M", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","fac")
    Term_n <- c(0,0)
    tform <- c("loglin","loglin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$beta_0,c(-0.0996,-0.05697),tolerance=1e-2)
})
test_that("Coxph loglin_plin_M", {
    fname <- 'l_pl_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","fac")
    Term_n <- c(0,0)
    tform <- c("loglin","plin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,0.5)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$beta_0,c(0.1747772,0.75),tolerance=1e-2)
})
test_that("Coxph loglin_plin_A", {
    fname <- 'l_pl_A_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","fac")
    Term_n <- c(0,1)
    tform <- c("loglin","plin")
    keep_constant <- c(0,0)
    a_n <- c(0.01,0.5)
    modelform <- "A"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 100,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$beta_0,c(0.11, 1.01),tolerance=1e-2)
})

test_that("Coxph dose list", {
    fname <- 'dose.csv'
    colTypes=c("double","double","double","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","dose","dose","dose","dose","dose","dose","dose","dose","dose","dose")
    Term_n <- c(0,0,0,0,0,0,0,0,0,0,0)
    tform <- c("loglin_top","lin_slope","lin_int","quad_slope","step_slope","step_int","lin_quad_slope","lin_quad_int","lin_exp_slope","lin_exp_int","lin_exp_exp_slope")
    keep_constant <- c(0,0,0,0,0,0,0,0,0,0,0)
    tform <- c("loglin_top","lin_slope","lin_int","quad_slope","step_slope","step_int","lin_quad_slope","lin_quad_int","lin_exp_slope","lin_exp_int","lin_exp_exp_slope")
    a_n <-   c(-0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    expect_equal(e$beta_0,c(-0.10032969, -0.09962793,  0.99931507, -0.10823105,  1.04233436,  2.01958463, 0.32083364, 1.74980209, 0.20085333, 0.88407633, 1.00020558),tolerance=1e-2)
})

test_that("Coxph fixed intercept", {
    fname <- 'dose.csv'
    colTypes=c("double","double","double","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose","dose","dose","dose","dose")
    Term_n <- c(0,0,0,0,0)
    tform <- c("loglin","lin_slope","lin_int","step_slope","step_int")
    keep_constant <- c(0,0,1,0,1)
    a_n <-   c(-0.1          ,0.1       ,-1,0.1,-1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    expect_no_error(RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})





