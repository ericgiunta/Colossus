#test_that("Coxph loglin_M Basic", {
#    fname <- 'll_0.csv'
#    colTypes=c("double","double","double","integer","integer")
#    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
#    time1 <- "t0"
#    time2 <- "t1"
#    event <- "lung"
#    names <- c("dose","fac")
#    Term_n <- c(0,0)
#    tform <- c("loglin","loglin")
#    keep_constant <- c(0,0)
#    a_n <- c(0.01,0.1)
#    modelform <- "M"
#    fir <- 0
#    der_iden <- 0
#    control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#    e <- RunCoxRegression_Basic(df, time1, time2, event, names, keep_constant, a_n, der_iden, control)
#    expect_equal(e$beta_0,c(-0.0996,-0.05697),tolerance=1e-2)
#})
test_that("Coxph loglin_M Strata", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0.01)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression_STRATA(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,"fac")
    expect_equal(e$beta_0,c(-0.106),tolerance=1e-2)
})
test_that("Coxph loglin_M Single", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
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
    control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression_Single(df, time1, time2, event, names, Term_n, tform, a_n, modelform, fir, control)
    expect_equal(e$AIC,2109,tolerance=1e-2)
})
test_that("Coxph loglin_M Null", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=detectCores()-1,data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
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
    control=list('lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxNull(df, time1, time2, event, control)
    expect_equal(e$AIC,1052,tolerance=1e-2)
})




