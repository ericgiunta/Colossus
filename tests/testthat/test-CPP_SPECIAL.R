test_that("Coxph loglin_M Basic", {
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
    e <- RunCoxRegression_Basic(df, time1, time2, event, names, keep_constant, a_n, der_iden, control)
    expect_equal(e$beta_0,c(-0.09962, -0.05698),tolerance=1e-2)
})

test_that("Coxph censoring weight", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    df$censor <- (df$lung==0)
    event <- "censor"
    names <- c("dose","fac")
    Term_n <- c(0,0)
    tform <- c("loglin","loglin")
    keep_constant <- c(1,0)
    a_n <- c(0,0)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = -1,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    plot_options <- list("name"=paste(tempfile(),"run",sep=""),"verbose"=TRUE,"studyID"="studyID","age_unit"="years")
    dft <- GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)
    #
    #
    t_ref <- dft$t
    surv_ref <- dft$surv
    t_c <- df$t1
    cens_weight <- approx(t_ref, surv_ref, t_c,rule=2)$y
    message(sum(cens_weight))
    #
    event <- "lung"
    a_n <- c(-0.1,-0.1)
    keep_constant <- c(0,0)
    e0 <- RunCoxRegression_CR(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,cens_weight)

    a_n <- c(-0.1,-0.1)
    keep_constant <- c(0,0)
    e1 <- RunCoxRegression(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control)
    #
    expect_equal(e0$LogLik - e1$LogLik,8.680511,tolerance=1e-2)
})


