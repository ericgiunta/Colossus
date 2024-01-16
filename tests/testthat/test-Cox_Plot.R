
######################
## PLOTTING CHECKS
######################

test_that("Coxph plot time column missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a_bad"
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
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=FALSE,"surv_curv"=FALSE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=FALSE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot no type", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("Martingale"=FALSE,"surv_curv"=FALSE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot strata col not in data", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=FALSE,"strat_haz"=TRUE,"Strat_Col"="e", "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot strata col not given", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=FALSE,"strat_haz"=TRUE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot martingale dose col not in data", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=TRUE,"surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot martingale dose col not given", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=TRUE,"surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="a",'verbose'=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot Survival ID col not in data", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="f",'verbose'=TRUE,"KM"=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot Survival ID col not given", {
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV","run"),"Martingale"=FALSE,"surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE,'verbose'=TRUE,"KM"=TRUE)
    expect_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})

test_that("Coxph plot no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,0,1,0,1,0,0)
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")), "surv_curv"=TRUE,"studyID"="a",'verbose'=TRUE)
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph plot stratafied no error", {
    a <- c(0,0,0,0,0,0,1)
	b <- c(2,3,4,2,3,2,4)
	c <- c(1,0,1,0,1,1,0)
	d <- c(3,4,5,6,7,8,9)
	e <- c(1,1,0,0,1,0,1)
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
	plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")), "surv_curv"=TRUE,"strat_haz"=TRUE,"Strat_Col"="e","studyID"="a",'verbose'=TRUE)

    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})

test_that("Coxph risk no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,0,1,0,1,0,0)
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("RISK",paste(tempfile(),"run",sep="")),"studyID"="a",'verbose'=TRUE)
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})

test_that("Coxph schoenfeld no error", {
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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SCHOENFELD",paste(tempfile(),"run",sep="")),"studyID"="t0",'verbose'=TRUE)
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
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
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")),"Martingale"=TRUE,"cov_cols"="d","surv_curv"=FALSE,"strat_haz"=FALSE, "smooth_haz"=FALSE, "studyID"="e",'verbose'=TRUE)
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
test_that("Coxph KM no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,0,1,0,1,0,0)
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
    der_iden <- 0
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options=list("type"=c("SURV",paste(tempfile(),"run",sep="")), "surv_curv"=TRUE,"studyID"="a",'verbose'=TRUE,"KM"=TRUE)
    expect_no_error(RunCoxPlots(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options))
})
