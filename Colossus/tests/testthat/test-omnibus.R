test_that("Coxph strata_basic_single_null", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
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
    verbose <- FALSE

    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 1,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    for (i in c(TRUE,FALSE)){
        for (j in c(TRUE,FALSE)){
            for (k in c(TRUE,FALSE)){
                for (l in c(TRUE,FALSE)){
                    model_control=list('strata'=i, 'basic'=j, 'single'=k, 'null'=l)
                    if (verbose){print(model_control)}
                    a_n <- c(0.01)
                    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=c("loglin"), keep_constant=keep_constant, a_n=a_n, modelform="M", fir=0, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control))
                    if (verbose){print("---------------")}
                }
            }
        }
    }
})
test_that("Coxph basic_single_null match", {
    fname <- 'll_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    time1 <- "t0"
    time2 <- "t1"
    event <- "lung"
    names <- c("dose")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(0.0)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    verbose <- FALSE

    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(-1,-1),'halfmax' = -1,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    model_control=list('strata'=FALSE, 'basic'=FALSE, 'single'=FALSE, 'null'=FALSE)
    e0 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=c("loglin"), keep_constant=keep_constant, a_n=a_n, modelform="M", fir=0, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control)
    for (j in c(TRUE,FALSE)){
        for (k in c(TRUE,FALSE)){
            for (l in c(TRUE,FALSE)){
                model_control=list('strata'=FALSE, 'basic'=j, 'single'=k, 'null'=l)
                if (verbose){print(model_control)}
                a_n <- c(0.0)
                e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=c("loglin"), keep_constant=keep_constant, a_n=a_n, modelform="M", fir=0, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control)
                expect_equal(e0$LogLik,e1$LogLik,tolerance=1e-2)
                if (verbose){print("---------------")}
            }
        }
    }
})
test_that("Coxph strata_basic_single_CR", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)

    df$rand <- floor(runif(nrow(df), min=0, max=5))

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
    control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
    plot_options <- list("name"="run_04","verbose"=FALSE,"studyID"="studyID","age_unit"="years")
    dft <- GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)
    #
    file.remove('weight_surv_plot_run_04.jpeg')
    #
    t_ref <- dft$t
    surv_ref <- dft$surv
    t_c <- df$t1
    cens_weight <- approx(t_ref, surv_ref, t_c,rule=2)$y
    #
    #    cens_weight <- dft$surv
    event <- "lung"
    a_n <- c(-0.1,-0.1)
    keep_constant <- c(0,0)

    #stop()
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)

    verbose <- FALSE

    for (i in c(TRUE,FALSE)){
        for (j in c(TRUE,FALSE)){
            for (k in c(TRUE,FALSE)){
                for (l in c(TRUE,FALSE)){
                    model_control=list('strata'=i, 'basic'=j, 'single'=k, 'CR'=l)
                    if (verbose){print(model_control)}
                    a_n <- c(-0.1,-0.1)
                    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
                    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="rand", model_control=model_control, cens_weight=cens_weight))
                    a_n <- c(-0.1,-0.1)
                    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='efron','double_step'=0)
                    expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="rand", model_control=model_control, cens_weight=cens_weight))
                    if (verbose){print("---------------")}
                }
            }
        }
    }
})

test_that("Pois strata_single", {
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
    
    for (i in c(TRUE,FALSE)){
        for (j in c(TRUE,FALSE)){
            model_control=list('strata'=i, 'single'=j)
            if (verbose){print(model_control)}
            a_n <- c(0.01,0.1,0.1)
            modelform <- "PAE"
            expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
            a_n <- c(0.01,0.1,0.1)
            modelform <- "A"
            expect_no_error(RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col,model_control))
            if (verbose){print("---------------")}
        }
    }
})
