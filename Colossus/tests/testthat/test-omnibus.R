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
    plot_options <- list("name"=paste(tempfile(),"run",sep=""),"verbose"=FALSE,"studyID"="studyID","age_unit"="years")
    dft <- GetCensWeight(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, control, plot_options)
    #
    #
    t_ref <- dft$t
    surv_ref <- dft$surv
    t_c <- df$t1
    cens_weight <- approx(t_ref, surv_ref, t_c,rule=2)$y
    #
    event <- "lung"
    a_n <- c(-0.1,-0.1)
    keep_constant <- c(0,0)

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

test_that("risk check omnibus plain", {
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
    names <- c("dose","fac","dose","fac","rand")
    Term_n <- c(0,0,1,1,1)
    tform <- c("loglin","lin","lin","plin","loglin")
    keep_constant <- c(0,0,0,0,0)
    a_n <- c(-0.1,0.1,0.2,0.3,-0.5)
    modelform <- "M"
    fir <- 0
    der_iden <- 0

    cens_weight <- c(0)
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)

    verbose <- FALSE

    df_order <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant, "a_n"=a_n, "names"=names, "order"=1:5)

    model_list <- c('M','A','PA','PAE')
    means <- c(0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.223545565747452, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0079328174864236, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.0532705372693687, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931, 0.186140663450931)
    medians <- c(0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0.155704743824735, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.0822977918796923, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932, 0.112941187026932)
    sums <- c(53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 223.545565747452, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 7.9328174864236, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 53.2705372693687, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931, 186.140663450931)

    for (model_i in 1:4){
        modelform <- model_list[model_i]
        for (fir in c(0,1)){
            for(i in 1:5){
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
                e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, control=control, model_control=model_control)$Risk

                a_i <- (model_i-1)*10 + fir*5 + i
                expect_equal(mean(e),means[a_i],tolerance=1e-2)
                expect_equal(median(e),medians[a_i],tolerance=1e-2)
                expect_equal(sum(e),sums[a_i],tolerance=1e-2)
            }
        }
    }
})

test_that("risk check omnibus gmix", {
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

    means <- c(0.6918334,0.6918334,0.6918334,2.9080098,3.8839872,2.7077608,0.6918334,1.2728663,1.8214142,2.9080098,0.6918334,1.4807626,1.8214142,3.8839872,0.6918334,1.4807626,1.2728663,2.7077608)
    medians <- c(0.5871897,0.5871897,0.5871897,2.8144283,3.7521540,2.4285966,0.5871897,1.0226682,1.7472830,2.8144283,0.5871897,1.2200643,1.7472830,3.7521540,0.5871897,1.2200643,1.0226682,2.4285966)
    sums <- c(691.8334,691.8334,691.8334,2908.0098,3883.9872,2707.7608,691.8334,1272.8663,1821.4142,2908.0098,691.8334,1480.7626,1821.4142,3883.9872,691.8334,1480.7626,1272.8663,2707.7608)


    count <- 0
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
                    e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, control=control, model_control=model_control)$Risk
                    count <- count + 1
                    expect_equal(mean(e),means[count],tolerance=1e-2)
                    expect_equal(median(e),medians[count],tolerance=1e-2)
                    expect_equal(sum(e),sums[count],tolerance=1e-2)
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
                e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, control=control, model_control=model_control)$Risk
                count <- count + 1
                expect_equal(mean(e),means[count],tolerance=1e-2)
                expect_equal(median(e),medians[count],tolerance=1e-2)
                expect_equal(sum(e),sums[count],tolerance=1e-2)
            }
        }
    }
})

test_that("risk check omnibus dose", {
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
    names <- c("dose","fac","dose","fac","rand")
    Term_n <- c(0,0,1,1,1)
    tform <- c("loglin","lin","lin","plin","loglin")
    keep_constant <- c(0,0,0,0,0)
    a_n <- c(-0.1,0.1,0.2,0.3,-0.5)
    modelform <- "M"
    fir <- 0
    der_iden <- 0

    cens_weight <- c(0)
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)

    verbose <- FALSE

    df_order <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant, "a_n"=a_n, "names"=names, "order"=1:5)

    model_list <- c('M','A','PA','PAE')
    names <- c("dose","dose","dose","dose","dose","dose","dose","dose","dose","dose","dose",  "fac","fac","fac","fac","fac","fac","fac","fac","fac","fac","fac")
    Term_n <- c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1)
    tform <- c("loglin_top","lin_slope","lin_int","quad_slope","step_slope","step_int","lin_quad_slope","lin_quad_int","lin_exp_slope","lin_exp_int","lin_exp_exp_slope","loglin_top","lin_slope","lin_int","quad_slope","step_slope","step_int","lin_quad_slope","lin_quad_int","lin_exp_slope","lin_exp_int","lin_exp_exp_slope")
    keep_constant <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    a_n <-   c(-0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1, -0.1          ,-0.1       ,1        ,-0.1        ,1           ,2         ,0.3             ,1.5           ,0.2            ,0.7          ,1)

    df_order <- data.table("Term_n"=Term_n, "tform"=tform, "keep_constant"=keep_constant, "a_n"=a_n, "names"=names, "order"=1:22)
    means <- c(3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 2.67297456164882, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 1.75338940554459, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 3.26883019325272, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529, 2.91092317948529)
    medians <- c(2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 2.55867858245872, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 1.6364085377938, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.89841334635291, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858, 2.90408823467858)
    sums <- c(3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 2672.97456164882, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 1753.38940554459, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 3268.83019325272, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529, 2910.92317948529)

    for (model_i in 1:4){
        modelform <- model_list[model_i]
        for (fir in c(0,1)){
            for(i in 1:22){
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
                e <- Cox_Relative_Risk(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, control=control, model_control=model_control)$Risk
                a_i <- (model_i-1)*44 + fir*22 + i
                expect_equal(mean(e),means[a_i],tolerance=1e-2)
                expect_equal(median(e),medians[a_i],tolerance=1e-2)
                expect_equal(sum(e),sums[a_i],tolerance=1e-2)
            }
        }
    }
})

test_that("check deviation calc", {
    fname <- 'll_comp_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    set.seed(3742)
    df$rand <- floor(runif(nrow(df), min=0, max=5))

    time1 <- "t0"
    time2 <- "t1"


    #
    event <- "lung"
    names <- c("dose","fac","rand")
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
    for (i in 1:3){
        a_n <- c(0.6465390, 0.4260961, 0.1572781)
        keep_constant <- c(0,0,0)
        keep_constant[i] <- 1
        #
        control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
        e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control)
        devs <- c(devs, sum(e$Standard_Deviation))
    }
    a_n <- c(0.6465390, 0.4260961, 0.1572781)
    keep_constant <- c(0,0,0)
    #
    control=list("Ncores"=2,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="fac", model_control=model_control)
    devs <- c(devs, sum(e$Standard_Deviation))
    expect_equal(devs,c(0.6091269, 0.5356671, 0.7385757, 0.9448081),tolerance=1e-4)
})

test_that("check Linear Constraints", {
    fname <- 'l_pl_0.csv'
    colTypes=c("double","double","double","integer","integer")
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    df$pyr <- df$t1 - df$t0
    time1 <- "t0"
    time2 <- "t1"
    pyr <- 'pyr'
    event <- "lung"
    names <- c("dose","fac")
    Term_n <- c(0,0)
    tform <- c("loglin","plin")
    keep_constant <- c(0,0)
    model_control=list('strata'=F, 'basic'=F, 'single'=F, 'null'=F,'constraint'=T)
    Constraint_Matrix <- matrix(c(1,-1),nrow=1)
    Constraint_const  <- c(0.0)
    for (i in 1:20){
        a_n <- 2*runif(2)-1
        del <- abs(a_n[1]-a_n[2])
        a_n0 <- rep(sum(a_n)/2,2)
        a_n <- a_n0 - c(-del/2,del/2)
        modelform <- "M"
        fir <- 0
        der_iden <- 0
        control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
        e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col="fac", model_control=model_control,Cons_Mat=Constraint_Matrix, Cons_Vec=Constraint_const)
        expect_equal(e$beta_0,c(0.357333, 0.357333),tolerance=1e-2)
    }
    for (i in 1:20){
        a_n <- 2*runif(2)-1
        del <- abs(a_n[1]-a_n[2])
        a_n0 <- rep(sum(a_n)/2,2)
        a_n <- a_n0 + c(-del/2,del/2)
        modelform <- "M"
        fir <- 0
        der_iden <- 0
        control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
        e <- RunPoissonRegression_Omnibus(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,Strat_Col="fac",model_control=model_control,Cons_Mat=Constraint_Matrix, Cons_Vec=Constraint_const)
        expect_equal(e$beta_0,c(-0.472812, -0.472812),tolerance=1e-2)
    }
})

