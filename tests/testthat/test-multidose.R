test_that( "Coxph multidose", {
    fname <- 'll_comp_0.csv'
    colTypes <- c( "double", "double", "double", "integer", "integer" )
    df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
    set.seed(3742)
    df$rand <- floor(runif(nrow(df), min=0, max=5))
    df$rand0 <- floor(runif(nrow(df), min=0, max=5))
    df$rand1 <- floor(runif(nrow(df), min=0, max=5))
    df$rand2 <- floor(runif(nrow(df), min=0, max=5))
    print( "starting" )
    time1 <- "t0"
    time2 <- "t1"
    # df$censor <- (df$lung==0)
    df$lung <- (df$lung>0)
    # event <- "censor"
    names <- c( "dose", "rand" )
    term_n <- c(0,0)
    tform <- c( "loglin", "loglin" )
    realization_columns = matrix(c( "rand0", "rand1", "rand2" ),nrow=1)
    realization_index=c( "rand" )
    keep_constant <- c(1,0)
    a_n <- c(0,0)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    cens_weight <- c(0)
    #
    event <- "lung"
    a_n <- c(-0.1,-0.1)
    keep_constant <- c(0,0)

    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = 10, 'halfmax' = 2, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=3, 'ties'='breslow', 'double_step'=1)

    verbose <- FALSE
    j_iterate <- 1
    # LL_comp <- c(-69.51585, -69.51585, -77.97632, -77.97632, -59.95167, -60.05273, -75.34028, -75.3691, -69.51585, -69.51585, -77.97632, -77.97632, -59.95167, -60.05273, -75.34028, -75.3691, -111.3009, -111.3009, -119.9814, -119.9814, -100.8329, -101.007, -117.0147, -117.0539, -111.3009, -111.3009, -119.9814, -119.9814, -100.8329, -101.007, -117.0147, -117.0539)
    for (i in c(FALSE,TRUE)){
        for (j in c(FALSE,TRUE)){
            model_control <- list( 'strata'=i, 'basic'=j)
            if (verbose){print(model_control)}
            a_n <- c(-0.1,-0.1)
            #expect_equal(0,0)
            control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = 10, 'halfmax' = 2, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=3, 'ties'='breslow', 'double_step'=1)
            expect_no_error(RunCoxRegression_Omnibus_Multidose(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, realization_columns=realization_columns, realization_index=realization_index, control=control,strat_col="fac", model_control=model_control, cens_weight=cens_weight))
        }
    }
})
