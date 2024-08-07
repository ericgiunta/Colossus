test_that( "Coxph time column missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
    time1 <- "a_bad"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that( "Coxph no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
# test_that( "Coxph no error", {
#     a <- c(0,1,2,3,4,5,6)
#     b <- c(1,2,3,4,5,6,7)
#     c <- c(0,1,0,0,0,1,0)
#     d <- c(3,4,5,6,7,8,9)
#     df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
#     time1 <- "a"
#     time2 <- "b"
#     event <- "c"
#     names <- c( "d" )
#     term_n <- c(0)
#     tform <- c( "loglin" )
#     keep_constant <- c(0)
#     a_n <- c(-0.1)
#     modelform <- "M"
#     fir <- 0
#     der_iden <- 0
#     control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
#     expect_no_error(RunCoxRegression(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
# })

test_that( "Coxph_strata time column missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,0,0,1,0,1)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d, "e"=e)
    time1 <- "a_bad"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    strat_col <- "e"
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression_STRATA(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col))
})
test_that( "Coxph_strata no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,0,0,1,0,1)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d, "e"=e)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    strat_col <- "e"
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression_STRATA(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col))
})
test_that( "Coxph_strata no strata", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,0,0,1,0,1)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d, "e"=e)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    strat_col <- "e_bad"
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression_STRATA(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col))
})
test_that( "Coxph_strata strata with no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    e <- c(1,1,0,0,1,0,1)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d, "e"=e)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    names <- c( "d" )
    strat_col <- "e_bad"
    term_n <- c(0)
    tform <- c( "loglin" )
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxRegression_STRATA(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col))
})
# test_that( "Coxph_strata no error", {
#     a <- c(0,1,2,3,4,5,6)
#     b <- c(1,2,3,4,5,6,7)
#     c <- c(0,1,0,0,0,1,0)
#     d <- c(3,4,5,6,7,8,9)
#     e <- c(1,1,0,0,1,0,1)
#     df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d, "e"=e)
#     time1 <- "a"
#     time2 <- "b"
#     event <- "c"
#     names <- c( "d" )
#     strat_col <- "e"
#     term_n <- c(0)
#     tform <- c( "loglin" )
#     keep_constant <- c(0)
#     a_n <- c(-0.1)
#     modelform <- "M"
#     fir <- 0
#     der_iden <- 0
#     control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
#     expect_no_error(RunCoxRegression_STRATA(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,strat_col))
# })

# test_that( "Coxph relative risk no error", {
#     a <- c(0,1,2,3,4,5,6)
#     b <- c(1,2,3,4,5,6,7)
#     c <- c(0,1,0,0,0,1,0)
#     d <- c(3,4,5,6,7,8,9)
#     df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
#     time1 <- "a"
#     time2 <- "b"
#     event <- "c"
#     names <- c( "d" )
#     term_n <- c(0)
#     tform <- c( "loglin" )
#     keep_constant <- c(0)
#     a_n <- c(-0.1)
#     modelform <- "M"
#     fir <- 0
#     control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
#     expect_no_error(Cox_Relative_Risk(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, fir, control))
# })

test_that( "Coxph null time column missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
    time1 <- "a_bad"
    time2 <- "b"
    event <- "c"
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxNull(df, time1, time2, event, control))
})
test_that( "Coxph null no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
    expect_error(RunCoxNull(df, time1, time2, event, control))
})
# test_that( "Coxph null no error", {
#     a <- c(0,1,2,3,4,5,6)
#     b <- c(1,2,3,4,5,6,7)
#     c <- c(0,1,0,0,0,1,0)
#     d <- c(3,4,5,6,7,8,9)
#     df <- data.table( "a"=a, "b"=b, "c"=c, "d"=d)
#     time1 <- "a"
#     time2 <- "b"
#     event <- "c"
#     control <- list( "ncores"=2, 'lr' = 0.75, 'maxiter' = -1, 'halfmax' = 5, 'epsilon' = 1e-9,  'deriv_epsilon' = 1e-9, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1)
#     expect_no_error(RunCoxNull(df, time1, time2, event, control))
# })

