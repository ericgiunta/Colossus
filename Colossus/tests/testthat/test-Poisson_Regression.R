test_that("Poisson time column missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a_bad"
    time2 <- "b"
    event <- "c"
    pyr <- "a_bad"
    names <- c("d")
    Term_n <- 1
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 1,'epsilon' = 0,'dbeta_max' = 0.5,'deriv_epsilon' = 0, 'abs_max'=1.0,'change_all'=FALSE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='efron')
    expect_error(RunPoissonRegression(df, time1, time2, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a_bad"
    time2 <- "b"
    event <- "c"
    pyr <- "a"
    names <- c("d")
    Term_n <- 1
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 1,'epsilon' = 0,'dbeta_max' = 0.5,'deriv_epsilon' = 0, 'abs_max'=1.0,'change_all'=FALSE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='efron')
    expect_error(RunPoissonRegression(df, time1, time2, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
