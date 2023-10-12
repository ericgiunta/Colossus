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
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("Ncores"=2,'lr' = 0.95,'maxiter' = -1,'halfmax' = 1,'epsilon' = 1e-9,'dbeta_max' = 1.0,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=1.0,'verbose'=FALSE, 'double_step'=1)
    expect_error(RunPoissonRegression(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    pyr <- "a"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("Ncores"=2,'lr' = 0.95,'maxiter' = -1,'halfmax' = 1,'epsilon' = 1e-9,'dbeta_max' = 1.0,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=1.0,'verbose'=FALSE, 'double_step'=1)
    expect_error(RunPoissonRegression(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson no events", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    pyr <- "a"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("Ncores"=2,'lr' = 0.95,'maxiter' = -1,'halfmax' = 1,'epsilon' = 1e-9,'dbeta_max' = 1.0,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=1.0,'verbose'=FALSE, 'double_step'=1)
    expect_error(RunPoissonRegression(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
test_that("Poisson no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,1,0,0,0,1,0)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    time1 <- "a"
    time2 <- "b"
    event <- "c"
    pyr <- "a"
    names <- c("d")
    Term_n <- c(0)
    tform <- c("loglin")
    keep_constant <- c(0)
    a_n <- c(-0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("Ncores"=2,'lr' = 0.95,'maxiter' = -1,'halfmax' = 1,'epsilon' = 1e-9,'dbeta_max' = 1.0,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=1.0,'verbose'=FALSE, 'double_step'=1)
    expect_no_error(RunPoissonRegression(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control))
})
