library(Colossus)
library(data.table)
library(parallel)

if (identical(Sys.getenv("TESTTHAT"), "true")){
} else {
    set.seed(12132024)
#    subjID <- c()
#    t0 <- c()
#    t1 <- c()
#    dose <- c()
#    sex <- c()
#    p <- c()
#    event <- c()
#
#    param_back <- c(1/700000, -0.03)
#    param_front <- c(-0.95, 0.1)
#
#    for (i in 1:10000){
#        t <- 18
#        d <- 0
#        s <- as.integer(runif(1) * 2)
#        while (t < 99){
#            #
#            t0 <- c(t0, t)
#            d_int <- runif(1)/60
#            dose <- c(dose, d_int + d)
#            d <- d_int + d
#            sex <- c(sex, s)
#            subjID <- c(subjID, 1000000 + i)
#            #
#            cov_rate <- (1 + param_front[1] * d) * exp(s * param_front[2])
#            back_rate <- param_back[1] * t^2 * exp(param_back[2] * t)
#            #
#            p_event <- cov_rate * back_rate
#            p <- c(p, p_event)
#            U <- runif(1)
#            #
#            if (U < p_event){
#                t1 <- c(t1, t + U/p_event)
#                event <- c(event, 1)
#                t <- 100
#            } else {
#                t1 <- c(t1, t + 1)
#                event <- c(event, 0)
#                t <- t + 1
#            }
#        }
#    }
#
#    df <- data.table(
#        "SubjectID"=subjID,
#        "t0"=t0,
#        "t1"=t1,
#        "dose"=dose,
#        "sex"=sex,
#        "event"=event
#    )
#    fwrite(df, "excess_test.csv")
#    stop()
    df <- fread("excess_test/excess_test.csv")
    df <- df[event < 2]
    print(sum(df$event))
    print(nrow(df))
#    stop()
    time1 <- "age_entry"
    time2 <- "exit_adjust"
    event <- "event"
    names <- c("dose", "sex")
    term_n <- c(0, 0)
    tform <- c("loglin", "loglin")
    keep_constant <- c(0, 0)
    a_n <- c(0.1, 0.1)
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("lr" = 0.75, "maxiters" = c(100, 100), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 3, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control = control)
    Interpret_Output(e)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 0, "manual" = FALSE, "maxstep" = 50, "search_mult" = 2.0)
    e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 1, "manual" = FALSE, "maxstep" = 50, "search_mult" = 2.0)
    e2 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e2)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 0, "manual" = TRUE, "maxstep" = 50, "search_mult" = 2.0)
    e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 1, "manual" = TRUE, "maxstep" = 50, "search_mult" = 2.0)
    e2 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e2)

    time1 <- "age_entry"
    time2 <- "exit_adjust"
    event <- "event"
    names <- c("dose", "sex")
    term_n <- c(0, 0)
    tform <- c("plin", "loglin")
    keep_constant <- c(0, 0)
    a_n <- list(c(0.1, -0.1),c(0.3, -0.1),c(0.5, -0.1),c(0.7, -0.1),c(0.9, -0.1),c(1.1, -0.1))
    modelform <- "M"
    fir <- 0
    der_iden <- 0
    control <- list("lr" = 0.75, "maxiters" = c(10, 10,10,10,10,10,100), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 3, "ties" = "breslow", "double_step" = 1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, a_n, modelform, control = control)
    Interpret_Output(e)
    control <- list("lr" = 0.75, "maxiters" = c(100), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 3, "ties" = "breslow", "double_step" = 1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 0, "manual" = FALSE, "maxstep" = 50, "search_mult" = 1.5)
    e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 1, "manual" = FALSE, "maxstep" = 50, "search_mult" = 2.0)
    e2 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e2)
    control <- list("lr" = 0.75, "maxiters" = c(100), "halfmax" = 5, "epsilon" = 1e-6, "deriv_epsilon" = 1e-6, "abs_max" = 1.0, "change_all" = TRUE, "dose_abs_max" = 100.0, "verbose" = 3, "ties" = "breslow", "double_step" = 1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 0, "manual" = TRUE, "maxstep" = 50, "search_mult" = 1.5)
    e1 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e1)
    model_control <- list("basic" = FALSE, "log_bound" = TRUE, "alpha" = 0.05, "para_number" = 1, "manual" = TRUE, "maxstep" = 50, "search_mult" = 2.0)
    e2 <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n, tform, keep_constant, e$beta_0, modelform, control = control, model_control = model_control)
    Interpret_Output(e2)
    warnings()
}
