#install.packages(c("fuzzr","Colossus","data.table"))

library(fuzzr)
library(Colossus)
library(data.table)


df <- data.table::data.table(
  UserID = c(112, 114, 213, 214, 115, 116, 117),
  Starting_Age = c(18, 20, 18, 19, 21, 20, 18),
  Ending_Age = c(30, 45, 57, 47, 36, 60, 55),
  Cancer_Status = c(0, 0, 1, 0, 1, 0, 0),
  a = c(0, 1, 1, 0, 1, 0, 1),
  b = c(1, 1.1, 2.1, 2, 0.1, 1, 0.2),
  c = c(10, 11, 10, 11, 12, 9, 11),
  d = c(0, 0, 0, 1, 1, 1, 1),
  e = c(0, 0, 1, 0, 0, 0, 1)
)
control <- list(
  ncores = 1, lr = 0.75, maxiters = c(1, 1),
  halfmax = 1
)
formula <- Cox(Starting_Age, Ending_Age, Cancer_Status) ~
  loglinear(a, b, c, 0) + plinear(d, 0) + multiplicative()

#test_all(), test_char(), test_int(), test_dbl(), test_lgl(), test_fctr(), test_date(), test_raw(), test_df(), test_null()

#gradient_control <- "a"
#
#e <- CoxRun(df=df, model=formula, gradient_control = gradient_control, verbose = 0)
#print(e)
#
#stop()


### ------------------- Model testing ------------------------ ###
fr <- fuzz_function(CoxRun, "df", model = formula, ncores = 2, tests = test_all(), check_args = FALSE)
for (i in 1:length(fr)){
    res <- fr[[i]]$test_result$errors
    test_name <- fr[[i]]$test_name$x
    if (is.null(res)){ stop(test_name) }
}

### ------------------- Model testing ------------------------ ###
fr <- fuzz_function(CoxRun, "model", df=df, ncores = 2, tests = test_all(), check_args = FALSE)
for (i in 1:length(fr)){
    res <- fr[[i]]$test_result$errors
    test_name <- fr[[i]]$test_name$x
    if (is.null(res)){ stop(test_name) }
}

### ------------------- Numeric arguement ------------------------ ###
for (num_arg in c("a_n", "keep_constant")) {
    # Definately failing
    for (test_func in c("test_char","test_date","test_df")){
        fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        for (i in 1:length(fr)){
            res <- fr[[i]]$test_result$errors
            if (is.null(res)){
                stop(paste(num_arg,test_func, sep = " "))
            }
        }
    }
    # Could be passing
    for (test_func in c("test_lgl", "test_null")){
        fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        for (i in 1:length(fr)){
            res <- fr[[i]]$test_result$errors
            arg_in <- fr[[i]]$test_result$call$args
            arg_in <- arg_in[[num_arg]]
            if ((!is.null(res)) && (all(!is.na(arg_in)))) {
                stop(paste(num_arg,test_func, sep = " "))
            }
        }
    }
    # Remaining options to check for memory failures
    for (test_func in c("test_int", "test_dbl", "test_fctr", "test_raw")){
        fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
    }
}

for (num_arg in c("cons_mat", "cons_vec", "lr", "maxiter", "maxiters", "halfmax", "epsilon", "ll_epsilon", "deriv_epsilon", "step_max", "thres_step_max")) {
    # Definately failing
    for (test_func in c("test_char", "test_fctr","test_date","test_df", "test_null")){
        if (num_arg == "cons_vec") {
            fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, cons_mat = c(1,-1), ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        } else {
            fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        }
        for (i in 1:length(fr)){
            res <- fr[[i]]$test_result$errors
            test_name <- fr[[i]]$test_name
            test_name <- test_name[[num_arg]]
            if ((is.null(res)) && (!test_name %in% c("char_empty","fctr_empty", "df_one_col"))){
                print(paste(num_arg,test_func, sep = " "))
                print(fr[[i]]$test_result)
                print(test_name)
                stop(paste(num_arg,test_func, sep = " "))
            }
        }
    }
    # Remaining options to check for memory failures
    for (test_func in c("test_int", "test_dbl", "test_lgl", "test_raw")){
        if (num_arg == "cons_vec") {
            fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, cons_mat = c(1,-1), ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        } else {
            fr <- fuzz_function(CoxRun, num_arg, model = formula, df=df, ncores = 2, tests = do.call(test_func, args = list()), check_args = FALSE)
        }
    }
}

### ------------------- List arguement ------------------------ ###
for (cntrl_arg in c("control","gradient_control")) {
    fr <- fuzz_function(CoxRun, cntrl_arg, model = formula, df=df, ncores = 2, tests = test_all(), check_args = FALSE)
    for (i in 1:length(fr)){
        res <- fr[[i]]$test_result$errors
        test_name <- fr[[i]]$test_name
        test_name <- test_name[[cntrl_arg]]
        if ((is.null(res)) && (!test_name %in% c("df_empty", "char_empty","fctr_empty"))){
            print(fr[[i]])
            stop(paste(cntrl_arg,test_name, sep = " "))
        }
    }
}
