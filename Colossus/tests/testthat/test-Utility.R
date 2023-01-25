
test_that("No truncation columns", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_equal(Check_Trunc(df,c("time0","time1"))$ce, c("time0","time1"))
})
test_that("Right truncation columns", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_equal(Check_Trunc(df,c("%trunc%","time1"))$ce, c("right_trunc","time1"))
})
test_that("Left truncation columns", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_equal(Check_Trunc(df,c("time0","%trunc%"))$ce, c("time0","left_trunc"))
})
test_that("Truncation no column error", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_error(Check_Trunc(df,c()))
})
test_that("Truncation left column not in df error", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_error(Check_Trunc(df,c("timebad","%trunc%")))
})
test_that("Truncation right column not in df error", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_error(Check_Trunc(df,c("%trunc%","timebad")))
})
test_that("Truncation both sides", {
    df <- data.table("time0"=c(0,1,2,3,4,5,6),"time1"=c(1,2,3,4,5,6,7),"dummy"=c(0,0,1,1,0,1,0))
    expect_error(Check_Trunc(df,c("%trunc%","%trunc%")))
})


test_that("No dupe columns", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","d"),c(0,0,0,0),FALSE), c("a","b","c","d"))
})
test_that("No columns", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_equal(Check_Dupe_Columns(df,c(),c(),FALSE), c())
})
test_that("One column with varying", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_equal(Check_Dupe_Columns(df,c("a"),c(0),FALSE), c("a"))
})
test_that("One column with constant", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_equal(Check_Dupe_Columns(df,c("c"),c(0),FALSE), c("c"))
})
test_that("One duplicate column", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=a)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","d","e"),c(0,0,0,0,0),FALSE), c("a","b","c","d"))
})
test_that("One duplicate column, different term", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=a)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","d","e"),c(0,0,0,1,1),FALSE), c("a","b","c","d","e"))
})
test_that("Multiple duplicate columns", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d,"e"=a,"f"=b)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","e","f"),c(0,0,0,0,0),FALSE), c("a","b","c"))
})
test_that("All duplicate columns, different terms", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=a,"c"=a,"d"=a,"e"=a,"f"=a)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","e","f"),c(0,1,2,3,4),FALSE), c("a","b","c","e","f"))
})
test_that("Repeated duplicate columns", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=a,"e"=a,"f"=a)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c","d","f"),c(0,0,0,0,0),FALSE), c("a","b","c"))
})
test_that("All but one duplicate column with varying", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=a,"c"=a)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c"),c(0,0,0),FALSE), c("a"))
})
test_that("All but one duplicate column with constant", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=c,"b"=c,"c"=c)
    expect_equal(Check_Dupe_Columns(df,c("a","b","c"),c(0,0,0),FALSE), c())
})
test_that("Duplicate with column not in df error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=c,"b"=c,"c"=c)
    expect_error(Check_Dupe_Columns(df,c("a","b","c","e"),c(0,0,0,0),FALSE))
})

test_that("Improve Ratio test", {
    a <- list("LogLik"=-400)
    b <- list("LogLik"=-350)
    expect_equal(Likelihood_Ratio_Test(b,a), 100)
})
test_that("Worse Ratio test", {
    a <- list("LogLik"=-300)
    b <- list("LogLik"=-350)
    expect_equal(Likelihood_Ratio_Test(b,a), -100)
})
test_that("Same Ratio test", {
    a <- list("LogLik"=-300)
    b <- list("LogLik"=-300)
    expect_equal(Likelihood_Ratio_Test(a,b), 0)
})
test_that("No Data Ratio test", {
    a <- list("baditem"=-300)
    b <- list("LogLik"=-300)
    expect_error(Likelihood_Ratio_Test(a,b))
})

test_that("Iteract no dupes", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?+?b","a?*?b")
    new_names <- c("","")
    expect_equal(interact_them(df,interactions,new_names,FALSE)$cols, c("a+b","a*b"))
})
test_that("Iteract no dupes with rename", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?+?b","a?*?b")
    new_names <- c("","formtemp")
    expect_equal(interact_them(df,interactions,new_names,FALSE)$cols, c("a+b","formtemp"))
})
test_that("Iteract with direct dupes", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?+?b","a?*?b","a?+?b","a?+?a")
    new_names <- c("","","","")
    expect_equal(interact_them(df,interactions,new_names,FALSE)$cols, c("a*b","a+b","a+a"))
})
test_that("Iteract with reverse dupes", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?+?b","a?*?b","b?+?a","a?+?a")
    new_names <- c("","","","")
    expect_equal(interact_them(df,interactions,new_names,FALSE)$cols, c("a*b","b+a","a+a"))
})
test_that("Iteract formula long error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?+?b?+c","a?*?b")
    new_names <- c("","")
    expect_error(interact_them(df,interactions,new_names,FALSE))
})
test_that("Iteract formula operation error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=c,"b"=c,"c"=c)
    interactions <- c("a?++?b","a?*?b")
    new_names <- c("","")
    expect_error(interact_them(df,interactions,new_names,FALSE))
})

#######################################
## FACTORING
#######################################

test_that("Factorize factor", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("c")
    expect_equal(factorize(df,col_list)$cols, c("c_1"))
})
test_that("Factorize discrete", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("a")
    expect_equal(factorize(df,col_list)$cols, c("a_0","a_1","a_2","a_3","a_4","a_5","a_6"))
})
test_that("Factorize missing", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("d")
    expect_error(factorize(df,col_list))
})

test_that("Factorize parallel", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,0,2,0,1,2,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("c")
    expect_no_error(factorize_parallel(df,col_list))
})
test_that("Factorize factor parallel", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("c")
    expect_equal(factorize_parallel(df,col_list)$cols, c("c_1"))
})
test_that("Factorize discrete parallel", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("a")
    expect_equal(factorize_parallel(df,col_list)$cols, c("a_0","a_1","a_2","a_3","a_4","a_5","a_6"))
})
test_that("Factorize missing parallel", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(0,0,0,0,0,0,0)
    df <- data.table("a"=a,"b"=b,"c"=c)
    col_list <- c("d")
    expect_error(factorize_parallel(df,col_list))
})


#######################################
## Time Dependent Cov gens
#######################################

test_that("Gen_time_dep time error", {
    a <- c(20,20,5,10,15)
    b <- c(1,2,1,1,2)
    c <- c(0,0,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    #
    time1="%trunc%"
    time2="a_bad"
    event="c"
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    grt_f <- function(df,time_col){
        return ((df[,"b"] * df[,get(time_col)])[[1]])
    }
    func_form <- c("lin")
    
    #
    expect_error(gen_time_dep(df,time1,time2,event,TRUE,0.01,c("grt"),c(),c(grt_f),paste("test","_new.csv",sep=""), func_form))
})
test_that("Gen_time_dep event error", {
    a <- c(20,20,5,10,15)
    b <- c(1,2,1,1,2)
    c <- c(0,0,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    #
    time1="%trunc%"
    time2="a"
    event="c_bad"
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    grt_f <- function(df,time_col){
        return ((df[,"b"] * df[,get(time_col)])[[1]])
    }
    func_form <- c("lin")
    
    #
    expect_error(gen_time_dep(df,time1,time2,event,TRUE,0.01,c("grt"),c(),c(grt_f),paste("test","_new.csv",sep=""), func_form))
})
test_that("Gen_time_dep function error", {
    a <- c(20,20,5,10,15)
    b <- c(1,2,1,1,2)
    c <- c(0,0,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    #
    time1="%trunc%"
    time2="a"
    event="c_bad"
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    grt_f <- function(df,time_col){
        stop()
        return ((df[,"b"] * df[,get(time_col)])[[1]])
    }
    func_form <- c("lin")
    
    #
    expect_error(gen_time_dep(df,time1,time2,event,TRUE,0.01,c("grt"),c(),c(grt_f),paste("test","_new.csv",sep=""), func_form))
})

test_that("Gen_time_dep no error", {
    a <- c(20,20,5,10,15)
    b <- c(1,2,1,1,2)
    c <- c(0,0,1,1,1)
    df <- data.table("a"=a,"b"=b,"c"=c)
    #
    time1="%trunc%"
    time2="a"
    event="c"
    control <- list('lr' = 0.75,'maxiter' = -1,'halfmax' = 5,'epsilon' = 1e-9,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-9, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
    grt_f <- function(df,time_col){
        return ((df[,"b"] * df[,get(time_col)])[[1]])
    }
    func_form <- c("lin")
    
    #
    expect_no_error(gen_time_dep(df,time1,time2,event,TRUE,0.01,c("grt"),c(),c(grt_f),paste("test","_new.csv",sep=""), func_form))
})

test_that("linked quad negative slope error", {
    tforms <- list("first"="quad")
    paras  <- list("first"=c(-0.1,10))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked quad string slope error", {
    tforms <- list("first"="quad")
    paras  <- list("first"=c("a",10))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked quad string threshold error", {
    tforms <- list("first"="quad")
    paras  <- list("first"=c(0.1,"a"))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked quad no error", {
    tforms <- list("first"="quad")
    paras  <- list("first"=c(0.1,10))
    expect_no_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked exp negative slope error", {
    tforms <- list("first"="exp")
    paras  <- list("first"=c(-0.1,10,5))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked exp string slope error", {
    tforms <- list("first"="exp")
    paras  <- list("first"=c("a",10,5))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked exp string threshold error", {
    tforms <- list("first"="exp")
    paras  <- list("first"=c(0.1,"a",5))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked exp string exp slope error", {
    tforms <- list("first"="exp")
    paras  <- list("first"=c(0.1,10,"a"))
    expect_error(Linked_Dose_Formula(tforms,paras,FALSE))
})
test_that("linked exp no error", {
    tforms <- list("first"="exp")
    paras  <- list("first"=c(0.1,10,5))
    expect_no_error(Linked_Dose_Formula(tforms,paras,FALSE))
})

test_that("linked exp parameter low goal error", {
    y=10
    a0=1
    a_goal=5
    expect_error(Linked_Lin_Exp_Para(y,a0,a_goal,FALSE))
})
test_that("linked exp parameter negative slope error", {
    y=10
    a0=-0.1
    a_goal=5
    expect_error(Linked_Lin_Exp_Para(y,a0,a_goal,FALSE))
})
test_that("linked exp parameter no error", {
    y=10
    a0=0.1
    a_goal=5
    expect_no_error(Linked_Lin_Exp_Para(y,a0,a_goal,FALSE))
})

test_that("Missing Value missing column error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_error(Replace_Missing(df,c("a","e"),0.0))
})
test_that("Missing Value NA replacement error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_error(Replace_Missing(df,c("a","b","c","d"),NA))
})
test_that("Missing Value no error", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_no_error(Replace_Missing(df,c("a","b","c","d"),0.0))
})
test_that("Missing Value checked replaced 0", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(NA,0,0,1,0,0,1)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    #
    df0 <- Replace_Missing(df,c("a","b"),0.0)
    expect_equal(c(sum(df0$a),sum(df0$b)),c(sum(df$a),2))
})
test_that("Missing Value checked replaced 1", {
    a <- c(0,1,2,3,4,5,6)
    b <- c(NA,0,0,1,0,0,1)
    c <- c(1,1,1,1,1,1,1)
    d <- c(3,4,5,6,7,8,9)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    #
    df0 <- Replace_Missing(df,c("a","b"),1.0)
    expect_equal(c(sum(df0$a),sum(df0$b)),c(sum(df$a),3))
})

test_that("Build Cluster, no error", {
    a <- c(1,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(10,10.5,11,12,15,20,22)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_no_error(Check_Cluster(df,c("a","b"),"d"))
})
test_that("Build Cluster, Real Example", {
    a <- c(2.89590128732351,2.89561870735048,2.88729426293386,2.88632109713572,2.88682401866894,2.88814871952403,2.86044875402007,2.83712327958272,2.81656778565829,2.71489702667994,2.71047710238058,2.70809002989701,2.33321629448343,2.37336593993027,1.44369764007809)
    b <- c(5.76915890997471,5.77293068875401,5.76924070126138,5.76813274714431,5.76940952447491,5.78663890952283,5.81092132443917,5.81053871936002,5.83965169734053,5.98042757357243,6.00265348503193,6.07703757564256,6.6155048455783,6.44020280315287,7.83732630045058)
    c <- c(0.960955880367663,0.407925244234502,-6.37203860087871,0.418717290348606,0.625047282009153,0.531352484038798,0.615650905764196,1.13425252520829,0.885062725279015,0.577218915238278,0.227952190984506,0.350085715007503,0.465717396184068,0.276913823611103,0.491846674636864)
    d <- c(742.267573600719,742.325896394396,742.488463168289,742.52882105719,742.538183765492,742.675989286731,743.76102685265,745.162337013951,747.19634619851,760.337874781641,761.940876952424,766.359832891647,823.84782476707,867.630547786468,2242.53849745672)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_no_error(Check_Cluster(df,c("a","b","c"),"d"))
})
test_that("Build Cluster, Real Example Check a_n", {
    a <- c(2.89590128732351,2.89561870735048,2.88729426293386,2.88632109713572,2.88682401866894,2.88814871952403,2.86044875402007,2.83712327958272,2.81656778565829,2.71489702667994,2.71047710238058,2.70809002989701,2.33321629448343,2.37336593993027,1.44369764007809)
    b <- c(5.76915890997471,5.77293068875401,5.76924070126138,5.76813274714431,5.76940952447491,5.78663890952283,5.81092132443917,5.81053871936002,5.83965169734053,5.98042757357243,6.00265348503193,6.07703757564256,6.6155048455783,6.44020280315287,7.83732630045058)
    c <- c(0.960955880367663,0.407925244234502,-6.37203860087871,0.418717290348606,0.625047282009153,0.531352484038798,0.615650905764196,1.13425252520829,0.885062725279015,0.577218915238278,0.227952190984506,0.350085715007503,0.465717396184068,0.276913823611103,0.491846674636864)
    d <- c(742.267573600719,742.325896394396,742.488463168289,742.52882105719,742.538183765492,742.675989286731,743.76102685265,745.162337013951,747.19634619851,760.337874781641,761.940876952424,766.359832891647,823.84782476707,867.630547786468,2242.53849745672)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    #
    e <- Check_Cluster(df,c("a","b","c"),"d")
    expect_equal(e$a_n,c(2.89590128732351,5.76915890997471,0))
})
test_that("Build Cluster, Real Example Check guess_constant", {
    a <- c(2.89590128732351,2.89561870735048,2.88729426293386,2.88632109713572,2.88682401866894,2.88814871952403,2.86044875402007,2.83712327958272,2.81656778565829,2.71489702667994,2.71047710238058,2.70809002989701,2.33321629448343,2.37336593993027,1.44369764007809)
    b <- c(5.76915890997471,5.77293068875401,5.76924070126138,5.76813274714431,5.76940952447491,5.78663890952283,5.81092132443917,5.81053871936002,5.83965169734053,5.98042757357243,6.00265348503193,6.07703757564256,6.6155048455783,6.44020280315287,7.83732630045058)
    c <- c(0.960955880367663,0.407925244234502,-6.37203860087871,0.418717290348606,0.625047282009153,0.531352484038798,0.615650905764196,1.13425252520829,0.885062725279015,0.577218915238278,0.227952190984506,0.350085715007503,0.465717396184068,0.276913823611103,0.491846674636864)
    d <- c(742.267573600719,742.325896394396,742.488463168289,742.52882105719,742.538183765492,742.675989286731,743.76102685265,745.162337013951,747.19634619851,760.337874781641,761.940876952424,766.359832891647,823.84782476707,867.630547786468,2242.53849745672)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    #
    e <- Check_Cluster(df,c("a","b","c"),"d")
    expect_equal(e$guess_constant,c(1,1,0))
})
test_that("Build Cluster, missing dev", {
    a <- c(1,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(10,10.5,11,12,15,20,22)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_error(Check_Cluster(df,c("a","b"),"e"))
})
test_that("Build Cluster, missing name", {
    a <- c(1,1,2,3,4,5,6)
    b <- c(1,2,3,4,5,6,7)
    c <- c(1,1,1,1,1,1,1)
    d <- c(10,10.5,11,12,15,20,22)
    df <- data.table("a"=a,"b"=b,"c"=c,"d"=d)
    expect_error(Check_Cluster(df,c("a","e"),"d"))
})



