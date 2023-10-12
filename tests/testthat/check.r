library(Colossus)
library(data.table)
library(parallel)

Term_n <- c(0,1,1,0,0)
tform <- c("loglin",'quad_slope','lin', "lin_int", "lin_slope")
keep_constant <- c(0,0,0,1,0)
a_n <- list(c(1,2,3,4,5))
names <- c("a","a","a","a","a")
Cons_Mat <- matrix(c(1,2,3,4,4,2,3,1,1,3,2,4),nrow=3,byrow=T)
Cons_Vec <- c(1,0,-1)

Cons_Mat

val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Cons_Mat, Cons_Vec)
Cons_Mat <- val$Cons_Mat

Cons_Mat





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
a_n <- 2*runif(2)-1

Constraint_Matrix
val <- Correct_Formula_Order(Term_n, tform, keep_constant, a_n, names,0, Constraint_Matrix, Constraint_const)
Constraint_Matrix
stop()




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
