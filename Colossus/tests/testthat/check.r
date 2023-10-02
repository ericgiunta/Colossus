library(data.table)
library(parallel)
library(Colossus)

fname <- 'MULTI_COV.csv'
colTypes=c("double","double","integer","integer","integer")
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
pyr <- "t1"
event <- "lung"
names <- c("a","b")
Term_n <- c(0,1)
tform <- c("loglin","loglin")
keep_constant <- c(0,0)
a_n <- c(0.01,-15)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
guesses_control=list("Iterations"=2,"guesses"=2,"lin_min"=0.001,"lin_max"=1,"loglin_min"=-1, "loglin_max"=1,"lin_method"="uniform", "loglin_method"="uniform",'strata'=FALSE,'term_initial' = c(0,1),'verbose'=TRUE)
RunPoissonRegression_Guesses_CPP(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,guesses_control,Strat_Col)





