library(data.table)
library(parallel)
library(Colossus)

fname <- 'base_example.csv'
df <- fread(fname)

pyr <- "exit"   
event <- "event"
names <- c("dose0","dose1")
term_n <- c(0,1)
tform <- c("loglin","lin")
keep_constant <- c(0,0)
a_n <- c(-2.917, 0.06526)
modelform <- "M"
fir <- 0
der_iden <- 0
#
model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=FALSE, 'alpha'=0.1)
control <- list("ncores"=2,'lr' = 0.75,'maxiters' = c(10,10),'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=4, 'ties'='breslow','double_step'=1, 'guesses'=10)

alpha <- 0.005
a_n <- c(-2.917, 0.06526)
model_control <- list( 'basic'=FALSE, 'maxstep'=5, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=FALSE)
e <- RunPoissonRegression_Omnibus(df,pyr, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="rand", model_control=model_control)
a_n <- c(-2.917, 0.06526)
model_control <- list( 'basic'=FALSE, 'maxstep'=5, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=TRUE)
e <- RunPoissonRegression_Omnibus(df,pyr, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="rand", model_control=model_control)