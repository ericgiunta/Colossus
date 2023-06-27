library(data.table)
library(Colossus)
library(parallel)

fname <- 'll_0.csv'
colTypes=c("double","double","double","integer","integer")
df <- fread(fname,nThread=min(c(detectCores(),2)),data.table=TRUE,header=TRUE,colClasses=colTypes,verbose=FALSE,fill=TRUE)
time1 <- "t0"
df$pyr <- df$t1-df$t0
pyr <- "pyr"
event <- "lung"
names <- c("dose")
Term_n <- c(0)
tform <- c("loglin")
keep_constant <- c(0)
a_n <- c(0.01)
modelform <- "M"
fir <- 0
der_iden <- 0
control=list("Ncores"=2,'lr' = 0.75,'maxiter' = 20,'halfmax' = 5,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
e <-RunPoissonRegression_STRATA(df, pyr, event, names, Term_n, tform, keep_constant, a_n, modelform, fir, der_iden, control,c("fac"))

print(e)
