library(data.table)
library(parallel)
library(Colossus)

fname <- 'base_example.csv'
df <- fread(fname)

time1 <- "entry"
time2 <- "exit"
event <- "event"
names <- c("dose0","dose1")
Term_n <- c(0,0)
tform <- c("loglin","loglin")
keep_constant <- c(0,0)
a_n <- c(0,0)
modelform <- "M"
fir <- 0
der_iden <- 0
#

#a_n <- c(-0.6067, 5.019)
#model_control=list( 'basic'=TRUE, 'Log_Bound'=TRUE, 'alpha'=0.1)
#control=list("Ncores"=10,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#
#alphas <- c(0.005)
#for (alpha_i in 1:length(alphas)){
#    a_n <- c(-0.6067, 5.019)
#    model_control=list( 'basic'=TRUE, 'Log_Bound'=TRUE, 'alpha'=alphas[alpha_i], 'para_number'=0)
#    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
#}

time1 <- "entry"
time2 <- "exit"
event <- "event"
names <- c("dose0","dose1","dose0")
Term_n <- c(0,0,1)
tform <- c("loglin","loglin","lin")
keep_constant <- c(0,0,0)
#a_n <- c(0.2462, 5.020, -0.5909)
a_n <- c(0.2462, 5.020,-0.7)
modelform <- "M"
fir <- 0
der_iden <- 0
model_control=list( 'basic'=FALSE, 'Log_Bound'=FALSE)
control=list("Ncores"=18,'lr' = 0.75,'maxiters' = c(1,100),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=0.1,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=TRUE, 'ties'='breslow','double_step'=1)
#e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
#
#a <- e$LogLik
#beta <- e$beta_0[3]
#print(paste(beta, a, sep=" "))

#for (i in 0:100){
#    beta <- (0.5+0.82)/100 * i -0.82
#    a_n <- c(0.2462, 5.020,beta)
#    keep_constant <- c(0,0,1)
#    model_control=list( 'basic'=FALSE, 'Log_Bound'=FALSE)
#    control=list("Ncores"=18,'lr' = 0.75,'maxiters' = c(1,100),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)
#    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
#    a <- e$LogLik
#    print(paste(beta, a, sep=" "))
#}
for (alpha in c(0.005)){#0.75, 0.5, 1-0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)){
    a_n <- c(0.2462, 5.020,-0.599)
    model_control=list( 'basic'=FALSE, 'Log_Bound'=TRUE, 'alpha'=alpha, 'para_number'=0, 'maxstep' = 100)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
    a <- e$Parameter_Limits
    print(paste((1-alpha)*100, a[1], a[2], sep=" "))
}
