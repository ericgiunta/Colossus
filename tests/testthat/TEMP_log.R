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

a_n <- c(-0.6067, 5.019)
model_control=list( 'basic'=TRUE, 'Log_Bound'=TRUE, 'alpha'=0.1)
control=list("Ncores"=10,'lr' = 0.75,'maxiters' = c(1,1),'halfmax' = 2,'epsilon' = 1e-6,'dbeta_max' = 0.5,'deriv_epsilon' = 1e-6, 'abs_max'=1.0,'change_all'=TRUE,'dose_abs_max'=100.0,'verbose'=FALSE, 'ties'='breslow','double_step'=1)

v_lower <- c(-0.6305960, -0.6572672, -0.6817293, -0.6929630, -0.7300938, -0.7537744, -0.7749381, -0.8001031, -0.8175117)
v_upper <- c(-0.5828725, -0.5562505, -0.5318645, -0.5206756, -0.4837373, -0.4602148, -0.4392159, -0.4142752, -0.3970399)

alphas <- c(0.75, 0.5, 1-0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
for (alpha_i in 1:length(alphas)){
    a_n <- c(-0.6067, 5.019)
    model_control=list( 'basic'=TRUE, 'Log_Bound'=TRUE, 'alpha'=alphas[alpha_i], 'para_number'=0)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
    a <- e$Parameter_Limits
    print(c(a[1]-v_lower[alpha_i], a[2]-v_upper[alpha_i]))
}

v_lower <- c(4.981497, 4.939337, 4.900838, 4.883211, 4.825191, 4.788380, 4.755608, 4.716794, 4.690041)
v_upper <- c(5.057414, 5.100032, 5.139239, 5.157283, 5.217094, 5.255376, 5.289680, 5.330581, 5.358945)
for (alpha_i in 1:length(alphas)){
    a_n <- c(-0.6067, 5.019)
    model_control=list( 'basic'=TRUE, 'Log_Bound'=TRUE, 'alpha'=alphas[alpha_i], 'para_number'=1)
    e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, Term_n=Term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,Strat_Col="nan", model_control=model_control)
    a <- e$Parameter_Limits
    print(c(a[1]-v_lower[alpha_i], a[2]-v_upper[alpha_i]))
}
