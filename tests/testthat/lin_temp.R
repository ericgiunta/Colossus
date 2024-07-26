library(data.table)
library(parallel)
library(Colossus)

fname <- 'base_example.csv'
df <- fread(fname)

time1 <- "entry"
time2 <- "exit"
event <- "event"
names <- c( "dose0", "dose1", "dose0" )
term_n <- c(0,0,1)
tform <- c( "loglin", "loglin", "lin" )
keep_constant <- c(0,0,0)
#a_n <- c(0.2462, 5.020, -0.5909)
a_n <- c(0.2462, 5.020,-0.7)
modelform <- "M"
fir <- 0
der_iden <- 0
#
model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=FALSE, 'alpha'=0.1)
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(10,10), 'halfmax' = 5, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=3, 'ties'='breslow', 'double_step'=1, 'guesses'=10)

alpha <- 0.005
a_n <- c(0.2462, 5.020,-0.599)
model_control <- list( 'basic'=FALSE, 'maxstep'=5, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=FALSE)
#expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control))
a_n <- c(0.2462, 5.020,-0.599)
model_control <- list( 'basic'=FALSE, 'maxstep'=5, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=1, 'manual'=TRUE)
#expect_no_error(RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control))
alpha_list <- c(0.75, 0.5, 1-0.683, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005)
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(1,1), 'halfmax' = 2, 'epsilon' = 1e-6,  'deriv_epsilon' = 1e-6, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1, 'guesses'=10)
control <- list( "ncores"=2, 'lr' = 0.75, 'maxiters' = c(10,10), 'halfmax' = 5, 'epsilon' = 1e-4,  'deriv_epsilon' = 1e-3, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=0, 'ties'='breslow', 'double_step'=1, 'guesses'=10)

v_lower <- c(-0.643365949558998, -0.677336655540846, -0.706075250211414, -0.718165409196492, -0.753647332793819, -0.773208334303991, -0.789018704115451, -0.806061085000755, -0.816875114954096)
v_upper <- c(-0.521472203247917, -0.444964438813732, -0.327862977142017, -0.235044092073815, 2.91573713669059, 3.21014641617297, 3.48490803194128, 3.82648584413642, 4.07272009904963)
for (alpha_i in 6:6){
   alpha <- alpha_list[alpha_i]
   a_n <- c(0.2462, 5.020,-0.599)
   control <- list( "ncores"=18, 'lr' = 0.75, 'maxiters' = c(10,10), 'halfmax' = 5, 'epsilon' = 1e-4,  'deriv_epsilon' = 1e-3, 'abs_max'=1.0, 'change_all'=TRUE, 'dose_abs_max'=100.0, 'verbose'=4, 'ties'='breslow', 'double_step'=1, 'guesses'=20)
   model_control <- list( 'basic'=FALSE, 'maxstep'=100, 'log_bound'=TRUE, 'alpha'=alpha, 'para_number'=2, 'manual'=TRUE)
   e <- RunCoxRegression_Omnibus(df, time1, time2, event, names, term_n=term_n, tform=tform, keep_constant=keep_constant, a_n=a_n, modelform=modelform, fir=fir, der_iden=der_iden, control=control,strat_col="nan", model_control=model_control)
   a <- e$Parameter_Limits
   print(a[2])
#   expect_equal(a[1], v_lower[alpha_i],tolerance=1e-4)
#   expect_equal(a[2], v_upper[alpha_i],tolerance=1e-4)
}
